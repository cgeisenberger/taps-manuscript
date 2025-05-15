#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict, Counter

def parse_args():
    """Handle command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute conversion efficiency of CG dinucleotides "
                    "with consensus-based deduplication."
    )
    parser.add_argument(
        "-r", "--reference",
        required=True,
        help="Path to the reference genome FASTA file (indexed)."
    )
    parser.add_argument(
        "-b", "--bam",
        required=True,
        help="Path to the mapped, sorted, and indexed BAM file."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output file."
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=0,
        help="Minimum mapping quality to keep a read (default: 0)."
    )
    parser.add_argument(
        "--max_insert_size",
        type=int,
        default=None,
        help="Maximum allowed insert size (optional)."
    )
    return parser.parse_args()

def find_cg_positions(reference_fasta):
    """
    Scan the reference FASTA and return a dict[chrom] -> list of CG positions.
    Each position is the index of the 'C' in 'CG'.
    """
    ref = pysam.FastaFile(reference_fasta)
    cg_dict = defaultdict(list)

    for chrom in ref.references:
        seq = ref.fetch(chrom)
        idx = 0
        while True:
            idx = seq.find("CG", idx)
            if idx == -1:
                break
            cg_dict[chrom].append(idx)
            idx += 1  # move past this 'C'
    
    # Sort positions for each chrom
    for chrom in cg_dict:
        cg_dict[chrom].sort()
    ref.close()
    return cg_dict

def get_base_context(reference_fasta, chrom, position):
    """
    Return the +/-1 base context around the CG site for output.
    position is the index of 'C'.
    For safety, handle edges by fetching at most from (position-1) to (position+2).
    """
    ref = pysam.FastaFile(reference_fasta)
    start = max(0, position - 1)
    end = position + 3  # end is exclusive
    subseq = ref.fetch(chrom, start, end)
    ref.close()
    return subseq

def majority_vote_base(bases):
    """
    Given a list of single-character bases, return the majority base.
    If there's a tie, return 'N'.
    """
    # Count each base
    counts = Counter(bases)
    # Find the top count
    max_count = max(counts.values())  # e.g. 5 if base 'T' appears 5 times
    # Get all bases with that count
    major_bases = [b for b, c in counts.items() if c == max_count]
    if len(major_bases) == 1:
        return major_bases[0]
    else:
        # tie
        return 'N'

def main():
    args = parse_args()
    ref_fasta = args.reference
    bam_file = args.bam
    output_file = args.output
    min_mapq = args.min_mapq
    max_insert_size = args.max_insert_size

    # 1) Identify all CG positions from the reference
    cg_dict = find_cg_positions(ref_fasta)  # dict: chrom -> list of positions

    # 2) We'll store final counts in data[(chrom, cg_pos)] = [converted_count, unconverted_count].
    data = defaultdict(lambda: [0, 0])

    # 3) Open BAM for reading
    bam = pysam.AlignmentFile(bam_file, "rb")

    # 4) Process each chromosome in the reference
    #    We'll fetch all reads from that chromosome (bam.fetch(chrom)),
    #    then group them by (r1_start, CB, RX).
    #    For each group, we collect the 2-bp calls at each CG position they cover,
    #    and then produce a consensus for that group at that position.
    #
    #    We then do a single converted/unconverted increment for that group at that CG.
    
    for chrom in cg_dict.keys():
        cg_positions = cg_dict[chrom]
        if not cg_positions:
            continue  # no CGs found on this chrom

        # Data structure for reads grouped by (r1_start, CB, RX)
        # For each group, we track coverage of CG positions:
        #   read_groups[(r1_start, cb, rx)][cg_pos] = list of (base1, base2) from reads that cover it
        read_groups = defaultdict(lambda: defaultdict(list))

        # ---- STEP A: FETCH & FILTER READS ----
        for read in bam.fetch(chrom):
            # 1) Filter out reads that are unmapped or their mate is unmapped
            #    or that are not properly paired.
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if not read.is_paired or not read.is_proper_pair:
                continue

            # 2) Filter by mapping quality
            if read.mapping_quality < min_mapq:
                continue

            # 3) Filter by insert size if requested
            if max_insert_size is not None:
                if abs(read.template_length) > max_insert_size:
                    continue

            # If the read does not have the required tags, treat them as empty
            cb = read.get_tag("CB") if read.has_tag("CB") else ""
            rx = read.get_tag("RX") if read.has_tag("RX") else ""

            # Deduplication group key: (read1_position, CB, RX)
            # read1_position is read.reference_start if read.is_read1
            # otherwise read.next_reference_start if read.is_read2
            if read.is_read1:
                r1_start = read.reference_start
            else:
                r1_start = read.next_reference_start

            # This group will store the 2-bp read calls per CG site
            group_key = (r1_start, cb, rx)

            # 4) Identify which CGs this read covers
            #    We'll use read.get_aligned_pairs(matches_only=True) to see
            #    the alignment from read base to reference base.
            aligned_pairs = read.get_aligned_pairs(matches_only=False)

            # Safely build a lookup
            ref_to_query = {}
            for qpos, rpos in aligned_pairs:
                if rpos is not None and qpos is not None:
                    if 0 <= qpos < read.query_length:
                        ref_to_query[rpos] = qpos
            # We'll iterate over the CG positions in this chromosome,
            # checking if the read covers cg_pos and cg_pos+1.  
            # (Alternatively, we could iterate over the aligned_pairs directly,
            #  but we already have the cg_positions sorted.  For large data, you might
            #  want a more efficient approach, e.g., searching only relevant sub-range.)
            
            # A small optimization: skip all CGs outside read.reference_start..read.reference_end
            read_start = read.reference_start
            read_end = read.reference_end  # 1-based exclusive in pysam

            # We'll do a quick bounding approach: 
            # find the first CG pos >= read_start, then proceed while cg_pos < read_end
            # This can be done via a simple while or binary search. For clarity, do a while:
            # We'll do a manual index approach for cg_positions
            # (This is optional, but speeds up for large data.)
            
            # find the starting index for CG >= read_start
            # a naive approach is to move an index forward until cg_positions[idx] >= read_start
            # to keep the code simpler, just do a for with a break condition.
            for cg_pos in cg_positions:
                if cg_pos < read_start:
                    continue
                if cg_pos >= read_end - 1:
                    break  # no need to check further if the C is beyond the read coverage

                # Check that both cg_pos and cg_pos+1 are in ref_to_query
                if cg_pos in ref_to_query and (cg_pos+1) in ref_to_query:
                    qpos_c = ref_to_query[cg_pos]
                    qpos_g = ref_to_query[cg_pos+1]
                    base_c = read.query_sequence[qpos_c]
                    base_g = read.query_sequence[qpos_g]
                    read_groups[group_key][cg_pos].append((base_c, base_g))

        # ---- STEP B: BUILD CONSENSUS PER GROUP & UPDATE COUNTS ----
        # For each group, we go through each CG position it covers and build a consensus call.
        # Then we increment converted/unconverted counts in data.
        for group_key, cg_calls_dict in read_groups.items():
            for cg_pos, base_pairs_list in cg_calls_dict.items():
                # base_pairs_list is e.g. [('C','G'), ('T','G'), ('C','G'), ...]
                # We'll find the consensus for position cg_pos (the 'C') and cg_pos+1 (the 'G')
                
                if not base_pairs_list:
                    continue  # no coverage?

                # Separate out the first base calls and second base calls
                base1_list = [bp[0] for bp in base_pairs_list]
                base2_list = [bp[1] for bp in base_pairs_list]

                consensus_base1 = majority_vote_base(base1_list)
                consensus_base2 = majority_vote_base(base2_list)
                consensus_pair = consensus_base1 + consensus_base2

                # Classify as converted / unconverted
                # "TG" or "CA" => converted, "CG" => unconverted, else skip
                if consensus_pair in ["TG", "CA"]:
                    data[(chrom, cg_pos)][0] += 1  # converted
                elif consensus_pair == "CG":
                    data[(chrom, cg_pos)][1] += 1  # unconverted
                else:
                    # e.g. tie => "NG" or "TN" or partial mismatch => skip
                    pass

    # Close the BAM file
    bam.close()

    # ---- STEP C: WRITE OUTPUT ----
    # We omit sites with zero coverage (converted+unconverted=0).
    with open(output_file, "w") as out:
        header = ["chromosome", "coordinate", "num_converted", "num_unconverted", "context"]
        out.write("\t".join(header) + "\n")
        
        for (chrom, pos), (conv_count, unconv_count) in sorted(data.items()):
            if (conv_count + unconv_count) == 0:
                continue  # omit zero coverage
            context_seq = get_base_context(ref_fasta, chrom, pos)
            out.write(
                f"{chrom}\t{pos}\t{conv_count}\t{unconv_count}\t{context_seq}\n"
            )

if __name__ == "__main__":
    main()

