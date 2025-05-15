import pysam
import argparse
import os

def split_bam_by_tag(input_bam_path, tag, output_dir):
    """
    Splits a BAM file into multiple files based on the specified tag.

    Parameters:
        input_bam_path (str): Path to the input BAM file.
        tag (str): The tag to split the BAM file by (e.g., BC, CB, RG).
        output_dir (str): Directory to save the split BAM files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Open the input BAM file
    input_bam = pysam.AlignmentFile(input_bam_path, "rb")

    # Create a dictionary to hold BAM writers for each tag value
    writers = {}

    for read in input_bam:
        # Skip reads without the specified tag
        if not read.has_tag(tag):
            continue

        # Get the value of the tag
        tag_value = read.get_tag(tag)

        # Create a BAM writer for this tag value if not already created
        if tag_value not in writers:
            output_bam_path = os.path.join(output_dir, f"{tag}_{tag_value}.bam")
            writers[tag_value] = pysam.AlignmentFile(output_bam_path, "wb", header=input_bam.header)

        # Write the read to the appropriate BAM file
        writers[tag_value].write(read)

    # Close all file handles
    input_bam.close()
    for writer in writers.values():
        writer.close()

    print(f"BAM file split by tag '{tag}' completed. Output files are in: {output_dir}")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Split a BAM file based on a user-defined tag.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input BAM file.")
    parser.add_argument("-t", "--tag", required=True, help="Tag to split the BAM file by (e.g., BC, CB, RG).")
    parser.add_argument("-o", "--output_dir", default="split_bam_output", help="Directory to save the split BAM files (default: split_bam_output).")

    args = parser.parse_args()

    # Call the function with parsed arguments
    split_bam_by_tag(args.input, args.tag, args.output_dir)
