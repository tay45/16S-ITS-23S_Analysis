#!/usr/bin/env python3

import argparse
import subprocess
import os
import glob
import logging
import sys
import shutil
import csv
import re
import pandas as pd

def setup_logging(log_file):
    """
    Sets up logging to both file and stdout with the specified log level.
    """
    logging.basicConfig(
        level=logging.INFO,  # Set to INFO to reduce verbosity
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

def run_command(cmd, description):
    """
    Runs a shell command and logs its output.
    """
    logging.info(f"Starting: {description}")
    try:
        result = subprocess.run(
            cmd, shell=True, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        if result.stdout.strip():
            logging.info(result.stdout.strip())
        if result.stderr.strip():
            logging.warning(result.stderr.strip())
        logging.info(f"Completed: {description}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed: {description}")
        logging.error(e.stderr.strip())
        raise

def skera_split(input_bam, adapters_fasta, output_bam):
    """
    Runs the Skera Split step.
    """
    cmd = f"skera split {input_bam} {adapters_fasta} {output_bam}"
    run_command(cmd, "Skera Split")

def lima_demultiplex(input_bam, barcodes_fasta, output_prefix, barcode_type):
    """
    Runs the Lima Demultiplexing step.
    """
    option = "--same" if barcode_type.lower() == "symmetric" else "--different"
    output_bam = f"{output_prefix}.bam"
    cmd = f"lima {input_bam} {barcodes_fasta} {output_bam} {option} --split-bam-named --min-score 26"
    run_command(cmd, "Lima Demultiplexing")

def extract_primer_directions(bam_name):
    """
    Extracts the barcode from the BAM filename.
    Expected pattern: *.Kinnex16S_Fwd_**--Kinnex16S_Rev_**.bam
    Returns:
        barcode (str): Extracted barcode
        status (str): "Valid" or specific issue
    """
    # Updated regex to capture the barcode regardless of prefixes
    pattern = r'(Kinnex16S_Fwd_\d+--Kinnex16S_Rev_\d+)\.bam$'
    match = re.search(pattern, bam_name)
    if not match:
        logging.debug(f"Filename does not match expected pattern: {bam_name}")
        return "", "Incorrect filename format"
    
    barcode = match.group(1)
    logging.debug(f"Extracted barcode: {barcode}")
    
    # Validate primer directions within the barcode
    if "Fwd" in barcode and "Rev" in barcode:
        return barcode, "Valid"
    else:
        return barcode, "Incorrect primer directions"

def move_invalid_files(bam_files, invalid_dir, report_file):
    """
    Moves invalid BAM files to the invalid directory and logs issues.
    Returns:
        valid_bams (list): List of valid BAM file paths.
    """
    os.makedirs(invalid_dir, exist_ok=True)
    valid_bams = []
    
    with open(report_file, 'w', newline='') as csvfile:
        fieldnames = ['BAM File', 'Issue']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
    
        for bam in bam_files:
            bam_name = os.path.basename(bam)
            barcode, status = extract_primer_directions(bam_name)
            size_zero = os.path.getsize(bam) == 0
            issues = []
            if status != "Valid":
                issues.append(status)
            if size_zero:
                issues.append("Zero file size")
            
            if issues:
                writer.writerow({'BAM File': bam, 'Issue': '; '.join(issues)})
                target_subdir = "_".join([issue.replace(" ", "_") for issue in issues])
                target_dir = os.path.join(invalid_dir, target_subdir)
                os.makedirs(target_dir, exist_ok=True)
                
                # Move BAM file and related files
                for file_ext in ['', '.pbi']:
                    file_pattern = f"{bam}{file_ext}"
                    if os.path.exists(file_pattern):
                        try:
                            shutil.move(file_pattern, target_dir)
                            logging.info(f"Moved {file_pattern} to {target_dir}")
                        except shutil.Error as move_err:
                            logging.error(f"Error moving {file_pattern} to {target_dir}: {move_err}")
                
                # Handle XML files separately
                xml_pattern = f"{bam.replace('.bam', '')}*.xml"
                for file in glob.glob(xml_pattern):
                    if os.path.exists(file):
                        try:
                            shutil.move(file, target_dir)
                            logging.info(f"Moved {file} to {target_dir}")
                        except shutil.Error as move_err:
                            logging.error(f"Error moving {file} to {target_dir}: {move_err}")
            else:
                logging.info(f"Valid BAM file: {bam}")
                valid_bams.append(bam)
    return valid_bams

def convert_bam_to_fastx(bam_files, output_dir, convert_types=['fastq', 'fasta'], compression_level=1, uncompressed=False):
    """
    Converts BAM files to FASTA and/or FASTQ using bam2fasta and bam2fastq.
    """
    for bam in bam_files:
        bam_basename = os.path.splitext(os.path.basename(bam))[0]
        output_prefix = os.path.join(output_dir, bam_basename)
        if 'fastq' in convert_types:
            cmd = f"bam2fastq -o \"{output_prefix}\" \"{bam}\""
            if uncompressed:
                cmd += " -u"
            else:
                cmd += f" -c {compression_level}"
            run_command(cmd, f"Converting {bam} to FASTQ")
        
        if 'fasta' in convert_types:
            cmd = f"bam2fasta -o \"{output_prefix}\" \"{bam}\""
            if uncompressed:
                cmd += " -u"
            else:
                cmd += f" -c {compression_level}"
            run_command(cmd, f"Converting {bam} to FASTA")

def generate_fasta_index(fasta_files):
    """
    Generates FASTA index files using samtools faidx.
    """
    for fasta in fasta_files:
        if not os.path.exists(fasta):
            logging.error(f"FASTA file not found: {fasta}")
            continue
        if os.path.getsize(fasta) == 0:
            logging.warning(f"FASTA file is empty, skipping indexing: {fasta}")
            continue
        cmd = f"samtools faidx \"{fasta}\""
        run_command(cmd, f"Generating FASTA index for {fasta}")

def ensure_tab_delimited(csv_path, converted_csv_path):
    """
    Ensures the CSV file is tab-delimited and correctly formatted. Converts it if necessary.
    Args:
        csv_path (str): Path to the original CSV.
        converted_csv_path (str): Path to save the converted tab-delimited CSV.
    Returns:
        None
    """
    try:
        # Open the original CSV, strip leading/trailing quotes from each line, and write to a cleaned temp file
        with open(csv_path, 'r') as f_in, open(converted_csv_path, 'w') as f_out:
            for line in f_in:
                # Strip leading and trailing whitespace and quotes
                line = line.strip().strip('"')
                f_out.write(line + '\n')

        # Now read the cleaned CSV with pandas, assuming it's tab-delimited
        df = pd.read_csv(converted_csv_path, delimiter='\t')

        # Check if required columns exist (case-insensitive)
        columns_lower = [col.lower() for col in df.columns]
        if 'barcode' not in columns_lower or 'sample name' not in columns_lower:
            logging.error("CSV must contain 'Barcode' and 'Sample Name' columns.")
            sys.exit(1)
        
        # Rename columns to standard names regardless of original case
        df.columns = [col.strip().lower() for col in df.columns]
        df.rename(columns={'barcode': 'Barcode', 'sample name': 'Sample Name'}, inplace=True)
        
        # Save the cleaned DataFrame as tab-delimited without quotes
        df.to_csv(converted_csv_path, sep='\t', index=False, quoting=csv.QUOTE_NONE, escapechar='\\')
        logging.info(f"Tab-delimited CSV saved as: {converted_csv_path}")
    except pd.errors.ParserError:
        logging.error("Failed to parse the CSV file. Please ensure it is correctly formatted.")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error processing CSV file: {e}")
        sys.exit(1)

def combine_fasta(combined_fasta, group_file, barcode_to_sample, valid_fasta_files, output_dir):
    """
    Combines individual FASTA files into a single combined FASTA file and creates a group file.
    Renames sequence headers to include the corresponding Bio Sample name while retaining original sequence ID.
    """
    logging.info("Combining FASTA files and generating group file...")
    
    if not valid_fasta_files:
        logging.error("No valid FASTA files to combine.")
        sys.exit(1)
    logging.info(f"Found {len(valid_fasta_files)} valid FASTA files to combine.")
    
    try:
        with open(combined_fasta, 'w') as outfile, open(group_file, 'w') as grp:
            grp.write("sequenceID\tgroup\n")  # Header for group file
            for fasta_file in valid_fasta_files:
                # Extract barcode from fasta filename
                # Assuming the filename contains the barcode as *.Kinnex16S_Fwd_01--Kinnex16S_Rev_13.fasta
                barcode_match = re.search(r'(Kinnex16S_Fwd_\d+--Kinnex16S_Rev_\d+)', fasta_file)
                if barcode_match:
                    barcode = barcode_match.group(1)
                    sample_name = barcode_to_sample.get(barcode, "Unknown")
                else:
                    logging.warning(f"Could not extract barcode from filename: {fasta_file}")
                    sample_name = "Unknown"

                if sample_name == "Unknown":
                    logging.warning(f"Sample key '{barcode}' not found in barcode_to_sample mapping. Assigning as 'Unknown'.")

                fasta_path = os.path.join(output_dir, fasta_file)
                logging.info(f"Processing FASTA file: {fasta_path} | Sample Name: {sample_name}")
                try:
                    with open(fasta_path, 'r') as infile:
                        seq_counter = 1
                        for line in infile:
                            if line.startswith('>'):
                                original_seq_id = line[1:].strip()
                                # Modify the sequence ID by appending Bio Sample name as postfix
                                # Example: >original_seq_id_Bio_Sample_12_seq1
                                # Replace spaces with underscores in sample name
                                sample_name_clean = sample_name.replace(' ', '_')
                                unique_seq_id = f"{original_seq_id}_{sample_name_clean}_seq{seq_counter}"
                                grp.write(f"{unique_seq_id}\t{sample_name}\n")
                                outfile.write(f">{unique_seq_id}\n")
                                seq_counter += 1
                            else:
                                outfile.write(line)
                except Exception as e:
                    logging.error(f"Error processing {fasta_path}: {e}")
                    sys.exit(1)
    except Exception as e:
        logging.error(f"Error combining FASTA files: {e}")
        sys.exit(1)
    logging.info(f"Combined FASTA file created at: {combined_fasta}")
    logging.info(f"Group file created at: {group_file}")

def validate_output_files(combined_fasta, group_file):
    """
    Validates the existence and non-emptiness of combined FASTA and group files.
    """
    logging.info("Validating combined FASTA and group files...")
    if not os.path.exists(combined_fasta):
        logging.error(f"Combined FASTA file not found: {combined_fasta}")
        sys.exit(1)
    if not os.path.exists(group_file):
        logging.error(f"Group file not found: {group_file}")
        sys.exit(1)
    if os.path.getsize(combined_fasta) == 0:
        logging.error(f"Combined FASTA file is empty: {combined_fasta}")
        sys.exit(1)
    if os.path.getsize(group_file) == 0:
        logging.error(f"Group file is empty: {group_file}")
        sys.exit(1)
    logging.info("Combined FASTA and group files validated successfully.")

def run_pipeline(args):
    """
    Runs the entire 16S rRNA sequencing analysis pipeline.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    log_file = os.path.join(args.output_dir, "pipeline.log")
    setup_logging(log_file)

    # Skera Split
    if args.adapters_fasta and not args.skip_skera:
        skera_out = os.path.join(args.output_dir, "skera_split.bam")
        skera_split(args.input_bam, args.adapters_fasta, skera_out)
        input_bam = skera_out
    else:
        input_bam = args.input_bam
        logging.info("Skera Split skipped.")

    # Lima Demultiplexing
    lima_out_prefix = os.path.join(args.output_dir, "lima_output")
    lima_demultiplex(input_bam, args.barcodes_fasta, lima_out_prefix, args.barcode_type)

    # Find demultiplexed BAM files
    demux_bams = glob.glob(f"{lima_out_prefix}.*.bam")
    if not demux_bams:
        logging.error("No demultiplexed BAM files found.")
        sys.exit(1)

    # Move invalid files and generate report, get valid BAM files based on filename pattern and size
    invalid_dir = os.path.join(args.output_dir, "invalid_files")
    report_file = os.path.join(args.output_dir, "invalid_files_report.csv")
    valid_bams = move_invalid_files(demux_bams, invalid_dir, report_file)

    # Ensure sample-barcode.csv exists and is tab-delimited
    if args.sample_barcode_csv:
        converted_csv = os.path.join(args.output_dir, "sample-barcode_converted.tsv")
        ensure_tab_delimited(args.sample_barcode_csv, converted_csv)
        logging.info(f"Using converted sample-barcode CSV file: {converted_csv}")
    else:
        logging.error("Sample barcode CSV file is required.")
        sys.exit(1)

    # Read barcode CSV to map barcodes to sample names
    barcode_to_sample = {}
    try:
        with open(converted_csv, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')  # Now tab-delimited
            for row in reader:
                barcode = row['Barcode'].strip()
                sample_name = row['Sample Name'].strip()
                barcode_to_sample[barcode] = sample_name
        logging.info(f"Total Barcodes mapped: {len(barcode_to_sample)}")
    except KeyError as e:
        logging.error(f"Missing expected column in barcode CSV: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error reading converted barcode CSV: {e}")
        sys.exit(1)

    # Further filter valid BAMs based on barcodes present in sample-barcode.csv
    final_valid_bams = []
    for bam in valid_bams:
        bam_name = os.path.basename(bam)
        barcode, status = extract_primer_directions(bam_name)
        if barcode in barcode_to_sample:
            final_valid_bams.append(bam)
        else:
            # Move to invalid_files with reason 'Barcode not in sample-barcode.csv'
            reason = "Barcode not in sample-barcode.csv"
            logging.info(f"BAM file '{bam}' has barcode '{barcode}' not listed in sample-barcode.csv. Moving to invalid_files.")
            # Update report
            with open(report_file, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=['BAM File', 'Issue'])
                writer.writerow({'BAM File': bam, 'Issue': reason})
            # Move BAM and related files
            target_subdir = "Barcode_not_in_sample-barcode.csv"
            target_dir = os.path.join(invalid_dir, target_subdir)
            os.makedirs(target_dir, exist_ok=True)
            for file_ext in ['', '.pbi']:
                file_pattern = f"{bam}{file_ext}"
                if os.path.exists(file_pattern):
                    try:
                        shutil.move(file_pattern, target_dir)
                        logging.info(f"Moved {file_pattern} to {target_dir}")
                    except shutil.Error as move_err:
                        logging.error(f"Error moving {file_pattern} to {target_dir}: {move_err}")
            # Handle XML files separately
            xml_pattern = f"{bam.replace('.bam', '')}*.xml"
            for file in glob.glob(xml_pattern):
                if os.path.exists(file):
                    try:
                        shutil.move(file, target_dir)
                        logging.info(f"Moved {file} to {target_dir}")
                    except shutil.Error as move_err:
                        logging.error(f"Error moving {file} to {target_dir}: {move_err}")

    # Replace valid_bams with final_valid_bams
    valid_bams = final_valid_bams

    # Proceed only with valid BAM files
    if valid_bams:
        logging.info("Converting BAM files to FASTA/FASTQ...")
        convert_bam_to_fastx(
            valid_bams,
            output_dir=args.output_dir,
            convert_types=args.convert_types,
            compression_level=args.compression_level,
            uncompressed=args.uncompressed
        )
        logging.info("BAM to FASTA/FASTQ conversion completed.")

        # Generate FASTA index files
        if 'fasta' in args.convert_types:
            fasta_files = [os.path.join(args.output_dir, f"{os.path.splitext(os.path.basename(bam))[0]}.fasta") for bam in valid_bams]
            logging.info("Generating FASTA index files...")
            generate_fasta_index(fasta_files)
            logging.info("FASTA indexing completed.")

        # Combine FASTA files and generate group file
        combined_fasta = os.path.join(args.output_dir, "combined.fasta")
        combined_group_file = os.path.join(args.output_dir, "combined.groups")
        logging.info("Combining FASTA files and generating group file...")
        
        # List of valid FASTA files
        valid_fasta_files = [os.path.basename(fast) for fast in fasta_files if os.path.exists(fast)]
        
        # Combine FASTA files using the updated function
        combine_fasta(combined_fasta, combined_group_file, barcode_to_sample, valid_fasta_files, args.output_dir)

        # Validate combined FASTA and group files
        validate_output_files(combined_fasta, combined_group_file)
    else:
        logging.warning("No valid BAM files to process after barcode validation.")

    logging.info("Pipeline execution completed.")

def main():
    """
    Parses command-line arguments and initiates the pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Pipeline: Skera Split, Lima Demultiplexing, Validate and Move Invalid Files, Convert BAM to FASTA/FASTQ, Generate FASTA Index, and Combine FASTA Files for 16S rRNA Sequencing Analysis"
    )
    parser.add_argument('--input-bam', required=True, help="Input BAM file.")
    parser.add_argument('--adapters-fasta', help="Adapters FASTA for Skera Split.")
    parser.add_argument('--barcodes-fasta', required=True, help="Barcodes FASTA for Lima Demultiplexing.")
    parser.add_argument('--output-dir', required=True, help="Output directory.")
    parser.add_argument('--athena-db', required=True, help="Athena database directory.")
    parser.add_argument('--barcode-type', choices=['symmetric', 'asymmetric'], required=True, help="Barcode type for Lima.")
    parser.add_argument('--skip-skera', action='store_true', help="Skip Skera Split step.")

    # Conversion options
    parser.add_argument('--convert-types', nargs='+', choices=['fasta', 'fastq'], default=['fasta', 'fastq'], help="Types to convert BAM files into.")
    parser.add_argument('--compression-level', type=int, choices=range(1, 10), default=1, help="Compression level for output files.")
    parser.add_argument('--uncompressed', action='store_true', help="Output uncompressed FASTA/FASTQ files.")

    # Mothur processing parameters are removed from this script

    parser.add_argument('--sample-barcode-csv', required=True, help="CSV file with 'Barcode' and 'Sample Name' columns, can be comma or tab-delimited.")

    args = parser.parse_args()

    # Validate arguments
    if not args.sample_barcode_csv:
        parser.error("The --sample-barcode-csv argument is required.")

    try:
        run_pipeline(args)
        logging.info("Pipeline completed successfully.")
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
