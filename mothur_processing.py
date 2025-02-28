#!/usr/bin/env python3

import os
import sys
import argparse
import logging
import subprocess
import shutil

def setup_logging(log_file, log_level='INFO'):
    """
    Sets up logging to both file and stdout with the specified log level.
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        print(f"Invalid log level: {log_level}")
        sys.exit(1)
    
    # Ensure the directory for the log file exists
    log_dir = os.path.dirname(log_file)
    if log_dir and not os.path.exists(log_dir):
        try:
            os.makedirs(log_dir, exist_ok=True)
            print(f"Created log directory: {log_dir}")
        except Exception as e:
            print(f"Failed to create log directory {log_dir}: {e}")
            sys.exit(1)
    
    logging.basicConfig(
        level=numeric_level,
        format='[%(asctime)s] %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # Log the current PATH for debugging purposes
    logging.debug(f"Current PATH: {os.environ.get('PATH')}")

def run_mothur_command(mothur_path, command, description):
    """
    Runs a mothur command using subprocess and handles errors.
    """
    logging.info(f"Running mothur {description}...")
    logging.debug(f"Executing command: {command}")
    
    # Verify that mothur_path exists and is executable
    if not os.path.isfile(mothur_path):
        logging.error(f"mothur executable not found at specified path: {mothur_path}")
        sys.exit(1)
    if not os.access(mothur_path, os.X_OK):
        logging.error(f"mothur executable at {mothur_path} is not executable.")
        sys.exit(1)
    
    try:
        result = subprocess.run(
            command,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logging.debug(result.stdout)
        if result.stderr.strip():
            logging.warning(result.stderr.strip())
        logging.info(f"mothur {description} completed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"mothur {description} failed with error:\n{e.stderr}")
        sys.exit(1)

def summarize_seqs(mothur_path, combined_fasta, processors):
    """
    Runs mothur's summary.seqs command on the combined FASTA file.
    """
    summary_cmd = f'"{mothur_path}" "#summary.seqs(fasta={combined_fasta}, processors={processors})"'
    run_mothur_command(mothur_path, summary_cmd, "summary.seqs")

def classify_seqs(mothur_path, combined_fasta, combined_group_file, reference_fasta, taxonomy_file, method, numwanted, search, processors, output_dir):
    """
    Runs mothur's classify.seqs command on the combined FASTA file.
    """
    classify_cmd = (
        f'"{mothur_path}" "#classify.seqs(fasta={combined_fasta}, '
        f'group={combined_group_file}, '
        f'reference={reference_fasta}, '
        f'taxonomy={taxonomy_file}, '
        f'method={method}, '
        f'numwanted={numwanted}, '
        f'search={search}, '
        f'processors={processors})"'
    )
    run_mothur_command(mothur_path, classify_cmd, "classify.seqs")

def remove_lineage(mothur_path, classified_fasta, classified_taxonomy, lineage_exclude):
    """
    Runs mothur's remove.lineage command to filter unclassified sequences.
    """
    lineage_cmd = (
        f'"{mothur_path}" "#remove.lineage(fasta={classified_fasta}, '
        f'taxonomy={classified_taxonomy}, '
        f'lineage={lineage_exclude})"'
    )
    run_mothur_command(mothur_path, lineage_cmd, "remove.lineage")

def main():
    parser = argparse.ArgumentParser(description="Mothur Processing Script")
    parser.add_argument('--combined-fasta', required=True, help="Path to the combined FASTA file.")
    parser.add_argument('--combined-group', required=True, help="Path to the combined group file.")
    parser.add_argument('--output-dir', required=True, help="Output directory.")
    parser.add_argument('--reference-fasta', required=True, help="Path to reference FASTA file.")
    parser.add_argument('--taxonomy-file', required=True, help="Path to taxonomy file.")
    parser.add_argument('--method', default='knn', help="Method for classification.")
    parser.add_argument('--numwanted', type=int, default=1, help="Number of taxonomic classifications to keep per sequence.")
    parser.add_argument('--search', default='blastplus', help="Search algorithm for classify.seqs.")
    parser.add_argument('--processors', type=int, default=8, help="Number of processors to utilize.")
    parser.add_argument('--remove-lineage', action='store_true', help="Flag to execute remove.lineage command to filter unclassified sequences.")
    parser.add_argument('--lineage-exclude', default='lineage.exclude', help="Filename for lineage exclusion.")
    parser.add_argument('--log-file', required=True, help="Log file name (including path).")
    parser.add_argument('--mothur-path', default=None, help="Optional full path to mothur executable.")
    
    args = parser.parse_args()
    
    # Ensure output directory exists before setting up logging
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory set to: {args.output_dir}")
    
    # Setup logging after ensuring output directory exists
    setup_logging(args.log_file, 'DEBUG')  # Defaulting to DEBUG for detailed logs
    
    logging.info("Starting mothur_processing.py script.")
    
    # Determine mothur_path
    if args.mothur_path:
        mothur_path = args.mothur_path
        if not os.path.isfile(mothur_path):
            logging.error(f"Provided mothur path is not a file: {mothur_path}")
            sys.exit(1)
        if not os.access(mothur_path, os.X_OK):
            logging.error(f"mothur executable at {mothur_path} is not executable.")
            sys.exit(1)
        logging.info(f"Using specified mothur path: {mothur_path}")
    else:
        mothur_path = shutil.which("mothur")
        if not mothur_path:
            logging.error("mothur executable not found in PATH.")
            sys.exit(1)
        else:
            logging.info(f"mothur found at: {mothur_path}")
    
    # Summarize sequences
    summarize_seqs(mothur_path, args.combined_fasta, args.processors)
    
    # Classify sequences
    classify_seqs(
        mothur_path=mothur_path,
        combined_fasta=args.combined_fasta,
        combined_group_file=args.combined_group,
        reference_fasta=args.reference_fasta,
        taxonomy_file=args.taxonomy_file,
        method=args.method,
        numwanted=args.numwanted,
        search=args.search,
        processors=args.processors,
        output_dir=args.output_dir
    )
    
    # Optional: Remove lineage for unclassified sequences
    if args.remove_lineage:
        # Determine the output taxonomy file generated by classify.seqs
        # The output taxonomy file is usually named: <fasta_basename>.<taxonomy_basename>.<method>.taxonomy
        # For example: filtered_combined.athena_v2_2.knn.taxonomy
        fasta_basename = os.path.basename(args.combined_fasta).replace('.fasta', '')
        taxonomy_basename = os.path.basename(args.taxonomy_file).replace('.tax', '')
        method_basename = args.method
        taxonomy_output_filename = f"{fasta_basename}.{taxonomy_basename}.{method_basename}.taxonomy"
        classified_taxonomy = os.path.join(args.output_dir, taxonomy_output_filename)
        
        # Verify that the taxonomy file exists
        if not os.path.isfile(classified_taxonomy):
            logging.error(f"Expected taxonomy output file not found: {classified_taxonomy}")
            logging.error("Please verify that classify.seqs ran correctly and generated the taxonomy file.")
            sys.exit(1)
        else:
            logging.info(f"Found taxonomy file: {classified_taxonomy}")
        
        # Use the original combined_fasta for remove.lineage
        classified_fasta = args.combined_fasta
        
        # Run remove.lineage
        remove_lineage(
            mothur_path=mothur_path,
            classified_fasta=classified_fasta,
            classified_taxonomy=classified_taxonomy,
            lineage_exclude=args.lineage_exclude
        )
    
    logging.info("All mothur processing steps completed successfully.")

if __name__ == "__main__":
    main()
