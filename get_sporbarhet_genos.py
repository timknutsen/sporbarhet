#!/usr/bin/env python3
"""
Genotype extraction pipeline for salmon traceability (sporbarhet)

This pipeline:
1. Fetches genotypes from AquaGen database (via aqg_db2plink)
2. Identifies true duplicates using KING
3. Validates that all duplicates by ID are confirmed by genotype
4. Removes duplicates, keeping one per animal
5. Subsets to reduced SNP panel
6. Updates SNP IDs from Affx-* -> AX-*
7. Exports ped/map with sex=1 (male)

Usage:
    python pipeline.py --output <output_prefix> --parents <parent_file>

Examples:
    python pipeline.py --output 2024_Cermaq-Q3 --parents hann-foreldre
    python pipeline.py -o 2024_Cermaq-Q3 -p hann-foreldre
    python pipeline.py -o 2024_Cermaq-Q3 -p hann-foreldre --no-cleanup

Requirements:
    - aqg_db2plink (activate rust environment first)
    - PLINK 1.9
    - KING
    - pandas
"""

import argparse
import subprocess
import pandas as pd
import os
import sys
import shutil
import re
import glob
from datetime import datetime

# Fixed paths - DO NOT CHANGE
PLINK19 = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/tools/plink"
KING = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/tools/king"
SNPLIST = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/AquaGen_snps.reduced-set.affx"
AXIOM_ANNOTATION = "/mnt/efshome/aquagen/viktige_filer/Affy_annotation_files/Ssa70kv2/Axiom_Ssa70kv2_Annotation.r1.csv"

# Global variables set by main()
OUTFILES_CORE = None
PARENT_FILE = None
DO_CLEANUP = True


class Logger:
    """Handles logging to both terminal and file"""
    
    def __init__(self, log_file_path):
        self.log_file_path = log_file_path
        self.log_file = open(log_file_path, 'w')
        
    def log(self, msg, level="INFO"):
        """Print timestamped message to terminal and log file"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        formatted = f"[{timestamp}] [{level}] {msg}"
        print(formatted)
        self.log_file.write(formatted + '\n')
        self.log_file.flush()
    
    def info(self, msg):
        self.log(msg, "INFO")
    
    def warn(self, msg):
        self.log(msg, "WARN")
    
    def error(self, msg):
        self.log(msg, "ERROR")
    
    def print_raw(self, msg):
        """Print raw message without timestamp (for data output)"""
        print(msg)
        self.log_file.write(msg + '\n')
        self.log_file.flush()
    
    def close(self):
        self.log_file.close()


# Global logger instance
logger = None


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Genotype extraction pipeline for salmon traceability',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python pipeline.py --output 2024_Cermaq-Q3 --parents hann-foreldre
    python pipeline.py -o 2024_Cermaq-Q3 -p hann-foreldre
    python pipeline.py -o 2024_Cermaq-Q3 -p hann-foreldre --no-cleanup
    
Notes:
    - Parent file should contain one PIT tag per line
    - Output files will be created in current directory
    - A log file will be created: <output>.<timestamp>.log
    - Use --no-cleanup to keep intermediate files for debugging
        """
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file prefix (e.g., "2024_Cermaq-Q3")'
    )
    
    parser.add_argument(
        '-p', '--parents',
        required=True,
        help='Parent file with PIT tags (one per line)'
    )
    
    parser.add_argument(
        '--no-cleanup',
        dest='cleanup',
        action='store_false',
        default=True,
        help='Keep intermediate files for debugging'
    )
    
    return parser.parse_args()


def validate_inputs():
    """Validate all inputs before running pipeline"""
    global logger
    
    logger.info("Validating inputs...")
    errors = []
    warnings = []
    
    # Check parent file exists
    if not os.path.exists(PARENT_FILE):
        errors.append(f"Parent file not found: {PARENT_FILE}")
    else:
        # Check it's not empty
        with open(PARENT_FILE) as f:
            lines = [l.strip() for l in f if l.strip()]
        if len(lines) == 0:
            errors.append(f"Parent file is empty: {PARENT_FILE}")
        else:
            logger.info(f"Parent file: {PARENT_FILE} ({len(lines)} PIT tags)")
    
    # Check PLINK exists and is executable
    if not os.path.exists(PLINK19):
        errors.append(f"PLINK not found: {PLINK19}")
    elif not os.access(PLINK19, os.X_OK):
        errors.append(f"PLINK not executable: {PLINK19}")
    else:
        logger.info(f"PLINK: {PLINK19} [OK]")
    
    # Check KING exists and is executable
    if not os.path.exists(KING):
        errors.append(f"KING not found: {KING}")
    elif not os.access(KING, os.X_OK):
        errors.append(f"KING not executable: {KING}")
    else:
        logger.info(f"KING: {KING} [OK]")
    
    # Check aqg_db2plink is available
    result = subprocess.run("which aqg_db2plink", shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        errors.append("aqg_db2plink not found in PATH. Did you activate the rust environment?")
    else:
        logger.info(f"aqg_db2plink: {result.stdout.strip()} [OK]")
    
    # Check SNP list exists
    if not os.path.exists(SNPLIST):
        errors.append(f"SNP list not found: {SNPLIST}")
    else:
        with open(SNPLIST) as f:
            snp_count = sum(1 for _ in f)
        logger.info(f"SNP list: {SNPLIST} ({snp_count} SNPs) [OK]")
    
    # Check Axiom annotation exists
    if not os.path.exists(AXIOM_ANNOTATION):
        errors.append(f"Axiom annotation not found: {AXIOM_ANNOTATION}")
    else:
        logger.info(f"Axiom annotation: {AXIOM_ANNOTATION} [OK]")
    
    # Check output directory is writable
    output_dir = os.path.dirname(os.path.abspath(OUTFILES_CORE)) or '.'
    if not os.access(output_dir, os.W_OK):
        errors.append(f"Output directory not writable: {output_dir}")
    else:
        logger.info(f"Output directory: {output_dir} [OK]")
    
    # Check if output files already exist
    existing_files = []
    for ext in ['bed', 'bim', 'fam', 'ped', 'map']:
        f = f"{OUTFILES_CORE}.{ext}"
        if os.path.exists(f):
            existing_files.append(f)
    
    if existing_files:
        warnings.append(f"Output files already exist and will be overwritten: {', '.join(existing_files)}")
    
    # Print warnings
    for w in warnings:
        logger.warn(w)
    
    # Print errors and exit if any
    if errors:
        logger.error("Validation failed with the following errors:")
        for e in errors:
            logger.error(f"  - {e}")
        logger.error("Please fix the above errors and try again.")
        sys.exit(1)
    
    logger.info("All inputs validated successfully!")
    return True


def run_cmd(cmd, description=None):
    """Run shell command and check for errors"""
    if description:
        logger.info(description)
    logger.info(f"Running: {cmd}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    # Log stdout/stderr
    if result.stdout.strip():
        for line in result.stdout.strip().split('\n'):
            logger.print_raw(f"  {line}")
    
    if result.returncode != 0:
        logger.error(f"Command failed with return code {result.returncode}")
        if result.stderr.strip():
            logger.error("STDERR:")
            for line in result.stderr.strip().split('\n'):
                logger.error(f"  {line}")
        raise RuntimeError(f"Command failed with return code {result.returncode}")
    
    return result


def get_base_id(sample_id):
    """
    Extract base ID from sample ID.
    E.g., '041A559AC3_1' -> '041A559AC3'
          '041A559AC3_0' -> '041A559AC3'
          '041A559AC3' -> '041A559AC3'
    """
    # Remove trailing _N suffix (where N is a number)
    return re.sub(r'_\d+$', '', sample_id)


def step1_fetch_genotypes():
    """Fetch genotypes from database with pit_tag as output ID"""
    logger.info("=" * 60)
    logger.info("Step 1 - Fetching genotypes from database")
    logger.info("=" * 60)
    
    cmd = f"""aqg_db2plink \
        --species Ssa \
        --db v5 \
        --out {OUTFILES_CORE}_raw \
        --status OK \
        --names {PARENT_FILE} \
        --idtype pit_tag \
        --idtype_out pit_tag \
        --genome SIMONv31 \
        --join outer"""
    
    run_cmd(cmd)
    
    # Show results
    fam_file = f"{OUTFILES_CORE}_raw.plink.fam"
    if not os.path.exists(fam_file):
        raise RuntimeError(f"Expected output file not created: {fam_file}")
    
    with open(fam_file) as f:
        lines = f.readlines()
    logger.info(f"Samples in fam file: {len(lines)}")
    for line in lines:
        logger.print_raw(f"  {line.strip()}")


def step2_king_duplicates():
    """Run KING to identify true duplicates"""
    logger.info("=" * 60)
    logger.info("Step 2 - Running KING duplicate check")
    logger.info("=" * 60)
    
    cmd = f"""{KING} \
        -b {OUTFILES_CORE}_raw.plink.bed \
        --sexchr 35 \
        --prefix {OUTFILES_CORE}_raw._DUPcheck \
        --duplicate"""
    
    run_cmd(cmd)
    
    con_file = f"{OUTFILES_CORE}_raw._DUPcheck.con"
    if os.path.exists(con_file):
        logger.info("Duplicate pairs found:")
        with open(con_file) as f:
            for line in f:
                logger.print_raw(f"  {line.strip()}")
    else:
        logger.info("No duplicates found (no .con file created)")


def step3_validate_and_create_removal_list():
    """
    Validate that all ID-based duplicates are confirmed by KING.
    Then create removal list.
    """
    logger.info("=" * 60)
    logger.info("Step 3 - Validating duplicates and creating removal list")
    logger.info("=" * 60)
    
    # Read FAM file to find ID-based duplicates
    fam_file = f"{OUTFILES_CORE}_raw.plink.fam"
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None, names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    
    # Extract base IDs and find duplicates
    fam['base_id'] = fam['IID'].apply(get_base_id)
    
    # Group by base_id to find which ones have multiple samples
    id_counts = fam.groupby('base_id').size()
    duplicated_base_ids = id_counts[id_counts > 1].index.tolist()
    
    logger.info(f"Total samples in FAM: {len(fam)}")
    logger.info(f"Unique base IDs: {len(id_counts)}")
    logger.info(f"Base IDs with duplicates: {len(duplicated_base_ids)}")
    
    if len(duplicated_base_ids) == 0:
        logger.info("No ID-based duplicates found. Nothing to validate.")
        # Create empty removal file
        open("remove_dups.txt", 'w').close()
        return
    
    # Show ID-based duplicates
    logger.info("ID-based duplicate groups:")
    expected_pairs = set()
    for base_id in duplicated_base_ids:
        samples = fam[fam['base_id'] == base_id]['IID'].tolist()
        logger.print_raw(f"  {base_id}: {samples}")
        # Create expected pairs (all combinations should be confirmed by KING)
        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                expected_pairs.add((samples[i], samples[j]))
    
    logger.info(f"Expected duplicate pairs to be confirmed: {len(expected_pairs)}")
    
    # Read KING results
    con_file = f"{OUTFILES_CORE}_raw._DUPcheck.con"
    if not os.path.exists(con_file):
        logger.error("=" * 60)
        logger.error("CRITICAL ERROR: DUPLICATE VALIDATION FAILED")
        logger.error("=" * 60)
        logger.error(f"Found {len(duplicated_base_ids)} base IDs with multiple samples,")
        logger.error("but KING did not find ANY duplicate pairs!")
        logger.error("")
        logger.error("This means samples with the same base ID are NOT genetically identical.")
        logger.error("This could indicate:")
        logger.error("  - Sample mix-up or mislabeling")
        logger.error("  - Genotyping errors")
        logger.error("  - Wrong samples in the database")
        logger.error("")
        logger.error("Affected base IDs:")
        for base_id in duplicated_base_ids:
            samples = fam[fam['base_id'] == base_id]['IID'].tolist()
            logger.error(f"  {base_id}: {samples}")
        logger.error("")
        logger.error("Please investigate these samples before proceeding!")
        raise RuntimeError("Duplicate validation failed: ID-based duplicates not confirmed by genotypes")
    
    # Parse KING output
    king_df = pd.read_csv(con_file, sep='\t')
    confirmed_pairs = set()
    for _, row in king_df.iterrows():
        pair = tuple(sorted([row['ID1'], row['ID2']]))
        confirmed_pairs.add(pair)
    
    logger.info(f"KING confirmed duplicate pairs: {len(confirmed_pairs)}")
    
    # Check that all expected pairs are confirmed
    # Normalize expected pairs for comparison
    expected_pairs_normalized = set(tuple(sorted(p)) for p in expected_pairs)
    
    missing_pairs = expected_pairs_normalized - confirmed_pairs
    
    if missing_pairs:
        logger.error("=" * 60)
        logger.error("CRITICAL ERROR: DUPLICATE VALIDATION FAILED")
        logger.error("=" * 60)
        logger.error("The following ID-based duplicate pairs were NOT confirmed by KING:")
        logger.error("")
        for pair in missing_pairs:
            logger.error(f"  {pair[0]} <-> {pair[1]}")
        logger.error("")
        logger.error("This means these samples have the same base ID but different genotypes!")
        logger.error("This could indicate:")
        logger.error("  - Sample mix-up or mislabeling")
        logger.error("  - Genotyping errors")
        logger.error("  - Wrong samples in the database")
        logger.error("")
        logger.error("Please investigate these samples before proceeding!")
        raise RuntimeError(f"Duplicate validation failed: {len(missing_pairs)} ID-based duplicate pairs not confirmed by genotypes")
    
    logger.info("All ID-based duplicates confirmed by KING genotype comparison!")
    
    # Create removal list - keep first sample of each base_id, remove the rest
    to_remove = []
    for base_id in duplicated_base_ids:
        samples = fam[fam['base_id'] == base_id]['IID'].tolist()
        # Sort to ensure consistent ordering, keep first (e.g., _0 before _1)
        samples_sorted = sorted(samples, key=lambda x: (len(x), x))
        keep = samples_sorted[0]
        remove = samples_sorted[1:]
        logger.info(f"  {base_id}: keeping '{keep}', removing {remove}")
        for r in remove:
            to_remove.append({'FID': '0', 'IID': r})
    
    # Write removal file
    remove_df = pd.DataFrame(to_remove)
    remove_df.to_csv("remove_dups.txt", sep=' ', header=False, index=False)
    
    logger.info(f"Removal list created: {len(to_remove)} samples to remove")


def step4_subset_snps():
    """Subset SNPs and remove duplicate samples"""
    logger.info("=" * 60)
    logger.info("Step 4 - PLINK subset to reduced SNP set and remove duplicates")
    logger.info("=" * 60)
    
    cmd = f"""{PLINK19} \
        --bfile {OUTFILES_CORE}_raw.plink \
        --remove remove_dups.txt \
        --extract {SNPLIST} \
        --make-bed \
        --dog \
        --out {OUTFILES_CORE}"""
    
    run_cmd(cmd)
    
    # Show results
    logger.info("Samples after filtering:")
    with open(f"{OUTFILES_CORE}.fam") as f:
        for line in f:
            logger.print_raw(f"  {line.strip()}")


def step5_update_snp_ids():
    """Update SNP IDs in BIM file from Affx-* to AX-*"""
    logger.info("=" * 60)
    logger.info("Step 5 - Updating SNP IDs from Affx-* to AX-*")
    logger.info("=" * 60)
    
    # Load Axiom annotation (skip comment lines)
    logger.info(f"Loading Axiom annotation from {AXIOM_ANNOTATION}")
    axiom = pd.read_csv(AXIOM_ANNOTATION, comment='#', dtype=str, usecols=[0, 1])
    axiom.columns = ['AX', 'Affx']
    axiom['Affx'] = axiom['Affx'].str.strip()
    axiom['AX'] = axiom['AX'].str.strip()
    
    # Create mapping dict
    affx_to_ax = dict(zip(axiom['Affx'], axiom['AX']))
    logger.info(f"Loaded {len(affx_to_ax)} Affx -> AX mappings")
    
    # Read BIM file
    bim_file = f"{OUTFILES_CORE}.bim"
    bim_backup = f"{OUTFILES_CORE}.bim.backup_affx"
    
    bim = pd.read_csv(bim_file, sep='\t', header=None, 
                      names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2'])
    
    original_ids = bim['snp_id'].copy()
    
    # Map Affx -> AX
    bim['snp_id'] = bim['snp_id'].map(affx_to_ax).fillna(bim['snp_id'])
    
    # Count conversions
    converted = (original_ids != bim['snp_id']).sum()
    unchanged = (original_ids == bim['snp_id']).sum()
    
    logger.info(f"SNP IDs converted: {converted}")
    logger.info(f"SNP IDs unchanged: {unchanged}")
    
    if unchanged > 0:
        # Show sample of unchanged IDs
        unchanged_ids = original_ids[original_ids == bim['snp_id']].head(5).tolist()
        logger.warn(f"Sample unchanged IDs (not in mapping): {unchanged_ids}")
    
    # Backup original BIM
    logger.info(f"Backing up original BIM to {bim_backup}")
    shutil.copy2(bim_file, bim_backup)
    
    # Write updated BIM
    bim.to_csv(bim_file, sep='\t', header=False, index=False)
    logger.info(f"Updated BIM written to {bim_file}")
    
    # Show sample of converted IDs
    logger.info("Sample SNP ID conversions:")
    sample = pd.DataFrame({
        'Original': original_ids.head(5),
        'New': bim['snp_id'].head(5)
    })
    logger.print_raw(sample.to_string(index=False))


def step6_export_ped_map():
    """Export to ped/map format"""
    logger.info("=" * 60)
    logger.info("Step 6 - Writing ped and map")
    logger.info("=" * 60)
    
    cmd = f"""{PLINK19} \
        --bfile {OUTFILES_CORE} \
        --recode tabx \
        --dog \
        --out {OUTFILES_CORE}"""
    
    run_cmd(cmd)
    
    logger.info("ped and map written:")
    for ext in ['ped', 'map']:
        f = f"{OUTFILES_CORE}.{ext}"
        size = os.path.getsize(f)
        logger.print_raw(f"  {f}: {size/1024:.1f}K")


def step7_set_sex():
    """Set sex column to 1 (male) in ped file"""
    logger.info("=" * 60)
    logger.info("Step 7 - Setting sex column (5) to 1 (male) in ped")
    logger.info("=" * 60)
    
    ped_file = f"{OUTFILES_CORE}.ped"
    
    # Read ped file
    with open(ped_file, 'r') as f:
        lines = f.readlines()
    
    # Update sex column (column index 4, 0-based)
    updated_lines = []
    for line in lines:
        cols = line.strip().split('\t')
        cols[4] = '1'  # Set sex to male
        updated_lines.append('\t'.join(cols) + '\n')
    
    # Write back
    with open(ped_file, 'w') as f:
        f.writelines(updated_lines)
    
    logger.info("Sex column updated. Ped preview (first 10 columns):")
    for line in updated_lines[:2]:
        cols = line.strip().split('\t')[:10]
        logger.print_raw(f"  {'\t'.join(cols)}")


def cleanup():
    """Remove intermediate files"""
    logger.info("=" * 60)
    logger.info("Cleaning up intermediate files")
    logger.info("=" * 60)
    
    patterns = [
        # Raw plink files from aqg_db2plink
        f"{OUTFILES_CORE}_raw.plink.bed",
        f"{OUTFILES_CORE}_raw.plink.bim",
        f"{OUTFILES_CORE}_raw.plink.fam",
        f"{OUTFILES_CORE}_raw.plink.log",
        f"{OUTFILES_CORE}_raw.plink.nosex",
        f"{OUTFILES_CORE}_raw.db2plink.log",
        f"{OUTFILES_CORE}_raw.meta",
        # KING files (all possible outputs)
        f"{OUTFILES_CORE}_raw._DUPcheck.con",
        f"{OUTFILES_CORE}_raw._DUPcheck.kin",
        f"{OUTFILES_CORE}_raw._DUPcheck.kin0",
        f"{OUTFILES_CORE}_raw._DUPcheckallsegs.txt",
        f"{OUTFILES_CORE}_raw._DUPcheck.seg",
        # PLINK temp files
        f"{OUTFILES_CORE}.nosex",
        f"{OUTFILES_CORE}.log",
        # Other temp files
        "remove_dups.txt",
    ]
    
    # Also use glob to catch any other KING output files
    king_glob = glob.glob(f"{OUTFILES_CORE}_raw._DUPcheck*")
    
    all_files = set(patterns + king_glob)
    
    for f in all_files:
        if os.path.exists(f):
            os.remove(f)
            logger.info(f"Removed {f}")


def main():
    global OUTFILES_CORE, PARENT_FILE, DO_CLEANUP, logger
    
    # Parse command line arguments
    args = parse_args()
    OUTFILES_CORE = args.output
    PARENT_FILE = args.parents
    DO_CLEANUP = args.cleanup
    
    # Create log file with timestamp
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = f"{OUTFILES_CORE}.{timestamp}.log"
    logger = Logger(log_file)
    
    logger.info("=" * 60)
    logger.info("GENOTYPE EXTRACTION PIPELINE FOR SALMON TRACEABILITY")
    logger.info("=" * 60)
    logger.info(f"Output prefix: {OUTFILES_CORE}")
    logger.info(f"Parent file: {PARENT_FILE}")
    logger.info(f"Log file: {log_file}")
    logger.info(f"Cleanup: {'ON' if DO_CLEANUP else 'OFF'}")
    logger.info("=" * 60)
    
    try:
        # Validate inputs
        validate_inputs()
        
        # Run pipeline steps
        step1_fetch_genotypes()
        step2_king_duplicates()
        step3_validate_and_create_removal_list()
        step4_subset_snps()
        step5_update_snp_ids()
        step6_export_ped_map()
        step7_set_sex()
        
        if DO_CLEANUP:
            cleanup()
        else:
            logger.info("=" * 60)
            logger.info("Skipping cleanup (--no-cleanup specified)")
            logger.info("=" * 60)
        
        # Final summary
        logger.info("=" * 60)
        logger.info("PIPELINE FINISHED SUCCESSFULLY!")
        logger.info("=" * 60)
        logger.info("Final output files:")
        for ext in ['bed', 'bim', 'fam', 'ped', 'map']:
            f = f"{OUTFILES_CORE}.{ext}"
            if os.path.exists(f):
                size = os.path.getsize(f)
                logger.print_raw(f"  {f}: {size/1024:.1f}K")
        
        # Also note the backup
        backup = f"{OUTFILES_CORE}.bim.backup_affx"
        if os.path.exists(backup):
            logger.info(f"Original BIM (Affx IDs) backed up to: {backup}")
        
        logger.info(f"Full log saved to: {log_file}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        logger.error(f"Check log file for details: {log_file}")
        logger.close()
        sys.exit(1)
    
    logger.close()


if __name__ == "__main__":
    main()
