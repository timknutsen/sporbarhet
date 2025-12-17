# Sporbarhet Genotype Extraction Pipeline

Extract and process genotypes from AquaGen database for salmon traceability (sporbarhet).

## What it does

1. Fetches genotypes from AquaGen database (via `aqg_db2plink`)
2. Validates duplicates using KING (if any exist)
3. Subsets to reduced SNP panel (~20k SNPs)
4. Converts SNP IDs from Affx-* → AX-* format
5. Exports ped/map files with sex set to male

## Requirements

- `aqg_db2plink` (rust-185 conda environment)
- PLINK 1.9
- KING
- pandas

## Usage

```bash
# Activate rust environment
conda activate rust-185

# Run pipeline
python /mnt/efshome/aquagen/code/timknu/sporbarhet/get_sporbarhet_genos.py -o <output_prefix> -p <parent_file>
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-o, --output` | Output file prefix (e.g., "2024_Cermaq-Q3") |
| `-p, --parents` | File with PIT tags, one per line |
| `--no-cleanup` | Keep intermediate files for debugging |

### Examples

```bash
# Basic usage
python /mnt/efshome/aquagen/code/timknu/sporbarhet/get_sporbarhet_genos.py -o 2024_Cermaq-Q3 -p hann-foreldre

# Keep temp files
python /mnt/efshome/aquagen/code/timknu/sporbarhet/get_sporbarhet_genos.py -o 2024_Cermaq-Q3 -p hann-foreldre --no-cleanup
```

## Input

Parent file with one PIT tag per line:
```
041A559AC3
041A5591E1
```

## Output

| File | Description |
|------|-------------|
| `<prefix>.bed` | PLINK binary genotypes |
| `<prefix>.bim` | SNP information (AX-* IDs) |
| `<prefix>.fam` | Sample information |
| `<prefix>.ped` | Tab-delimited genotypes |
| `<prefix>.map` | SNP map file |
| `<prefix>.bim.backup_affx` | Original BIM with Affx-* IDs |
| `<prefix>.<timestamp>.log` | Full log file |

## Output format

PED file columns:
```
FID  IID         PAT  MAT  SEX  PHENO  GENOTYPES...
0    041A559AC3  0    0    1    -9     G  T  G  A  ...
```

- Sex is set to 1 (male) for all samples
- SNP IDs use AX-* format (Axiom standard)

## Configuration

Tool paths are configured at the top of the script:
```python
PLINK19 = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/tools/plink"
KING = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/tools/king"
SNPLIST = "/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/AquaGen_snps.reduced-set.affx"
AXIOM_ANNOTATION = "/mnt/efshome/aquagen/viktige_filer/Affy_annotation_files/Ssa70kv2/Axiom_Ssa70kv2_Annotation.r1.csv"
```

## Author

Tim Knutsen, AquaGen

