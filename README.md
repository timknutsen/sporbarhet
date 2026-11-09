# sporbarhet

Extract parent genotypes for salmon traceability (sporbarhet) deliveries on the reduced
~20k SNP panel, in the same ped/map format as the 2025 and 2026 deliveries.

```bash
sporbarhet INFILE OUT_PREFIX [--force] [--keep-temp]

sporbarhet hann-foreldre 2026_Cermaq-Q3-2025_to_Q3-2026
sporbarhet hann-foreldre leveranse/2026_Salmar-Q3-2025_to_Q3-2026   # creates the folder
sbatch --wrap "sporbarhet hann-foreldre 2026_Cermaq-Q3-2025_to_Q3-2026"
```

- `INFILE`: one PIT tag per line. Excel/Windows line endings, a BOM, trailing spaces,
  blank lines and duplicate lines are cleaned automatically.
- `OUT_PREFIX`: output prefix, optionally with a folder. Letters, digits and `. _ + -` only.
- Output: `OUT_PREFIX.{ped,map,bed,bim,fam,bim.backup_affx}` and `OUT_PREFIX.<timestamp>.log`.
- Exit codes: 0 ok, 1 input or environment problem, 2 pipeline failed, 3 format check failed.

## Setup (once per user)

1. Add the command to your PATH in `~/.bashrc`:
   ```bash
   export PATH=/mnt/efshome/aquagen/code/timknu/sporbarhet/bin:$PATH
   ```
2. You need `~/.my.cnf` with read access to the genotype database (genodb_v5).

No `module load` or `conda activate` is needed. The command uses a fixed Python
environment and ignores whatever modules or conda environments are loaded.

## What it does

`bin/sporbarhet` is a wrapper around `pipeline_sporbarhet_v2.py`:

1. `aqg_db2plink` fetches genotypes by PIT tag (status OK and LOWCR, SIMONv31 map, outer join of chips)
2. KING `--duplicate` confirms that samples sharing a PIT tag are genetically identical; stops if not
3. one sample per PIT tag is kept; PLINK subsets to `AquaGen_snps.reduced-set.affx`
4. SNP IDs are converted from Affx-* to AX-* (Axiom Ssa70kv2 annotation)
5. PLINK writes tab-separated ped/map; sex is set to 1

The wrapper adds: a pinned interpreter, a clean environment, input cleaning, a private
work folder per run (parallel runs in one folder are safe), refusal to overwrite, and a
format check of the result (tab ped/map, FID/PAT/MAT=0, SEX=1, PHENO=-9, AX-* IDs,
chromosomes 1-29, every input PIT tag present exactly once).

## Environment requirements

| Requirement | Default | Override |
|---|---|---|
| Python >= 3.12 with pandas and aqg_tools | `/mnt/efshome/applications/miniconda3/envs/rust-185` | `SPORBARHET_ENV` |
| DB credentials | `~/.my.cnf` of the user running it | none |
| PLINK 1.9, KING, SNP list | `/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet/` | edit the pipeline |
| Axiom annotation | `/mnt/efshome/aquagen/viktige_filer/Affy_annotation_files/Ssa70kv2/` | edit the pipeline |

Python 3.7 and 3.9 cannot run the pipeline (f-string syntax needs 3.12).

`rust-185` is a symlink to a personal conda environment (`/mnt/efshome/home/fabgra/.conda/envs/rust-185`)
and `aqg_tools` is installed in editable mode from a personal source tree. If either moves or
breaks, point `SPORBARHET_ENV` at any environment with Python >= 3.12, pandas and aqg_tools.

## Tests

`tests/run_regression.sh` reruns the 2025 Mowi and 2026 Cermaq lists in a clean environment,
a polluted environment, with messy input and with an unknown ID, and compares the output to
the delivered files. It prints counts only.
