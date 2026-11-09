#!/usr/bin/env bash
# Regression + environment tests for the sporbarhet CLI. Prints summaries only (no IDs, no genotypes).
# usage: run_regression.sh [WORKDIR]
set -uo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
S=$(readlink -f "$here/../bin/sporbarhet")
CMP="/mnt/efshome/applications/miniconda3/envs/rust-185/bin/python $here/compare_plink.py"
TR=/mnt/efshome/aquagen/projects/6136_TRACK/sporbarhet
W=${1:-$here/run_regression}
rm -rf "$W"; mkdir -p "$W"; cd "$W" || exit 1

show() { grep -h '^sporbarhet:' "$1" | sed 's/^/    /'; }

# R1: completely clean environment (no profile, minimal PATH) - like cron or a bare SLURM job
echo "R1 clean env (env -i), 2026 Cermaq list"
env -i HOME="$HOME" USER="$USER" PATH=/usr/bin:/bin "$S" "$TR/2026/2026_Cermaq_Q3-2025_to_Q3-2026/hann-foreldre" r1/2026_Cermaq > r1.out 2>&1; echo "    rc=$?"; show r1.out

# R2: poisoned environment - old Python module, its PYTHONPATH, conda base active
echo "R2 poisoned env (python/3.7_MKL + PYTHONPATH + conda base), 2025 Mowi list"
bash -lc "module load mariadb/10.5.10 python/3.7_MKL plink/1.90beta king/2.2.8 >/dev/null 2>&1
  source /mnt/efshome/applications/miniconda3/etc/profile.d/conda.sh && conda activate base
  export PYTHONPATH=/mnt/efshome/applications/MKL/intelpython/python3.7/lib/python3.7/site-packages
  echo \"    shell python: \$(python --version 2>&1), PYTHONPATH set\"
  '$S' '$TR/2025_Mowi-Q1-2021_to_Q3-2025/hann-foreldre' r2/2025_Mowi" > r2.out 2>&1; echo "    rc=$?"; grep -h 'shell python' r2.out; show r2.out

# R3 + R4 run concurrently in the same folder (tests isolation of the pipeline's fixed-name temp files)
mkdir -p r34
{ printf '\xEF\xBB\xBF'; tr -d '\r' < "$TR/2025_Mowi-Q1-2021_to_Q3-2025/hann-foreldre" | awk 'NF{print $1"  \r"} NR==2{print ""; print $1"\r"}'; } > r34/messy_ids.txt
{ cat "$TR/2025_Mowi-Q1-2021_to_Q3-2025/hann-foreldre"; echo "FFFFFFFFFF"; } > r34/with_unknown.txt
echo "R3 messy input (BOM, CRLF, trailing spaces, blank line, duplicate) || R4 unknown ID, same folder, concurrently"
"$S" r34/messy_ids.txt r34/2025_Mowi_messy > r3.out 2>&1 & p3=$!
"$S" r34/with_unknown.txt r34/2025_Mowi_unknown > r4.out 2>&1 & p4=$!
wait $p3; echo "    R3 rc=$?"; show r3.out
wait $p4; echo "    R4 rc=$? (expect 2)"; show r4.out

echo
echo "== Compare R1 with Hanne's 29 Oct 2026 delivery"; $CMP r1/2026_Cermaq "$TR/2026/2026_Cermaq_Q3-2025_to_Q3-2026/2026_Cermaq_Q3-2025_to_Q3-2026" | sed 's/^/    /'
echo "== Compare R1 with Tim's 28 Oct 2026 LOWCR_v2 run"; $CMP r1/2026_Cermaq "$TR/2026_Cermaq_Q3-2025_to_Q3-2026/LOWCR_v2/2026_Cermaq-Q3-2025_to_Q3-2026" | sed 's/^/    /'
echo "== Compare R2 with the delivered 2025 Mowi files (Dec 2025, v1 pipeline)"; $CMP r2/2025_Mowi "$TR/2025_Mowi-Q1-2021_to_Q3-2025/2025_Mowi-Q1-2021_to_Q3-2025" | sed 's/^/    /'
echo "== Compare R3 (messy input) with R2 (clean input)"; $CMP r34/2025_Mowi_messy r2/2025_Mowi | sed 's/^/    /'
echo "== Byte-identical? R3 vs R2: ped $(cmp -s r34/2025_Mowi_messy.ped r2/2025_Mowi.ped && echo yes || echo no), map $(cmp -s r34/2025_Mowi_messy.map r2/2025_Mowi.map && echo yes || echo no)"
echo "== Files delivered per run (expect ped map bed bim fam bim.backup_affx + one .log)"
for p in r1/2026_Cermaq r2/2025_Mowi r34/2025_Mowi_messy; do printf '    %-22s %s\n' "$p" "$(ls "$(dirname $p)" | grep "^$(basename $p)\." | sed "s/^$(basename $p)\.//; s/^[0-9_]*\.log$/<ts>.log/" | sort | tr '\n' ' ')"; done
echo "== Leftover temp folders: $(find . -maxdepth 2 -name '.sporbarhet_*' | wc -l) (expect 1, from the failed R4 run)"
