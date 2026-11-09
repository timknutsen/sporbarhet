"""Compare two tab ped/map sets keyed by individual and SNP ID. Prints counts only, never IDs or genotypes.

usage: python compare_plink.py NEW_PREFIX OLD_PREFIX
"""
import sys


def load(prefix):
    snps = [line.split('\t')[1] for line in open(prefix + '.map')]
    geno = {}
    fixed = set()
    for line in open(prefix + '.ped'):
        f = line.rstrip('\n').split('\t')
        fixed.add(tuple(f[0:1] + f[2:6]))
        g = f[6:]
        geno[f[1]] = {s: (g[2 * i], g[2 * i + 1]) for i, s in enumerate(snps)}
    return snps, geno, fixed


new_snps, new, new_fixed = load(sys.argv[1])
old_snps, old, old_fixed = load(sys.argv[2])
common_ind = new.keys() & old.keys()
common_snp = set(new_snps) & set(old_snps)
order_same = [s for s in new_snps if s in common_snp] == [s for s in old_snps if s in common_snp]
cells = diff = new_missing_only = old_missing_only = 0
for i in common_ind:
    a, b = new[i], old[i]
    for s in common_snp:
        cells += 1
        ga, gb = a[s], b[s]
        if sorted(ga) != sorted(gb):
            diff += 1
            if ga == ('0', '0'):
                new_missing_only += 1
            elif gb == ('0', '0'):
                old_missing_only += 1
print(f"individuals: new={len(new)} old={len(old)} common={len(common_ind)}")
print(f"snps: new={len(new_snps)} old={len(old_snps)} common={len(common_snp)} same_order_on_common={order_same}")
print(f"fixed columns (FID,PAT,MAT,SEX,PHENO) value sets equal: {new_fixed == old_fixed}")
print(f"genotype calls compared={cells} differing={diff} "
      f"(missing only in new={new_missing_only}, missing only in old={old_missing_only}, "
      f"discordant called={diff - new_missing_only - old_missing_only})")
