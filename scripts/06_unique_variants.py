#!/usr/bin/env python3
"""
  06_unique_variants.py VCF CALLER OUT_PREFIX [--min-dp 5] [--min-mq 50] [--max-fs 20] [--max-sor 3] [--min-af F]
"""
import sys, gzip, math, argparse

TABLE_SIZE = 200


def opn(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def fisher_two_sided(a, b, c, d):
    r1, c1, n = a + b, a + c, a + b + c + d
    lf = lambda k: math.lgamma(k + 1)
    const = lf(r1) + lf(n - r1) + lf(c1) + lf(n - c1) - lf(n)
    logp = lambda x: const - lf(x) - lf(r1 - x) - lf(c1 - x) - lf(n - r1 - c1 + x)
    obs = logp(a)
    p = sum(math.exp(logp(x)) for x in range(max(0, r1 + c1 - n), min(r1, c1) + 1) if logp(x) <= obs + 1e-7)
    return min(p, 1.0)


def fisher_strand(rf, rr, af, ar):
    t = rf + rr + af + ar
    if t > TABLE_SIZE:
        s = TABLE_SIZE / t
        rf, rr, af, ar = (int(x * s) for x in (rf, rr, af, ar))
    p = fisher_two_sided(rf, rr, af, ar)
    return -10 * math.log10(max(p, 1e-320))


def strand_odds_ratio(rf, rr, af, ar):
    rf, rr, af, ar = rf + 1, rr + 1, af + 1, ar + 1
    r = (rf * ar) / (rr * af)
    return math.log(r + 1 / r) + math.log(min(rf, rr) / max(rf, rr)) - math.log(min(af, ar) / max(af, ar))


def ints(x):
    return [int(v) if v not in ('.', '') else 0 for v in x.split(',')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf'); ap.add_argument('caller'); ap.add_argument('out')
    ap.add_argument('--min-dp', type=int, default=5)
    ap.add_argument('--min-mq', type=float, default=50)
    ap.add_argument('--max-fs', type=float, default=20)
    ap.add_argument('--max-sor', type=float, default=3)
    ap.add_argument('--min-af', type=float, default=None)
    a = ap.parse_args()

    n = {'sites': 0, 'site_filter': 0, 'no_carrier': 0, 'shared': 0, 'unique': 0}
    with opn(a.vcf) as vcf, open(a.out + '.tsv', 'w') as tsv, open(a.out + '.sites.vcf', 'w') as sites:
        tsv.write('\t'.join(['CALLER', 'CHROM', 'POS', 'REF', 'ALT', 'LINE', 'DP', 'AD_ALT', 'AF', 'MQ', 'FS', 'SOR']) + '\n')
        for l in vcf:
            if l.startswith('##'):
                if l.startswith('##contig') or l.startswith('##fileformat'):
                    sites.write(l)
                continue
            f = l.rstrip('\n').split('\t')
            if l.startswith('#CHROM'):
                samples = f[9:]
                sites.write('##INFO=<ID=LINE,Number=1,Type=String,Description="Line carrying the variant">\n')
                sites.write('##INFO=<ID=AF_LINE,Number=1,Type=Float,Description="Alternative read fraction in that line">\n')
                sites.write('\t'.join(f[:8]) + '\n')
                continue
            if ',' in f[4]:
                sys.exit('Multiallelic record at %s:%s, split it with bcftools norm -m -any first' % (f[0], f[1]))
            n['sites'] += 1
            info = dict(x.split('=', 1) if '=' in x else (x, True) for x in f[7].split(';'))
            keys = f[8].split(':')
            gts = [dict(zip(keys, s.split(':'))) for s in f[9:]]

            mq = float(info['MQ']) if 'MQ' in info else None
            if 'FS' in info and 'SOR' in info:
                fs, sor = float(info['FS']), float(info['SOR'])
            else:
                rf = rr = af_ = ar = 0
                for g in gts:
                    if 'ADF' in g and 'ADR' in g:
                        F, R = ints(g['ADF']), ints(g['ADR'])
                        rf, af_, rr, ar = rf + F[0], af_ + F[1], rr + R[0], ar + R[1]
                fs, sor = fisher_strand(rf, rr, af_, ar), strand_odds_ratio(rf, rr, af_, ar)
            if (mq is not None and mq < a.min_mq) or fs > a.max_fs or sor > a.max_sor:
                n['site_filter'] += 1
                continue

            carriers = []
            for s, g in zip(samples, gts):
                ad = ints(g.get('AD', '0,0'))
                dp = sum(ad) if 'DP' not in g or g['DP'] in ('.', '') else int(g['DP'])
                frac = ad[1] / sum(ad) if len(ad) > 1 and sum(ad) else 0.0
                if dp < a.min_dp:
                    continue
                if a.min_af is not None:
                    alt = len(ad) > 1 and ad[1] > 0 and frac >= a.min_af
                else:
                    alt = '1' in g.get('GT', '.').replace('|', '/').split('/')
                if alt:
                    carriers.append((s, dp, ad[1] if len(ad) > 1 else 0, frac))
            if not carriers:
                n['no_carrier'] += 1
                continue
            if len(carriers) > 1:
                n['shared'] += 1
                continue
            n['unique'] += 1
            s, dp, alt, frac = carriers[0]
            mqs = '%.1f' % mq if mq is not None else 'NA'
            tsv.write('\t'.join([a.caller, f[0], f[1], f[3], f[4], s, str(dp), str(alt), '%.4f' % frac, mqs, '%.2f' % fs, '%.3f' % sor]) + '\n')
            sites.write('\t'.join([f[0], f[1], '.', f[3], f[4], '.', '.', 'LINE=%s;AF_LINE=%.4f' % (s, frac)]) + '\n')
    print('%s: %d sites, %d removed by MQ/FS/SOR, %d without a carrier at DP >= %d, %d shared by several lines, %d unique'
          % (a.caller, n['sites'], n['site_filter'], n['no_carrier'], a.min_dp, n['shared'], n['unique']))


if __name__ == '__main__':
    main()
