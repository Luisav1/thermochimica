#!/usr/bin/env python3
import json
import argparse


def load(path):
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def get_liquid_mu(block):
    sp = block.get('solution phases', {}).get('LIQUID', {}).get('species', {})
    return {k: v.get('chemical potential') for k, v in sp.items() if isinstance(v, dict)}


def main():
    ap = argparse.ArgumentParser(description='Diff two thermochimica_out.json files by state index.')
    ap.add_argument('baseline')
    ap.add_argument('candidate')
    args = ap.parse_args()

    b = load(args.baseline)
    c = load(args.candidate)

    keys = sorted(set(b.keys()) & set(c.keys()), key=lambda x: int(x) if str(x).isdigit() else x)
    print('idx,temperature,delta_iter,delta_functional_norm,delta_integral_gibbs,delta_mu(species...)')
    for k in keys:
        bb = b[k]
        cc = c[k]
        temp = bb.get('temperature', cc.get('temperature'))
        di = cc.get('GEM iterations', 0) - bb.get('GEM iterations', 0)
        dfn = cc.get('functional norm', 0.0) - bb.get('functional norm', 0.0)
        dg = cc.get('integral Gibbs energy', 0.0) - bb.get('integral Gibbs energy', 0.0)

        bmu = get_liquid_mu(bb)
        cmu = get_liquid_mu(cc)
        mus = sorted(set(bmu) | set(cmu))
        mu_parts = []
        for s in mus:
            bv = bmu.get(s)
            cv = cmu.get(s)
            if bv is None or cv is None:
                mu_parts.append(f'{s}=NA')
            else:
                mu_parts.append(f'{s}={cv-bv:+.6e}')

        print(f"{k},{temp},{di},{dfn:+.6e},{dg:+.6e}," + ';'.join(mu_parts))


if __name__ == '__main__':
    main()
