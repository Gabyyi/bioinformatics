# DNA restriction enzymes - able to link to a specific part of the dna and cut the dna inside that pattern or one of the sides depending on the restriction enzyme type

# digest a piece of text based on a restriction enzyme 
# migrate those pieces to a gel created last lab

# 1. take an arbitary dna sequence from NCBI (either a bacterial genome or a viral genome to use 5 restriction enzymes)
# (enzyme name, recognized sequence, cleavage site)

#         5' GAATTC 5'---G     AATTC---3'
# EcoRI   3' CTTAAG 3'---CTTAA G---5'

#         5'GGATCC 5'---G     GATCC---3'
# BamHI   3'CCTAGG 3'---CCTAG G---5'

#         5'AAGCTT 5'---A     AGCTT---3'
# HindIII 3'TTCGAA 3'---TTCGA A---5'

# TagI    5'TCGA 5'---T    CGA---3'
#         3'AGCT 3'---AGC  T---5'

# HaeIII  5'GGCC 5'---GG CC---3'
#         3'CCGG 3'---CC GG---5'

# Input: recognized sequence for any restriction enzyme a DNA to be digested
# Output: 1. no of cleavages for any restriction enzymes
# 2. cleavage position and length

# the output should also take into consideration combinations of different restriction enzymes

# python lab7_1.py "GCF_000853865.1_ViralProj14960_genomic(1).fna"


import argparse
import csv
from itertools import combinations
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Recognition motif, cut offset on the 5' -> 3' reference strand.
ENZYMES = {
    'EcoRI': ('GAATTC', 1),
    'BamHI': ('GGATCC', 1),
    'HindIII': ('AAGCTT', 1),
    'TaqI': ('TCGA', 1),
    'HaeIII': ('GGCC', 2),
}

def find_default_fasta(base_dir):
    for pattern in ('*.fna', '*.fa', '*.fasta'):
        matches = sorted(base_dir.glob(pattern), key=lambda p: p.name.lower())
        if matches:
            return matches[0]
    raise FileNotFoundError(f'No FASTA file found in {base_dir}. Place a .fna/.fa/.fasta file there or pass one explicitly.')


def read_fasta(path):
    records, seq, header = [], [], None
    with open(path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq).upper()))
                header, seq = line[1:], []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, ''.join(seq).upper()))
    if len(records) != 1:
        raise ValueError(f'Expected one FASTA record, found {len(records)}')
    return records[0]

def enzyme_cuts(dna, motif, offset):
    # Report 0-based cut boundaries: cut p lies between bases p and p+1 (1-based).
    cuts = []
    for i in range(len(dna) - len(motif) + 1):
        if dna[i:i+len(motif)] == motif:
            cuts.append(i + offset)
    return cuts

def digest(dna, names):
    per_enzyme = {name: enzyme_cuts(dna, *ENZYMES[name]) for name in names}
    cuts = sorted(set(p for positions in per_enzyme.values() for p in positions))
    boundaries = [0] + cuts + [len(dna)]
    fragments = [(i + 1, boundaries[i], boundaries[i+1], boundaries[i+1] - boundaries[i])
                 for i in range(len(boundaries) - 1)]
    return per_enzyme, cuts, fragments

def write_csv(path, header, rows):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def plot_gel(results, output):
    fig, ax = plt.subplots(figsize=(max(8, 1.3 * len(results) + 2), 9))
    ax.set_facecolor('#151824')
    fig.patch.set_facecolor('white')
    max_size = max(max(lengths) for _, lengths in results)
    min_size = max(1, min(min(lengths) for _, lengths in results))
    low, high = np.log10(min_size) - .18, np.log10(max_size) + .18
    ax.set_ylim(low, high)
    ax.invert_yaxis()  # larger fragments remain nearer the top
    for lane, (name, lengths) in enumerate(results, start=1):
        counts = {}
        for size in lengths:
            counts[size] = counts.get(size, 0) + 1
        for size, count in counts.items():
            y = np.log10(size)
            ax.plot([lane - .30, lane + .30], [y, y], color='#76e9fc',
                    linewidth=min(2 + 1.8 * np.sqrt(count), 11), alpha=.88)
    ticks = [10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 30000]
    ticks = [t for t in ticks if min_size * .7 <= t <= max_size * 1.4]
    ax.set_yticks(np.log10(ticks), [str(t) for t in ticks], color='black')
    ax.set_xticks(range(1, len(results) + 1), [name.replace('+', '\n+') for name, _ in results], fontsize=8)
    ax.tick_params(axis='x', length=0, pad=10)
    ax.set_xlim(.4, len(results) + .6)
    ax.set_ylabel('Fragment length (bp; log scale)')
    ax.set_title('Simulated restriction digest gel (schematic)')
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description='Restriction digest and schematic gel simulation')
    parser.add_argument('fasta', type=Path, nargs='?', default=None,
                        help='Path to a FASTA file. If omitted, the script uses the first .fna/.fa/.fasta file in the current directory.')
    parser.add_argument('--output', type=Path, default=Path('digest_results'))
    args = parser.parse_args()

    fasta_path = args.fasta if args.fasta is not None else find_default_fasta(Path.cwd())
    if not fasta_path.exists():
        raise FileNotFoundError(f'FASTA file not found: {fasta_path}')

    args.output.mkdir(parents=True, exist_ok=True)
    header, dna = read_fasta(fasta_path)
    if set(dna) - set('ACGTN'):
        raise ValueError('Unexpected characters in FASTA sequence')
    names = list(ENZYMES)
    combos = [group for count in range(1, len(names) + 1) for group in combinations(names, count)]
    summary, cut_rows, fragment_rows, gel_lanes = [], [], [], []
    for group in combos:
        label = '+'.join(group)
        per_enzyme, cuts, fragments = digest(dna, group)
        lengths = [f[3] for f in fragments]
        summary.append([label, len(cuts), len(fragments), min(lengths), max(lengths), ','.join(map(str, cuts))])
        for name, positions in per_enzyme.items():
            for pos in positions:
                cut_rows.append([label, name, pos, pos, pos + 1])
        for num, start, end, length in fragments:
            fragment_rows.append([label, num, start + 1, end, length])
        if len(group) == 1 or len(group) == len(names):
            gel_lanes.append((label, lengths))
    write_csv(args.output / 'summary.csv', ['enzymes', 'unique_cuts', 'fragments', 'smallest_bp', 'largest_bp', 'cut_boundaries_0based'], summary)
    write_csv(args.output / 'cleavage_sites.csv', ['digest', 'enzyme', 'cut_boundary_0based', 'left_base_1based', 'right_base_1based'], cut_rows)
    write_csv(args.output / 'fragments.csv', ['digest', 'fragment_number', 'start_base_1based', 'end_base_1based', 'length_bp'], fragment_rows)
    plot_gel(gel_lanes, args.output / 'gel.png')
    print('Genome:', header, '| Length:', len(dna), 'bases')
    print('Enzymes (motif, cut offset):', ENZYMES)
    print(f'{"Digest":45} {"Cuts":>5} {"Fragments":>10}')
    for row in summary:
        print(f'{row[0]:45} {row[1]:5} {row[2]:10}')
    print('Results saved to', args.output.resolve())

if __name__ == '__main__':
    main()
