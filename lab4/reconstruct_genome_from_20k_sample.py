# split the genome into subsequences such that 20k samples are taken about the genome
# each sample must have a lenght betwenn 100 and 150 letters

# design an application that is able to use the 20k samples in order to reconstruct the original genome
# the case here is that you consider that you dont have access to the original genome

import csv
import random
from collections import defaultdict


def read_fasta_file(filepath: str) -> str:
    seq_chunks = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_chunks.append(line.upper())
    return "".join(seq_chunks)


def generate_and_save_samples(
    genome: str,
    output_csv: str = "coronavirus_20k_samples.csv",
    output_fasta: str = "coronavirus_20k_samples.fasta",
    num_samples: int = 20_000,
    min_len: int = 100,
    max_len: int = 150,
    seed: int = 42,
):
    random.seed(seed)
    genome_len = len(genome)

    with open(output_csv, "w", newline="", encoding="utf-8") as f_csv, \
         open(output_fasta, "w", encoding="utf-8") as f_fa:
        
        writer = csv.writer(f_csv)
        writer.writerow(["sample_id", "start_pos_0based", "end_pos", "length", "subsequence"])

        for i in range(1, num_samples + 1):
            s_len = random.randint(min_len, max_len)
            start = random.randint(0, genome_len - s_len)
            end = start + s_len
            subseq = genome[start:end]
            writer.writerow([i, start, end, s_len, subseq])
            f_fa.write(f">sample_{i} length={s_len} start={start} end={end}\n{subseq}\n")

    print(f"Successfully generated {num_samples:,} samples.")


def load_samples_from_csv(csv_path: str) -> list[str]:
    reads = []
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            reads.append(row["subsequence"])
    return reads


def assemble_reads_debruijn(reads: list[str], k: int = 41) -> str:
    kmers = set()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmers.add(read[i : i + k])

    graph = defaultdict(set)
    in_degrees = defaultdict(int)

    for kmer in kmers:
        u = kmer[:-1]
        v = kmer[1:]
        graph[u].add(v)
        in_degrees[v] += 1

    start_nodes = [u for u in graph if in_degrees[u] == 0]
    if not start_nodes:
        start_nodes = [next(iter(graph.keys()))]

    curr = start_nodes[0]
    assembled_seq = [curr]
    visited = set()

    # Traverse graph path to the sink node
    while curr in graph and curr not in visited:
        visited.add(curr)
        next_nodes = list(graph[curr])
        # Follow the unique path forward
        nxt = next_nodes[0]
        assembled_seq.append(nxt[-1])
        curr = nxt

    return "".join(assembled_seq)


if __name__ == "__main__":
    fasta_input = "GCF_000853865.1_ViralProj14960_genomic.fna"
    csv_file = "coronavirus_20k_samples.csv"
    fasta_output = "coronavirus_20k_samples.fasta"

    orig_genome = read_fasta_file(fasta_input)
    generate_and_save_samples(
        genome=orig_genome,
        output_csv=csv_file,
        output_fasta=fasta_output,
        num_samples=20_000,
        min_len=100,
        max_len=150,
        seed=42
    )

    input_reads = load_samples_from_csv(csv_file)
    print(f"\nLoaded {len(input_reads):,} read sequences for blind assembly.")

    assembled_genome = assemble_reads_debruijn(input_reads, k=41)

    print(f"\nAssembly Results:")
    print(f"  - Reconstructed Contig Length : {len(assembled_genome):,} bp")
    print(f"  - Original Genome Length      : {len(orig_genome):,} bp")
    
    is_exact_subsequence = assembled_genome in orig_genome
    print(f"  - 100% Sequence Accuracy Match: {is_exact_subsequence}")
    if is_exact_subsequence:
        offset = orig_genome.find(assembled_genome)
        print(f"  - Genome Span Recovered       : positions {offset} to {offset + len(assembled_genome)}")
