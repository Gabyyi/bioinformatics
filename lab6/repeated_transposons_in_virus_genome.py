
from collections import defaultdict

def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                sequence += line.upper()

    return sequence


def find_repeats(dna_sequence, min_window=6, max_window=10):
    repeats = []
    for window_size in range(min_window, max_window + 1):
        positions = defaultdict(list)
        for i in range(len(dna_sequence) - window_size + 1):
            window = dna_sequence[i:i + window_size]
            positions[window].append(i)

        for sequence, locations in positions.items():
            if len(locations) >= 2:
                repeats.append({
                    "sequence": sequence,
                    "length": window_size,
                    "positions": locations
                })

    return repeats


def main():
    filename = "GCF_000853865.1_ViralProj14960_genomic.fna"
    dna_sequence = read_fasta(filename)
    print("DNA sequence length:", len(dna_sequence))
    print("\nSearching for direct repeats...\n")
    repeats = find_repeats(dna_sequence)
    for repeat in repeats:
        print("Repeated sequence:", repeat["sequence"])
        print("Length:", repeat["length"])
        for start in repeat["positions"]:
            end = start + repeat["length"]
            print("Start:", start)
            print("End:", end)
        print("-----------------------------")
    print("Total repeated patterns found:", len(repeats))

if __name__ == "__main__":
    main()
