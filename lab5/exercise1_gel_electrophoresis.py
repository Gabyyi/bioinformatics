import random
import urllib.request
import urllib.error
import matplotlib.pyplot as plt


def fetch_dna_from_ncbi(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )
    request = urllib.request.Request(
        url, headers={"User-Agent": "DNA-Gel-Electrophoresis-Student-Project"}
    )

    with urllib.request.urlopen(request, timeout=20) as response:
        fasta = response.read().decode("utf-8")

    sequence = "".join(
        line.strip() for line in fasta.splitlines() if not line.startswith(">")
    ).upper()

    if not sequence or any(base not in "ATCGN" for base in sequence):
        raise ValueError("NCBI did not return a valid DNA sequence.")

    return sequence


def take_random_samples(sequence, number_of_samples=10):
    samples = []

    for _ in range(number_of_samples):
        fragment_length = random.randint(100, min(3000, len(sequence)))
        start = random.randint(0, len(sequence) - fragment_length)
        fragment = sequence[start:start + fragment_length]
        samples.append((start, start + fragment_length, fragment))

    return samples


def simulate_gel(samples):
    lengths = [len(fragment) for _, _, fragment in samples]
    max_length = max(lengths)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.set_facecolor("#17212b")

    for lane, length in enumerate(lengths, start=1):
        distance = 1.0 + 7.0 * (
            (max_length - length) / max_length
        )

        ax.plot(
            [lane - 0.32, lane + 0.32],
            [distance, distance],
            linewidth=5,
            color="#79e0c2",
            solid_capstyle="round",
        )
        ax.text(
            lane, 0.55, f"{length} bp",
            ha="center", va="bottom", color="white", fontsize=8
        )

    ax.set_xlim(0.5, len(samples) + 0.5)
    ax.set_ylim(0, 9)
    ax.invert_yaxis()
    ax.set_xticks(range(1, len(samples) + 1))
    ax.set_xticklabels([f"Sample {i}" for i in range(1, len(samples) + 1)])
    ax.set_ylabel("Migration distance (simulated)")
    ax.set_title("Simulated DNA Gel Electrophoresis", color="black")
    ax.grid(axis="y", alpha=0.2)
    plt.tight_layout()
    plt.show()


def main():
    accession = input("Enter an NCBI nucleotide accession (e.g. NC_000913): ").strip()

    try:
        dna = fetch_dna_from_ncbi(accession)
    except (urllib.error.URLError, TimeoutError, ValueError) as error:
        print("Could not retrieve DNA sequence:", error)
        return

    if not 1000 <= len(dna) <= 3000:
        print(f"Retrieved sequence has {len(dna)} bases.")
        print("The lab requires a sequence between 1000 and 3000 bases.")
        print("Choose an accession whose sequence is within that range.")
        return

    print(f"\nDNA sequence length: {len(dna)} bases")
    samples = take_random_samples(dna)

    print("\n10 randomly selected DNA fragments:")
    for number, (start, end, fragment) in enumerate(samples, start=1):
        print(
            f"Sample {number}: start={start}, end={end}, "
            f"length={len(fragment)} bp"
        )
        print(fragment)

    simulate_gel(samples)


if __name__ == "__main__":
    main()
