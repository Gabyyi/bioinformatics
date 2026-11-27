# Design a software application which considers 2 inputs: one fasta file containing bacterial genome and 
# a set of examples that are used in a machine learning to detect novel functional structures over the genome. 
# The set of examples may be either sequences from the exon-intorn border or intron-exon border.
# Your application must show a chart representing the scores for each sliding window over the genome.

import numpy as np
import matplotlib.pyplot as plt

FASTA_FILE = "GCF_000853865.1_ViralProj14960_genomic.fna"

NUCLEOTIDES = ["A", "T", "C", "G"]

SW = 9
STEP = 1
NULL_MODEL = 0.25
PSEUDOCOUNT = 0.01

# Known motif sequences from the blackboard
EXAMPLES = [
    "GAGGTAAAC",
    "TCCGTAAGT",
    "CAGGTTGGA",
    "ACAGTCAGT",
    "TAGGTCATT",
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC",
    "TGTGTGAGT",
    "AAGGTAAGT"
]


# 1. Read genome from FASTA file
def read_fasta(filename):
    genome = ""

    with open(filename, "r") as file:
        for line in file:
            if not line.startswith(">"):
                genome += line.strip().upper()

    return genome


# 2. Build Position Frequency Matrix
def build_pfm(examples):
    pfm = np.zeros((4, SW), dtype=int)

    for sequence in examples:
        for position, nucleotide in enumerate(sequence):
            row = NUCLEOTIDES.index(nucleotide)
            pfm[row, position] += 1

    return pfm


# 3. Build Position Weight Matrix
def build_pwm(pfm, number_of_examples):
    pwm = pfm / number_of_examples

    return pwm


# 4. Build Log-Likelihood Matrix
def build_llm(pfm, number_of_examples):
    # Apply pseudocount to avoid log(0)
    probabilities = (
        (pfm + PSEUDOCOUNT) /
        (number_of_examples + 4 * PSEUDOCOUNT)
    )

    llm = np.log(probabilities / NULL_MODEL)

    return llm


# 5. Calculate the score of one sliding window
def calculate_score(window, llm):
    score = 0.0

    for position, nucleotide in enumerate(window):

        if nucleotide not in NUCLEOTIDES:
            return np.nan

        row = NUCLEOTIDES.index(nucleotide)

        score += llm[row, position]

    return score


# 6. Scan the entire genome
def scan_genome(genome, llm):
    positions = []
    scores = []

    for i in range(0, len(genome) - SW + 1, STEP):

        window = genome[i:i + SW]

        score = calculate_score(window, llm)

        positions.append(i + 1)
        scores.append(score)

    return np.array(positions), np.array(scores)


# 7. Print a matrix
def print_matrix(name, matrix):
    print(f"\n{name}:")

    print("     ", end="")

    for i in range(SW):
        print(f"{i + 1:8}", end="")

    print()

    for nucleotide, row in zip(NUCLEOTIDES, matrix):
        print(f"{nucleotide:5}", end="")

        for value in row:
            print(f"{value:8.2f}", end="")

        print()


# 8. Display the highest-scoring sequences
def display_top_matches(genome, positions, scores):
    valid = np.where(np.isfinite(scores))[0]

    top_indices = valid[
        np.argsort(scores[valid])[::-1][:10]
    ]

    print("\nTop 10 potential functional structures:")

    for index in top_indices:
        position = positions[index]
        score = scores[index]

        sequence = genome[position - 1:position - 1 + SW]

        print(
            f"Position: {position:6} | "
            f"Sequence: {sequence} | "
            f"Score: {score:.4f}"
        )


# 9. Plot the sliding window scores
def plot_scores(positions, scores):
    plt.figure(figsize=(16, 6))

    plt.plot(
        positions,
        scores,
        color="blue",
        linewidth=0.6,
        label="Sliding window score"
    )

    plt.axhline(
        y=0,
        color="red",
        linestyle="--",
        label="Null model reference"
    )

    plt.title("Functional Structure Detection - PWM / LLM")
    plt.xlabel("Genome position")
    plt.ylabel("Log-likelihood score")

    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


# Main application
def main():
    genome = read_fasta(FASTA_FILE)

    print("Genome length:", len(genome))
    print("Number of examples:", len(EXAMPLES))
    print("Sliding window size:", SW)

    # Learn the motif from known examples
    pfm = build_pfm(EXAMPLES)
    pwm = build_pwm(pfm, len(EXAMPLES))
    llm = build_llm(pfm, len(EXAMPLES))

    # Display all three matrices
    print_matrix("PFM", pfm)
    print_matrix("PWM", pwm)
    print_matrix("LLM", llm)

    # Scan genome
    positions, scores = scan_genome(genome, llm)

    print("\nTotal sliding windows:", len(scores))

    display_top_matches(genome, positions, scores)

    # Generate chart
    plot_scores(positions, scores)


if __name__ == "__main__":
    main()