# Design a software app that is able to read a fasta file of one genome from either a virus 
# or a bacteria.  Measure the frequencies of the codons from the sequence of the genome. 
# In order to measure the frequency use a sliding window of 3 positions to slide it over 
# the genome and counting the frequency of each codon in the sliding window. 
# The app must show the highest frequency top 10 codons, also your app must show for each 
# top 10 codons the associated aminoacid(from the table). 

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

RNA_CODON_TABLE = {
    #U Row
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "Stop", "UAG": "Stop",
    "UGU": "Cys", "UGC": "Cys", "UGA": "Stop", "UGG": "Trp",
    
    #C Row
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    
    #A Row
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    
    #G Row
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

AMINO_ACID_TO_PROTEIN = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

DNA_CODON_TABLE = {codon.replace("U", "T"): amino_acid for codon, amino_acid in RNA_CODON_TABLE.items()}


def read_fasta(file_path: str | Path) -> tuple[str, str]:
    """Read a FASTA file and return the header and the joined sequence."""
    fasta_file = Path(file_path)
    if not fasta_file.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")

    header = ""
    sequence_parts: list[str] = []

    with fasta_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:].strip()
                continue
            sequence_parts.append(line.upper())

    if not header:
        raise ValueError(f"The file {fasta_file} is not a valid FASTA file.")

    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(f"The FASTA file {fasta_file} does not contain any sequence data.")

    return header, sequence


def count_codons(sequence: str) -> Counter[str]:
    """Count codon occurrences in a sliding 3-base window across the sequence."""
    normalized = "".join(base for base in sequence.upper() if base in {"A", "C", "G", "T", "U"})
    normalized = normalized.replace("T", "U")

    codon_counts: Counter[str] = Counter()
    for index in range(len(normalized) - 2):
        codon = normalized[index : index + 3]
        if len(codon) == 3 and set(codon) <= {"A", "C", "G", "U"}:
            codon_counts[codon] += 1

    return codon_counts


def get_top_codons(counts: dict[str, int] | Counter[str], top_n: int = 10) -> list[tuple[str, int]]:
    """Return the codons with the highest counts, sorted from most common to least common."""
    if not counts:
        return []

    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:top_n]


def get_amino_acid(codon: str) -> str:
    """Return the amino acid name associated with a codon from the standard genetic code table."""
    return RNA_CODON_TABLE.get(codon.upper(), "Unknown")


def calculate_top_codon_report(sequence: str, top_n: int = 10) -> list[dict[str, object]]:
    """Build a ranked codon report with count, frequency, and amino acid."""
    codon_counts = count_codons(sequence)
    total_codons = sum(codon_counts.values())
    top_codons = get_top_codons(codon_counts, top_n)

    report: list[dict[str, object]] = []
    for rank, (codon, count) in enumerate(top_codons, start=1):
        percentage = (count / total_codons * 100) if total_codons else 0.0
        report.append(
            {
                "rank": rank,
                "codon": codon,
                "count": count,
                "frequency_percent": round(percentage, 3),
                "amino_acid": get_amino_acid(codon),
            }
        )

    return report


def main() -> None:
    default_fasta = Path(__file__).with_name("GCF_000853865.1_ViralProj14960_genomic.fna")
    parser = argparse.ArgumentParser(description="Find the top codons in a bacterial or viral genome FASTA file.")
    parser.add_argument(
        "fasta_file",
        nargs="?",
        default=str(default_fasta),
        help="Path to the input FASTA file. Defaults to the coronavirus genome in this folder.",
    )
    parser.add_argument("--top", type=int, default=10, help="Number of codons to show in the report.")
    args = parser.parse_args()

    fasta_path = Path(args.fasta_file)
    header, sequence = read_fasta(fasta_path)
    top_report = calculate_top_codon_report(sequence, args.top)

    print(f"Genome name: {header}")
    print(f"Genome length: {len(sequence)} bases")
    print(f"Valid codons analyzed: {sum(item['count'] for item in top_report)}")
    print("\nTop codons by frequency (sliding window of 3 bases):")
    print("Rank | Codon | Count | Frequency % | Amino acid")
    print("-----|-------|-------|------------|-----------")

    for item in top_report:
        codon = str(item["codon"]).replace("T", "U")
        print(
            f"{item['rank']:>4} | {codon:<6} | {item['count']:>5} | {item['frequency_percent']:>10.3f} | {item['amino_acid']}"
        )


if __name__ == "__main__":
    main()
