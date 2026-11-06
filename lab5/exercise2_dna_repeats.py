import urllib.request
import urllib.error

def fetch_dna_from_ncbi(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    )
    request = urllib.request.Request(
        url, headers={"User-Agent": "DNA-Repetition-Detector-Student-Project"}
    )

    with urllib.request.urlopen(request, timeout=20) as response:
        fasta = response.read().decode("utf-8")

    sequence = "".join(
        line.strip() for line in fasta.splitlines() if not line.startswith(">")
    ).upper()

    if not sequence or any(base not in "ATCGN" for base in sequence):
        raise ValueError("NCBI did not return a valid DNA sequence.")

    return sequence


def detect_tandem_repeats(sequence):
    repeats = []
    sequence_length = len(sequence)

    for unit_length in range(3, 6):
        position = 0

        while position <= sequence_length - (unit_length * 2):
            unit = sequence[position:position + unit_length]
            end = position + unit_length * 2

            if sequence[position + unit_length:end] == unit:
                copies = 2

                while (
                    end + unit_length <= sequence_length
                    and sequence[end:end + unit_length] == unit
                ):
                    copies += 1
                    end += unit_length

                repeats.append({
                    "start": position,
                    "end": end,
                    "unit": unit,
                    "copies": copies,
                    "length": end - position,
                    "sequence": sequence[position:end],
                })

                position = end
            else:
                position += 1

    repeats.sort(key=lambda item: (item["start"], len(item["unit"])))

    unique_repeats = []
    seen = set()

    for repeat in repeats:
        key = (repeat["start"], repeat["end"], repeat["unit"])
        if key not in seen:
            seen.add(key)
            unique_repeats.append(repeat)

    return unique_repeats


def main():
    accession = input("Enter an NCBI nucleotide accession: ").strip()

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
    repeats = detect_tandem_repeats(dna)

    if not repeats:
        print("No tandem repeats were found.")
        return

    print(f"\nFound {len(repeats)} repeated regions:\n")

    for number, repeat in enumerate(repeats, start=1):
        print(f"Repeat {number}")
        print(f"  Start:   {repeat['start']}")
        print(f"  End:     {repeat['end']} (exclusive)")
        print(f"  Unit:    {repeat['unit']}")
        print(f"  Copies:  {repeat['copies']}")
        print(f"  Length:  {repeat['length']} bases")
        print(f"  Region:  {repeat['sequence']}\n")


if __name__ == "__main__":
    main()
