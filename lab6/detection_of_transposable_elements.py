# 1. make an artificial DNA sequence of 200 up to 400 bases in length in which to simulate 3 to 4 
# transposable elements
# 2. implement a software application to detect the positionsof these transposable elements(start, end) 
# within the created DNA sequence

        #            /////TRANSPOSON/////
        #         ^                       ^
        #         |                       |
        #         ----INVERTED REPEATS-----
        # MOST DNA .... TAG ATGCA CAG....
        # TARGET SITE...ATG TACGT GTC...
        #                     ||
        # TAC ATGCA       TRANSPOSON      ........CAC...
        # ATG.......                      TACGT GTC...
        #                     ||
        #                 TRANSPOSON
        # TAC ATGCA                       ATGCA CAG...
        # TAG ATCGA                       TACGT GTC...
        #       ^                           ^
        #       |                           |
        #       --------DIRECT REPEATS-------



import random

BASES = "ATCG"

LEFT_REPEAT = "AGTCGATA"
RIGHT_REPEAT = "TATCGACT"

MIN_BODY_LENGTH = 20
MAX_BODY_LENGTH = 35

MIN_SEQUENCE_LENGTH = 200
MAX_SEQUENCE_LENGTH = 400


def generate_random_dna(length):
    return "".join(random.choice(BASES) for _ in range(length))


def generate_transposable_element():
    body_length = random.randint(
        MIN_BODY_LENGTH, MAX_BODY_LENGTH
    )

    body = generate_random_dna(body_length)
    element = LEFT_REPEAT + body + RIGHT_REPEAT
    return element


def generate_dna_sequence():
    number_of_elements = random.randint(3, 4)
    elements = [
        generate_transposable_element()
        for _ in range(number_of_elements)
    ]

    total_length = random.randint(
        MIN_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH
    )

    elements_length = sum(len(element) for element in elements)
    background_length = total_length - elements_length
    cuts = sorted(
        random.sample(
            range(background_length + number_of_elements),
            number_of_elements
        )
    )

    background_lengths = []
    previous = 0
    for cut in cuts:
        background_lengths.append(cut - previous)
        previous = cut

    background_lengths.append(
        background_length - sum(background_lengths)
    )

    dna_sequence = ""
    for i in range(number_of_elements):
        dna_sequence += generate_random_dna(background_lengths[i])
        dna_sequence += elements[i]

    dna_sequence += generate_random_dna(background_lengths[-1])
    return dna_sequence


def detect_transposable_elements(dna_sequence):
    detected_elements = []
    element_min_length = (
        len(LEFT_REPEAT) + MIN_BODY_LENGTH + len(RIGHT_REPEAT)
    )

    element_max_length = (
        len(LEFT_REPEAT) + MAX_BODY_LENGTH + len(RIGHT_REPEAT)
    )

    i = 0
    while i <= len(dna_sequence) - element_min_length:
        if dna_sequence[i:i + len(LEFT_REPEAT)] == LEFT_REPEAT:
            for element_length in range(element_min_length, element_max_length + 1):
                end = i + element_length
                if end > len(dna_sequence):
                    break

                right_start = end - len(RIGHT_REPEAT)
                right_repeat = dna_sequence[
                    right_start:end
                ]

                if right_repeat == RIGHT_REPEAT:
                    detected_elements.append({
                        "start": i,
                        "end": end,
                        "length": element_length,
                        "sequence": dna_sequence[i:end]
                    })
                    i = end
                    break

        i += 1

    return detected_elements


if __name__ == "__main__":
    dna_sequence = generate_dna_sequence()
    print("DNA SEQUENCE:")
    print(dna_sequence)
    print("\nDNA sequence length:", len(dna_sequence))
    detected_elements = detect_transposable_elements(
        dna_sequence
    )
    print("\nDETECTED TRANSPOSABLE ELEMENTS:")
    for index, element in enumerate(detected_elements, start=1):
        print(f"\nTransposable Element {index}")
        print("Start:", element["start"])
        print("End:", element["end"])
        print("Length:", element["length"])
        print("Sequence:", element["sequence"])

    print("\nTotal elements detected:", len(detected_elements))
