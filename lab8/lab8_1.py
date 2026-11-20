# Implement a software app that alligns two DNA sequences by using Needdleman-Wunsch algorithm.

# S1 = "ACCGTGAAGCCAATAC"
# S2 = "AGCGTGCAGCCAATAC"

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap

    for j in range(1, m + 1):
        score[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = score[i-1][j-1] + (
                match if seq1[i-1] == seq2[j-1]
                else mismatch
            )
            up = score[i-1][j] + gap
            left = score[i][j-1] + gap

            score[i][j] = max(diagonal, up, left)

    aligned1 = ""
    aligned2 = ""
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diagonal_score = (
                match if seq1[i-1] == seq2[j-1]
                else mismatch
            )
            if score[i][j] == score[i-1][j-1] + diagonal_score:
                aligned1 = seq1[i-1] + aligned1
                aligned2 = seq2[j-1] + aligned2
                i -= 1
                j -= 1
                continue

        if i > 0 and score[i][j] == score[i-1][j] + gap:
            aligned1 = seq1[i-1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        else:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j-1] + aligned2
            j -= 1

    return aligned1, aligned2, score[n][m], score

def print_matrix(score, seq1, seq2):
    print("\nScoring Matrix:")
    print(f"{'':>5}{'-':>5}", end="")
    for char in seq2:
        print(f"{char:>5}", end="")
    print()
    for i, row in enumerate(score):
        label = "-" if i == 0 else seq1[i-1]
        print(f"{label:>5}", end="")

        for value in row:
            print(f"{value:>5}", end="")
        print()

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"
aligned1, aligned2, optimal_score, matrix = needleman_wunsch(S1, S2)
print("Needleman-Wunsch DNA Alignment")
print_matrix(matrix, S1, S2)
print("\nOptimal Alignment:")
print("Sequence 1:", aligned1)
print("            ", "".join(
    "|" if a == b else " "
    for a, b in zip(aligned1, aligned2)
))
print("Sequence 2:", aligned2)
print("\nOptimal Alignment Score:", optimal_score)
