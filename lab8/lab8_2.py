# Implement a software app that moks DNA pattern from DNA seq of DNA parameters
# Compute the (C+G)% and the information content values from each sliding window(SW)
# to ensure that the alg. implementation works correctly use the following steps:
# 1. Use S="CGGACTGATCTATCT[A(27)]CGTAGCATCTATCGATCTATCTAGCGATCTA"
# 2. Use SW=30
# 3. Build a function to process C+G content - the value that f. must return is CG~29.27
# 4. Build a function to process information content - the value that f. must return is IC~0.29
# 5. Plot the pattern on a chart
# 6. Compute the center of weight of the pattern
# 7. Take the center of each pattern and plot it on a second chart
# 8. Take the DNA sequence of a promoter and generate a pattern

import math
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
SW = 30

def gc_content(sequence):
    gc = sequence.count("G") + sequence.count("C")
    return (gc / len(sequence)) * 100 if sequence else 0

def information_content(sequence):
    if not sequence:
        return 0

    entropy = 0
    for nucleotide in "ACGT":
        p = sequence.count(nucleotide) / len(sequence)
        if p > 0:
            entropy -= p * math.log2(p)
    return 2 - entropy

def analyze_sequence(sequence, sw):
    sequence = sequence.upper().replace(" ", "")
    if not sequence or any(c not in "ACGT" for c in sequence):
        raise ValueError("Sequence must contain only A, C, G, T")

    if sw < 1 or sw > len(sequence):
        raise ValueError("Invalid sliding window size")

    positions = []
    gc_values = []
    ic_values = []
    for i in range(len(sequence) - sw + 1):
        window = sequence[i:i + sw]
        center = i + (sw + 1) / 2
        positions.append(center)
        gc_values.append(gc_content(window))
        ic_values.append(information_content(window))
    return positions, gc_values, ic_values

def center_of_weight(positions, weights):
    total_weight = sum(weights)
    if total_weight == 0:
        return None
    return sum(
        x * w for x, w in zip(positions, weights)
    ) / total_weight

def find_patterns(positions, values):
    peaks = []
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and abs(values[j + 1] - values[i]) < 1e-10:
            j += 1

        left = values[i - 1] if i > 0 else float("-inf")
        right = values[j + 1] if j + 1 < len(values) else float("-inf")
        if values[i] > left and values[j] > right:
            peak_center = (positions[i] + positions[j]) / 2
            peaks.append((peak_center, values[i]))

        i = j + 1

    return peaks

def plot_patterns(positions, gc_values, ic_values, peaks, cow):
    fig, axes = plt.subplots(2, 1, figsize=(12, 9))
    ax1 = axes[0]
    ax1.plot(
        positions, gc_values,
        color="green", label="GC Content (%)"
    )
    ax1.set_xlabel("DNA Position")
    ax1.set_ylabel("GC Content (%)", color="green")
    ax1.tick_params(axis="y", labelcolor="green")
    ax1b = ax1.twinx()
    ax1b.plot(
        positions, ic_values,
        color="blue", label="Information Content"
    )
    ax1b.set_ylabel("Information Content (bits)", color="blue")
    ax1b.tick_params(axis="y", labelcolor="blue")
    ax1.set_title("DNA Sliding Window Pattern")
    ax1.grid(True, alpha=0.3)

    ax2 = axes[1]
    ax2.plot(
        positions, ic_values,
        color="blue", alpha=0.5, label="IC Pattern"
    )
    for center, value in peaks:
        ax2.scatter(center, value, color="red", zorder=5)

    if peaks:
        ax2.scatter([], [], color="red", label="Pattern Centers")

    if cow is not None:
        ax2.axvline(
            cow, color="purple", linestyle="--",
            label=f"Center of Weight = {cow:.2f}"
        )

    ax2.set_title("Detected Pattern Centers")
    ax2.set_xlabel("DNA Position")
    ax2.set_ylabel("Information Content (bits)")
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def run_analysis(sequence, sw=30):
    positions, gc_values, ic_values = analyze_sequence(sequence, sw)
    cow = center_of_weight(positions, ic_values)
    peaks = find_patterns(positions, ic_values)
    print("\nDNA PATTERN ANALYSIS")
    print("-" * 55)
    print("Sequence length:", len(sequence))
    print("Sliding window:", sw)
    print("Number of windows:", len(positions))
    print(f"\nOverall GC Content: {gc_content(sequence):.2f}%")
    print(f"Overall IC: {information_content(sequence):.4f}")
    print("\nSliding Window Results:")
    print(f"{'Window':<10}{'Center':<12}{'GC (%)':<12}{'IC':<12}")
    for i, (pos, gc, ic) in enumerate(
        zip(positions, gc_values, ic_values), start=1
    ):
        print(f"{i:<10}{pos:<12.1f}{gc:<12.2f}{ic:<12.4f}")

    print("\nCenter of Weight:",
          f"{cow:.2f}" if cow is not None else "Undefined")

    print("\nDetected Pattern Centers:")
    for center, value in peaks:
        print(f"Position: {center:.2f}, IC: {value:.4f}")

    plot_patterns(positions, gc_values, ic_values, peaks, cow)


if __name__ == "__main__":
    run_analysis(S, SW)
    answer = input("\nAnalyze a promoter sequence? (y/n): ")

    if answer.lower() == "y":
        promoter = input("Enter promoter DNA sequence: ").strip().upper()
        run_analysis(promoter, SW)
