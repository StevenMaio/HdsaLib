import re
import numpy as np
import matplotlib.pyplot as plt



def parse_solver_data(filename):
    """
    Read the given text file and return four NumPy vectors:
      1. inverse_iterations
      2. inverse_seconds
      3. sqrt_iterations
      4. sqrt_seconds
    """

    inverse_prefix = "W_u_acute_plus_scalar_M_u_Inverse::Iterative_Linear_Solve"
    sqrt_prefix = "W_u_acute_plus_scalar_M_u_Sqrt::Matrix_Sqrt_Apply"

    pattern = re.compile(
        r"with\s+(\d+)\s+iterations\s+in\s+([0-9.eE+-]+)\s+seconds"
    )

    inverse_iterations = []
    inverse_seconds = []
    sqrt_iterations = []
    sqrt_seconds = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith(inverse_prefix):
                match = pattern.search(line)
                if match:
                    inverse_iterations.append(int(match.group(1)))
                    inverse_seconds.append(float(match.group(2)))

            elif line.startswith(sqrt_prefix):
                match = pattern.search(line)
                if match:
                    sqrt_iterations.append(int(match.group(1)))
                    sqrt_seconds.append(float(match.group(2)))

    inverse_iterations = np.array(inverse_iterations, dtype=int)
    inverse_seconds = np.array(inverse_seconds, dtype=float)
    sqrt_iterations = np.array(sqrt_iterations, dtype=int)
    sqrt_seconds = np.array(sqrt_seconds, dtype=float)

    return inverse_iterations, inverse_seconds, sqrt_iterations, sqrt_seconds


# Example usage:
filename = "out_bilaplacian.txt"
inverse_iterations_bilaplacian, inverse_seconds_bilaplacian, sqrt_iterations_bilaplacian, sqrt_seconds_bilaplacian = parse_solver_data(filename)

filename = "out_lumped_mass.txt"
inverse_iterations_lumped_mass, inverse_seconds_lumped_mass, sqrt_iterations_lumped_mass, sqrt_seconds_lumped_mass = parse_solver_data(filename)

# Create scatter plot
plt.figure(figsize=(8, 6))

# Colorblind-friendly colors
color_bilaplacian = "#0072B2"   # blue
color_lumped_mass = "#D55E00"   # vermillion

label1 = r"$W+\alpha M$"
label2 = r"$\overline{W}+\alpha M$"

# Bilaplacian data
plt.scatter(
    inverse_iterations_bilaplacian,
    inverse_seconds_bilaplacian,
    facecolors="none",
    edgecolors=color_bilaplacian,
    marker="o",
    s=70,
    label=label1 + " Inverse"
)

plt.scatter(
    sqrt_iterations_bilaplacian,
    sqrt_seconds_bilaplacian,
    color=color_bilaplacian,
    marker="x",
    s=70,
    label=label1 + " Sqrt"
)

# Lumped mass data
plt.scatter(
    inverse_iterations_lumped_mass,
    inverse_seconds_lumped_mass,
    facecolors="none",
    edgecolors=color_lumped_mass,
    marker="o",
    s=70,
    label=label2 + " Inverse"
)

plt.scatter(
    sqrt_iterations_lumped_mass,
    sqrt_seconds_lumped_mass,
    color=color_lumped_mass,
    marker="x",
    s=70,
    label=label2 + " Sqrt"
)

# Labels and styling
plt.xlabel("Iterations", fontsize=16)
plt.ylabel("Seconds", fontsize=16)
plt.title("Computational Cost Summary", fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
#plt.grid(True, alpha=0.3)
plt.tight_layout()

# Write plot to EPS file
plt.savefig("Performance_Results.eps", format="eps")

#plt.show()