import matplotlib.pyplot as plt

# Figure font parameters
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 11,
    'legend.fontsize': 11,
    'font.family': 'lmodern',
    'text.usetex': True,
    'text.latex.preamble': (
        r'\usepackage{lmodern}'
    )
})

def plot(method, basis, label):

    radii = []
    energies = []
    try:
        with open(f"{method}_{basis}.dat") as f:
            for line in f:
                radii.append(float(line.split()[0]))
                energies.append(float(line.split()[1]))
    except:
        raise ValueError(f"File {method}_{basis}.dat does not exist")

    plt.plot(radii, energies, markers[methods.index(method)], color=colors[bases.index(basis)], fillstyle="none", linewidth=0.0, label=label)


# Plot parameters
markers = ["s", "o", "^", "v", "d"]
colors = ["blue", "crimson", "tab:orange", "tab:purple", "tab:pink", "tab:brown", "gold", "cornflowerblue", "navy", "rosybrown", "palegoldenrod"]

# Calculation inputs
methods = ["rhf", "mp2", "ccsd", "ccsdt", "fci"]
bases = ["sto-3g", "cc-pvdz", "cc-pvtz", "aug-cc-pvdz", "aug-cc-pvtz"]

# Legends
method_labels = {"rhf": "RHF", "mp2": "MP2", "ccsd": "CCSD", "ccsdt": "CCSD(T)", "fci": "FCI"}
basis_labels = {"sto-3g": "STO-3G", "cc-pvdz": "cc-pVDZ", "aug-cc-pvdz": "aug-cc-pVDZ", "cc-pvtz": "cc-pVTZ", "aug-cc-pvtz": "aug-cc-pVTZ"}

# Plot RHF energies for different basis sets
fig = plt.figure()
for basis in bases:
    plot("rhf", basis, f"{basis_labels[basis]}")
plt.legend()
plt.xlabel("r/Angstrom")
plt.ylabel("E/Hartree")
plt.savefig("rhf.pdf")
plt.close()

for method in methods:
    plot(method, "cc-pvdz", f"{method_labels[method]}")
plt.legend()
plt.xlabel("r/Angstrom")
plt.ylabel("E/Hartree")
plt.savefig("methods.pdf")
plt.close()
