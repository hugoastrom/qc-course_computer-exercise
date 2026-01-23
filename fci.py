import pyscf
import numpy as np

if __name__ == "__main__":

    #for basis in ("sto-3g", "cc-pvdz", "aug-cc-pvdz", "cc-pvtz", "aug-cc-pvtz"):
    basis = "cc-pvdz"
    with open(f"fci_{basis}.dat", "w") as f:
        for r in np.arange(0.4, 4.1, 0.1):
            # Create molecule object
            mol = pyscf.gto.Mole()
            mol.atom = f"H 0 0 0; H {r} 0 0"
            mol.basis = basis
            mol.spin = 0
            mol.build()
            
            print("basis = %s r= %.2f" %(basis, r))
            rhf = mol.RHF().run()
            cisolver = pyscf.fci.FCI(rhf)
            e_tot = cisolver.kernel()[0]
            print('E(UHF-FCI) = %.12f' %e_tot )
            f.write("%.2f %.12f\n" %(r, e_tot))
