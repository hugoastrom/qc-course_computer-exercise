import pyscf
import numpy as np

if __name__ == "__main__":

    #for basis in ("sto-3g", "cc-pvdz", "aug-cc-pvdz", "cc-pvtz", "aug-cc-pvtz"):
    basis = "aug-cc-pvtz"
    with open(f"ccsd_{basis}.dat", "w") as f:
        for r in np.arange(0.4, 4.1, 0.1):
            # Create molecule object
            mol = pyscf.gto.Mole()
            mol.atom = f"H 0 0 0; H {r} 0 0"
            mol.basis = basis
            mol.spin = 0
            mol.build()
            
            print("basis = %s r = %.2f" %(basis, r))
            mf = mol.RHF().run()
            cc = pyscf.cc.CCSD(mf).run()
            print('E(CCSD) = %.12f' %cc.e_tot )
            f.write("%.2f %.12f\n" %(r, cc.e_tot))
