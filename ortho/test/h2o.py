#!/usr/bin/env python

from pyscf import gto, scf
from pyscf.tools import wfn_format

name = 'h2o'

mol = gto.Mole()
mol.atom = '''
 O  0.0000  0.0000  0.1173
 H  0.0000  0.7572 -0.4692
 H  0.0000 -0.7572 -0.4692
'''
mol.basis = 'aug-cc-pv5z'
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.kernel()

wfn_file = name + '.wfn'
nmo = mol.nelectron//2
coef = mf.mo_coeff[:,:nmo]
occ = mf.mo_occ[:nmo]
ene = mf.mo_energy[:nmo]
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, coef, mo_occ=occ, mo_energy=ene)

