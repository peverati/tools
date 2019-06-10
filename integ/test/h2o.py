#!/usr/bin/env python

import avas, numpy
from pyscf import gto, scf, mcscf, fci, ao2mo, lib
from pyscf.tools import wfn_format

name = 'h2o'

mol = gto.Mole()
mol.atom = '''
 O  0.0000  0.0000  0.1173
 H  0.0000  0.7572 -0.4692
 H  0.0000 -0.7572 -0.4692
'''
mol.basis = 'cc-pvtz'
mol.verbose = 4
mol.spin = 0
mol.symmetry = 1
mol.charge = 0
mol.build()

mf = scf.RHF(mol)
mf.chkfile = name+'.chk'
#mf.__dict__.update(scf.chkfile.load(name+'.chk', 'scf'))
#dm = mf.make_rdm1()
mf.kernel()

ncore = 1
aolst1 = ['O 2s', 'O 2p']
aolst = aolst1 
ncas, nelecas, mo = avas.kernel(mf, aolst, threshold_occ=0.1, threshold_vir=0.01, minao='minao', ncore=ncore)

mc = mcscf.CASSCF(mf, ncas, nelecas)
mc.max_cycle_macro = 250
mc.max_cycle_micro = 7
mc.chkfile = name+'.chk'
mc.fcisolver = fci.direct_spin0_symm.FCI(mol)
mc.fix_spin_(shift=.5, ss=0)
#mc.__dict__.update(scf.chkfile.load(name+'.chk', 'mcscf'))
#mo = lib.chkfile.load(name+'.chk', 'mcscf/mo_coeff')
mc.kernel(mo)

nmo = mc.ncore + mc.ncas
rdm1, rdm2 = mc.fcisolver.make_rdm12(mc.ci, mc.ncas, mc.nelecas) 
rdm1, rdm2 = mcscf.addons._make_rdm12_on_mo(rdm1, rdm2, mc.ncore, mc.ncas, nmo)

eri_mo = ao2mo.kernel(mf._eri, mc.mo_coeff[:,:nmo], compact=False)
eri_mo = eri_mo.reshape(nmo,nmo,nmo,nmo)
h1 = reduce(numpy.dot, (mc.mo_coeff[:,:nmo].T, mf.get_hcore(), \
                        mc.mo_coeff[:,:nmo]))
e1 = numpy.einsum('ij,ij->', h1, rdm1)
lib.logger.info(mc,"h1e energy : %.8f" % e1)    
e2 = numpy.einsum('ijkl,ijkl->', eri_mo, rdm2)*0.5 
lib.logger.info(mc,"h2e energy : %.8f" % e2)    
enuc = mf.mol.energy_nuc()  
lib.logger.info(mc,"NN energy : %.8f" % enuc)    
etot = e1 + e2 + enuc
lib.logger.info(mc,"* Energy with 1/2-RDM : %.8f" % etot)    

wfn_file = name + '.wfn'
coef = mc.mo_coeff[:,:nmo]
occ = mf.mo_occ[:nmo]
ene = mf.mo_energy[:nmo]
with open(wfn_file, 'w') as f2:
    wfn_format.write_mo(f2, mol, coef, mo_occ=occ, mo_energy=ene)
    f2.write('CCIQA\n')
    f2.write('1-RDM:\n')
    for i in range(nmo):
        for j in range(nmo):
            f2.write('%i %i %.10f\n' % ((i+1), (j+1), rdm1[i,j]))
    f2.write('2-RDM:\n')
    for i in range(nmo):
        for j in range(nmo):
            for k in range(nmo):
                for l in range(nmo):
                    if (abs(rdm2[i,j,k,l]) > 1e-10):
                        f2.write('%i %i %i %i %.10f\n' % ((i+1), \
                        (j+1), (k+1), (l+1), rdm2[i,j,k,l]))

