#
#  alm_opt.py
#
#  This is an example to run ALM in the fitting mode.
#
import h5py
from alm import ALM
import numpy as np

def write_fc_to_hdf5(alm, only_fc2=False):
    natom = len(alm.numbers)
    fc2 = np.zeros((natom, natom, 3, 3), dtype='double', order='C')
    for fc, indices in zip(*alm.get_fc(1, mode='all')):
        v1, v2 = indices // 3
        c1, c2 = indices % 3
        fc2[v1, v2, c1, c2] = fc

    if not only_fc2:
        fc3 = np.zeros((natom, natom, natom, 3, 3, 3),
                       dtype='double', order='C')
        for (fc, indices) in zip(*alm.get_fc(2, mode='all')):
            v1, v2, v3 = indices // 3
            c1, c2, c3 = indices % 3
            fc3[v1, v2, v3, c1, c2, c3] = fc
            fc3[v1, v3, v2, c1, c3, c2] = fc

    with h5py.File('fc.hdf5', 'w') as w:
        w.create_dataset('fc2', data=fc2, compression='gzip')
        if not only_fc2:
            w.create_dataset('fc3', data=fc3, compression='gzip')


def cv(alm):
    optcontrol = {'linear_model': 2,
                  'cross_validation': 4,
                  'num_l1_alpha': 50}
    alm.set_optimizer_control(optcontrol)
    alm.optimize()
    return alm.get_cv_l1_alpha()


def optimize(alm, cv_l1_alpha):
    optcontrol = {'linear_model': 2,
                  'cross_validation': 0,  # change 2 -> 0
                  'l1_alpha': cv_l1_alpha}
    alm.set_optimizer_control(optcontrol)
    alm.optimize()


#def main(filename_prefix, num_snapshots):
#    
#    #crystal, disp, force = get_inputs(num_snapshots)
#    cutoff_radii = [np.ones((2, 2)) * -1, np.ones((2, 2)) * 4.8]
#
#    with ALM(*crystal) as alm:
#        alm.set_verbosity(1)
#        # alm.set_output_filename_prefix(filename_prefix)
#        alm.displacements = disp
#        alm.forces = force
#        alm.define(2, cutoff_radii=cutoff_radii)
#        cv_l1_alpha = cv(alm)
#        optimize(alm, cv_l1_alpha)
#        write_fc_to_hdf5(alm)
#        print("Force constants were written to fc.hdf5.")
#
#
#if __name__ == '__main__':
#    main("AlN", 40)

