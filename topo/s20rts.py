import spharm_io
import matplotlib.pyplot as plt
from rotate_models import filter_sh
from plot import display_cartopy, display_subplots_cartopy
import os

def modify_model(path):
    coeffs_s20rts = spharm_io.read_sh(s20rts_path)
    coeffs_s20rts *= -1
    spharm_io.write_sh(path, coeffs_s20rts)

s20rts_path = 'resources/models/numpy/s20rts.npy'
figure_root = '../figures'

l_set = [2]

#modify_model(s20rts_path)

# s20rts
_, coeffs_s20rts = spharm_io.read_sh(s20rts_path)
_, coeffs_s20rts_filt = filter_sh(coeffs_s20rts, l_set)

coeffs_dict = {'S20RTS (2891 km)': coeffs_s20rts,
               'S20RTS (L=2)': coeffs_s20rts_filt}

# FIXME bug when displaying model with l=2
display_subplots_cartopy(coeffs_dict, nrows=2, ncols=1)

figure_path = os.path.join(figure_root, 's20rts_cmb.pdf')
plt.savefig(figure_path)
