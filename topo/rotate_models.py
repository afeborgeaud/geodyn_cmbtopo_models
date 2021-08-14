"""
2021.08.13
Anselme F.E. Borgeaud

Rotate CMB topography models specified by spherical harmonic coefficients
so  that their degree-2 pattern best-fit the degree-2 pattern of
S-velocities anomalies at the CMB from the global model S20RTS
"""

import pyshtools
import numpy as np
import os
os.environ['PROJ_LIB'] = r'/Users/navy/anaconda3/envs/ev02/share/proj'
import spharm_io
import matplotlib.pyplot as plt
import logging
from plot import display_subplots_cartopy, parse_title

logging.basicConfig(level=logging.INFO)


def find_best_angles(coeffs_r, coeffs_r_target, phis=np.arange(0., 360., 2.),
    thetas=np.arange(0., 180., 2.)):
    xx, yy = np.meshgrid(phis, thetas)
    angles = np.dstack((xx, yy)).reshape((len(phis)*len(thetas), 2))
    euler_angles = np.column_stack((angles, np.zeros(len(angles))))
    euler_angles = euler_angles * np.pi / 180.
    errors = np.apply_along_axis(rotate_and_compute_error, 1, euler_angles,
        coeffs_r, coeffs_r_target)
    imin = errors.argmin()
    best_angles = euler_angles[imin]
    return best_angles, errors


def error_function(y_pred, y):
    return .5 * (1. - np.corrcoef(y_pred.ravel(), y.ravel())[0,1])


def rotate_and_compute_error(euler_angles, coeffs_r, coeffs_r_target):
    angles = np.array(euler_angles)
    dj = pyshtools.rotate.djpi2(coeffs_r.shape[1]-1)
    coeffs_r_rot = pyshtools.rotate.SHRotateRealCoef(coeffs_r, angles, dj)
    error = error_function(coeffs_r_rot, coeffs_r_target)
    return error


def filter_sh(coeffs, l):
    lmax = coeffs.shape[1] - 1
    mask = np.isin(list(range(lmax+1)), l)
    mask = ~mask
    coeffs_filt = coeffs.copy()
    coeffs_filt[:, mask, :] = 0.
    return coeffs_filt


def truncate_sh(coeffs, lmax):
    coeffs_tr = np.zeros((2, lmax+1, lmax+1))
    coeffs_tr[...] = coeffs[:, :lmax+1, :lmax+1]
    return coeffs_tr


if __name__ == '__main__':
    model_root = 'resources/models/dat'
    models = [
        'CMBtopo-T1-pPv_eta1e-3_SHcoef.dat',
        'CMBtopo-T1-pPv_SHcoef.dat',
        'CMBtopo-T1_SHcoef.dat',
        'CMBtopo-TC1-pPv_eta1e-3_SHcoef.dat',
        'CMBtopo-TC1-pPv_SHcoef.dat',
        'CMBtopo-TC1_SHcoef.dat',
        'CMBtopo-TC4_SHcoef.dat',
        'CMBtopo-TC6_SHcoef.dat',
        'CMBtopo-TC7_SHcoef.dat'
    ]
    s20rts_path = 'resources/models/numpy/s20rts.npy'
    # output folder for figures
    figure_root = '../figures'

    lmax_filt = 2
    l_set = [2]

    # s20rts
    lmax_s20rts, coeffs_s20rts = spharm_io.read_sh(
        os.path.expanduser(s20rts_path))
    coeffs_s20rts_filt = truncate_sh(coeffs_s20rts, lmax_filt)
    coeffs_s20rts_filt = filter_sh(coeffs_s20rts_filt, l_set)

    best_angles = dict()
    coeffs_r_filt_best = dict()
    coeffs_r_best = dict()

    for model in models:
        print(model)
        model_path = os.path.join(
            os.path.expanduser(model_root), model)

        lmax, coeffs_r = spharm_io.read_sh_c_to_r(model_path)

        # filter to l=0–20 for specfem compatibility
        #coeffs_r = filter_sh(coeffs_r, list(range(21)))

        coeffs_r_filt = truncate_sh(coeffs_r, lmax_filt)
        coeffs_r_filt = filter_sh(coeffs_r_filt, l_set)

        model_key = parse_title(model)

        # Thermo-chemical models have correlated dlnvs and dtopo
        if 'TC' in model_key:
            best_angles[model_key], errors = find_best_angles(
                coeffs_r_filt, coeffs_s20rts_filt)
        # Thermal models have anti-correlated dlnvs and dtopo
        else:
            best_angles[model_key], errors = find_best_angles(
                coeffs_r_filt, -1*coeffs_s20rts_filt)

        coeffs_r_filt_notrunc = filter_sh(coeffs_r, l_set)

        dj_filt = pyshtools.rotate.djpi2(coeffs_r_filt_notrunc.shape[1] - 1)
        dj = pyshtools.rotate.djpi2(coeffs_r.shape[1] - 1)
        coeffs_r_filt_best[model_key] = pyshtools.rotate.SHRotateRealCoef(
            coeffs_r_filt_notrunc,
            best_angles[model_key], dj_filt
        )
        coeffs_r_best[model_key] = pyshtools.rotate.SHRotateRealCoef(
            coeffs_r,
            best_angles[model_key],
            dj
        )

        model_rot_path = os.path.join(
            os.path.expanduser(model_root),
            model.replace('.dat', '_rot.ylm')
        )
        model_rot_path_l20 = os.path.join(
            os.path.expanduser(model_root),
            model.replace('.dat', '_rot_l20.ylm')
        )
        model_rot_path_npy = model_rot_path.replace('.ylm', '')

        spharm_io.write_sh_specfem_block(
            model_rot_path, coeffs_r_best[model_key])
        logging.info(
            f'Wrote rotated SPH model in specfem format to {model_rot_path}')

        spharm_io.write_sh_specfem_block(
            model_rot_path_l20, truncate_sh(coeffs_r_best[model_key], 20))
        logging.info(
            f'Wrote rotated SPH model, truncated to degree 0-20 in '
            f'specfem format to {model_rot_path}')

        spharm_io.write_sh(model_rot_path_npy, coeffs_r_best[model_key])
        logging.info(
            f'Wrote SPH model to numpy object to'
            f' {model_rot_path}')

    #plot best models
    plt.clf()
    fig, axes = display_subplots_cartopy(
        coeffs_r_best, 3, 3, title='l=0-{:d}'.format(lmax),
        cbar_label='CMB topo (km)')
    fig_path = os.path.join(
        os.path.expanduser(figure_root),
        'models_geodynamics_rot_l' + str(lmax) + '.pdf'
    )
    plt.savefig(fig_path)
    logging.info(f'Saved figure to {fig_path}')

    plt.clf()
    fig_filt, axes_filt = display_subplots_cartopy(
        coeffs_r_filt_best, 3, 3, title='l=0-{:d}'.format(lmax_filt),
        cbar_label='CMB topo (km)'
    )
    fig_path_filt = os.path.join(
        os.path.expanduser(figure_root),
        'models_geodynamics_rot_l' + str(lmax_filt) + '.pdf'
    )
    plt.savefig(fig_path_filt)
    logging.info(f'Saved figure to {fig_path}')
