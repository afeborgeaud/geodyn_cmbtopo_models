"""
2021.08.14
Anselme F.E. Borgeaud

Utility functions to read and write CMB topography
models in spherical harmonic.
"""

import pyshtools
import numpy as np


def read_sh_c_to_r(file, nondim=True, lmaxin=-1):
    """Read complex spherical harmonics coefficients with format:
    row 0: lmin lmax
    rows 1, (lmax+1)**2: l m c_real c_imag
    m = -l, l
    nondim=True is for nondimentionalized coefficients using the mantle
    width=2891 km

    return: lmax, c[k, l, m]
    c.shape = (2, lmax+1, lmax+1)
    c[0, l, m] is the cosine part
    c[1, l, m] is the sine part
    normalizatino input: orthonormalized
    normalization output: 4-pi (geodesy)
    """

    to_float = lambda x: float(x.replace('D', 'e'))
    converters = {2: to_float, 3: to_float}
    with open(file) as f:
        lmin, lmax = [int(x) for x in f.readline().split()]
    arr = np.loadtxt(file, dtype=np.float64,
                     skiprows=1, converters=converters, encoding='utf-8')
    ls = arr[:, 0].astype(np.int32)
    ms = arr[:, 1].astype(np.int32)
    coeffs_01 = arr[:, 2]
    coeffs_02 = arr[:, 3]

    ls_pos = ls[ms >= 0]
    ms_pos = ms[ms >= 0]
    coeffs_1 = coeffs_01[ms >= 0]
    coeffs_2 = coeffs_02[ms >= 0]

    coeffs_fmt = np.zeros((2, lmax + 1, lmax + 1), dtype=np.float64)

    for l in ls_pos:
        coeffs_fmt[0, l, ms_pos[ls_pos == l]] = coeffs_1[ls_pos == l]
        coeffs_fmt[1, l, ms_pos[ls_pos == l]] = coeffs_2[ls_pos == l]

    coeffs_fmt[:, :] *= np.power(-1, list(range(lmax + 1)))

    if (nondim):
        coeffs_fmt *= 2891.

    # convention=2: orthonormalized to real geodesy 4-pi
    coeffs_r = pyshtools.shio.SHctor(coeffs_fmt, convention=2)

    if (lmaxin != -1):
        coeffs_r = coeffs_r[:, :lmaxin + 1, :lmaxin + 1]
    else:
        lmaxin = lmax

    return lmaxin, coeffs_r


def read_sh_r_to_r(file, nondim=True, lmaxin=-1):
    """Read real spherical harmonics coefficients with format:
    row 0: lmin lmax
    rows 1, (lmax+1)**2: l m c_real 0.
    m = -l, l
    nondim=True is for nondimentionalized coefficients using the mantle
    width=2891 km

    return: lmax, c[k, l, m]
    c.shape = (2, lmax+1, lmax+1)
    c[0, l, m] is the cosine part
    c[1, l, m] is the sine part
    normalizatino input: orthonormalized
    normalization output: 4-pi (geodesy)
    """

    to_float = lambda x: float(x.replace('D', 'e'))
    converters = {2: to_float, 3: to_float}
    with open(file) as f:
        lmin, lmax = [int(x) for x in f.readline().split()]
    arr = np.loadtxt(file, dtype=np.float64,
                     skiprows=1, converters=converters, encoding='utf-8')
    ls = arr[:, 0].astype(np.int32)
    ms = arr[:, 1].astype(np.int32)
    coeffs = arr[:, 2]

    ls_1 = ls[ms >= 0]
    ms_1 = ms[ms >= 0]
    coeffs_1 = coeffs[ms >= 0]

    ls_2 = ls[ms < 0]
    ms_2 = np.abs(ms[ms < 0])
    coeffs_2 = coeffs[ms < 0]

    coeffs_fmt = np.zeros((2, lmax + 1, lmax + 1), dtype=np.complex128)

    for l in ls_1:
        coeffs_fmt[0, l, ms_1[ls_1 == l]] = coeffs_1[ls_1 == l]
        coeffs_fmt[1, l, ms_2[ls_2 == l]] = coeffs_2[ls_2 == l]

    coeffs_fmt[:, :] *= np.power(-1, list(range(lmax + 1)))

    if (nondim):
        coeffs_fmt *= 2891.

    coeffs_r = pyshtools.shio.convert(coeffs_fmt,
                                      normalization_in='ortho',
                                      normalization_out='4pi')

    if (lmaxin != -1):
        coeffs_r = coeffs_r[:, :lmaxin + 1, :lmaxin + 1]
    else:
        lmaxin = lmax

    return lmaxin, coeffs_r


def read_sh_r_block_specfem(file):
    """ Read real spherical harmonic coefficients
    in block format used in specfem
    """

    with open(file, 'r') as f:
        lmax = int(f.readline())
        tmps = list()
        for line in f:
            for c in np.asarray(line.split(), dtype=float):
                tmps.append(c)

    coeffs = np.zeros((2, lmax + 1, lmax + 1))
    count = 0
    for l in range(lmax + 1):
        for im in range(2 * l + 1):
            m = abs(im2m(im))
            k = 1 if (im % 2 == 0 and im != 0) else 0
            coeffs[k, l, m] = tmps[count]
            count += 1

    # re-normalize from specfem convention (in particular, specfem takes
    # dDepth = - dr)
    coeffs[..., 1:] /= np.sqrt(2)
    coeffs[:, :] *= np.power(-1, list(range(lmax + 1)))
    coeffs /= -1 * np.sqrt(4 * np.pi)

    return lmax, coeffs


def im2m(im):
    if im == 0:
        return 0
    elif im % 2 == 1:
        return int((im + 1) / 2)
    elif im % 2 == 0:
        return -int(im / 2)
    else:
        return int(1e10)


def read_sh_r_specfem(file):
    with open(file, 'r') as f:
        lmax = int(f.readline())

    coeffs_r = np.zeros((2, lmax + 1, lmax + 1))
    with open(file, 'r') as f:
        file.readline()
        for l, line in enumerate(f):
            for im, x in enumerate(line.split()):
                m = im2m(im)
                k = 0 if im % 2 == 0 else 1
                coeffs_r[k, l, m] = x

    return coeffs_r


def read_sh(file_npy):
    coeffs = np.load(file_npy)
    lmax = coeffs.shape[1] - 1
    return lmax, coeffs


def write_sh(file_npy, coeffs):
    np.save(file_npy, coeffs)


def write_sh_specfem_block(file, coeff):
    def write_coeff(file, c, pos, lmax):
        if (pos % (lmax + 1) == 0) and (pos != 0):
            file.write('{:.5e}\n'.format(c))
        else:
            file.write('{:.5e} '.format(c))
        return pos + 1

    lmax = coeff.shape[1] - 1

    # normalize for specfem convention (in particular, specfem takes dDepth
    # = - dr)
    coeff_norm = coeff.copy()
    coeff_norm[..., 1:] *= np.sqrt(2)
    coeff_norm[:, :] *= np.power(-1, list(range(lmax + 1)))
    coeff_norm *= -1 * np.sqrt(4 * np.pi)

    pos = 1
    with open(file, 'w') as f:
        f.write('{}\n'.format(lmax))
        for l in range(coeff.shape[1]):
            pos = write_coeff(f, coeff_norm[0, l, 0], pos, lmax)
            for m in range(1, l + 1):
                pos = write_coeff(f, coeff_norm[0, l, m], pos, lmax)
                pos = write_coeff(f, coeff_norm[1, l, m], pos, lmax)
