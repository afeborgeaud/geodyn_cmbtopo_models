from scipy.special import lpmv
import numpy as np

def read_sph(filename):
    file = open(filename,'r')
    coeffs = list()
    for line in file:
        coeffs.append(np.asarray(line.split(), dtype=float))
    file.close()
    return coeffs

def fact(l, m):
    x = 1.
    for i in range(l-m+1, l+m+1):
        x *= i
    return 1./x

def ylm(l, m, phi, theta):
    y = (lpmv(abs(m), l, np.cos(theta)) * (-1)**abs(m)
         * np.sqrt((2 * l + 1) * fact(l, abs(m))))
    if (m < 0):
        y *= np.sqrt(2) * np.sin(abs(m)*phi) 
    elif (m > 0):
        y *= np.sqrt(2) * np.cos(abs(m)*phi) 
    return y

def eval(phi, theta, coeffs):
    y = 0
    for l,cs in enumerate(coeffs):
        for im,c in enumerate(cs):
            m = im - int(len(cs) / 2)
            y += c * ylm(l, m, phi, theta)
    return y
