# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 17:36:02 2013

@author: BRENNICH
"""

import numpy 
import scipy.integrate as scint

print("imported rti")

def calcVc(dat, Rg, dRg, I0, dI0, imin):
    """Calculates the Rambo-Tainer invariant Vc, including extrapolation to q=0

    Arguments: 
    @param dat:  data in q,I,dI format, cropped to maximal q that should be used for calculation (normally 2 nm-1)
    @param Rg,dRg,I0,dI0:  results from Guinier approximation/autorg
    @param imin:  minimal index of the Guinier range, below that index data will be extrapolated by the Guinier approximation
    @returns: Vc and an error estimate based on non-correlated error propagation
    """
    dq = dat[1, 0] - dat[0, 0]
    qmin = dat[imin, 0]
    qlow = numpy.arange(0, qmin, dq)

    lowqint = scint.trapz((qlow * I0 * numpy.exp(-(qlow * qlow * Rg * Rg) / 3.0)), qlow)
    dlowqint = scint.trapz(qlow * numpy.sqrt((numpy.exp(-(qlow * qlow * Rg * Rg) / 3.0) * dI0) ** 2 + ((I0 * 2.0 * (qlow * qlow) * Rg / 3.0) * numpy.exp(-(qlow * qlow * Rg * Rg) / 3.0) * dRg) ** 2), qlow)
    vabs = scint.trapz(dat[imin:, 0] * dat[imin:, 1], dat[imin:, 0])
    dvabs = scint.trapz(dat[imin:, 0] * dat[imin:, 2], dat[imin:, 0])
    vc = I0 / (lowqint + vabs)
    dvc = (dI0 / I0 + (dlowqint + dvabs) / (lowqint + vabs)) * vc
    return (vc, dvc)

def RamboTainerInvariant(dat, Rg, dRg, I0, dI0, imin, qmax=2):
    """calculates the invariants Vc and Qr from the Rambo&Tainer 2013 Paper,
    also the the mass estimate based on Qr for proteins

    Arguments: 
    @param dat: data in q,I,dI format, q in nm-1
    @parma Rg,dRg,I0,dI0: results from Guinier approximation
    @parma imin: minimal index of the Guinier range, below that index data will be extrapolated by the Guinier approximation
    @parma qmax: maximum q-value for the calculation in nm-1
    @return: dict with Vc, Qr and mass plus errors
    """
    scale_prot = 1.0 / 0.1231
    power_prot = 1.0

    imax = abs(dat[:, 0] - qmax).argmin()
    vc = calcVc(dat[:imax, :], Rg, dRg, I0, dI0, imin)

    qr = vc[0] ** 2 / (Rg)
    mass = scale_prot * qr ** power_prot

    dqr = qr * (dRg / Rg + 2 * ((vc[1]) / (vc[0])))
    dmass = mass * dqr / qr

    return {'Vc': vc[0], 'dVc': vc[1], 'Qr': qr, 'dQr': dqr, 'mass': mass, 'dmass': dmass}
