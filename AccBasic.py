#!/usr/bin/env python
# -*- coding: utf-8 -*-

#===============================================================================
#
#  Fundamental functions for the particle accelerator related calculations
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================


from GlobalUnit import *
from math import sqrt, exp
from random import gauss as gauss_shoot
from random import uniform as uniform_shoot
import numpy as np


# ==================================================================
def gamma(w):
    """
    w : kinetic energy [MeV/u]
    """
    return (1. + w / amu_c2)


def gammaT(w, e0):
    """
    w : kinetic energy [MeV]
    e0 : rest energy [MeV]
    """
    return (1. + w / e0)


# ==================================================================
def beta(w):
    """
    w : kinetic energy [MeV/u]
    """
    return sqrt(1. - 1. / (gamma(w) ** 2))


def betaT(w, e0):
    """
    w : kinetic energy [MeV]
    e0 : rest energy [MeV]
    """
    return sqrt(1. - 1. / (gammaT(w, e0) ** 2))


# ==================================================================
def betagamma(w):
    """
    w : kinetic energy [MeV/u]
    """
    return sqrt(gamma(w) ** 2 - 1.)


def betagammaT(w, e0):
    """
    w : kinetic energy [MeV]
    e0 : rest energy [MeV]
    """
    return sqrt(gammaT(w, e0) ** 2 - 1.)


# ==================================================================
def mom(w, A):
    """
    w : kinetic energy [MeV/u]
    A : atomic mass number [integer]
    """
    return betagamma(w) * A * amu * c_light


def momT(w, e0):
    """
    w : kinetic energy [MeV]
    e0 : rest energy [MeV]
    """
    return betagammaT(w, e0) * e0 / c_light


# ==================================================================
def brho(w, zoa):
    """
    w : kinetic energy [MeV/u]
    qm : Q/A charge to mass ratio [arb. unit]
    """
    return betagamma(w) * amu * c_light / float(zoa)


def brhoT(w, e0, q):
    """
    w : kinetic energy [MeV]
    e0 : rest energy [MeV]
    q : charge number [eplus]
    """
    return betagammaT(w, e0) * e0 / c_light / float(q)


# ==================================================================
def dpop(w0, dwow0):
    """
    w0    : kinetic energy [MeV/u]
    dwow0 : energy spread [ratio]
    """
    return betagamma(w0 * (1.0 + dwow0)) / betagamma(w0) - 1.0


def dpopT(w, dwow0, e0):
    """
    w     : kinetic energy [MeV]
    dwow0 : energy spread [ratio]
    e0     : rest energy [MeV]
    """
    return betagammaT(w * (1.0 + dwow0), e0) / betagammaT(w, e0) - 1.0


# ==================================================================
def dwow(w0, dpop0):
    """
    w0    : kinetic energy [MeV/u]
    dpop0 : momentum spread [ratio]
    """
    T0 = gamma(w0) - 1.0
    T1 = sqrt(((dpop0 + 1.0) * betagamma(w0)) ** 2 + 1.0) - 1.0
    return T1 / T0 - 1.0


def dwowT(w, dpop, e0):
    """
    w    : kinetic energy [MeV/u]
    dpop : momentum spread [ratio]
    e0    : rest energy [MeV]
    """
    T0 = gammaT(w, e0) - 1.0
    T1 = sqrt(((dpop + 1.0) * betagammaT(w, e0)) ** 2 + 1.0) - 1.0
    return T1 / T0 - 1.0


# ==================================================================
def pdf(x, mu, sigma):
    """
    Gaussian Standard Distribution Function
    """
    return exp(-0.5 * (x - mu) ** 2 / sigma ** 2) / (sigma * sqrt(twopi))


# ==================================================================
def gen_gauss_kernel(N, sigma):
    """
    This returns a Gaussian Kernel for gaussian smoothing
    N     : size of array. 2N+1 array will be returned
    sigma : standard deviation. if it is close to the middle,
            kernel will be like a delta function,
            if it is close to the end, kernel will be like
            linear average kernel.
    """
    kernel = [pdf(float(i), float(N), sigma) for i in range(2 * N + 1)]
    intg = sum(kernel)
    return [e / intg for e in kernel]


# ==================================================================
def erf(x):
    """
    Error Function
    """
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    # A&S formula 7.1.26
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x)
    return sign * y  # erf(-x) = -erf(x)


# ==================================================================
def cdf(s, x):
    """
    Gaussian cumulative distribution function

    s : standard deviation value
    x : maximum integration limit( from -inf to x )
    """
    return 0.5 * (1.0 + erf(x / sqrt(2.0 * s * s)))


# ==================================================================
def Pxdx(sigma):
    """
    1-Dim integration value of gaussian standard distribution

    sigma : min and max sigma value (Int_(-sigma)^(+sigma) GP(x) dx)

    (example)
    To get 2 sigma probability

      >>> Pxdx(2.0)
      0.9544998742254873

    """
    return 2.0 * (cdf(1., sigma) - 0.5)


def absmax(aa, bb):
    if abs(aa) > abs(bb):
        return aa
    else:
        return bb


def absmaxabs(aa, bb):
    return abs(absmax(aa, bb))


def absmaxlist(one_dim_list):
    abslist = [abs(e) for e in one_dim_list]
    return max(abslist)


def drange(init, final, N):
    delta = (final - init) / float(N - 1)
    return [init + float(i) * delta for i in xrange(N)]


def drange_center(init, final, N):
    delta = (final - init) / float(N)
    return [init + (float(i) + 0.5) * delta for i in xrange(N)]


is_even_number = lambda number: number % 2 == 0


# ==================================================================

class histo2d():
    def __init__(self, minx, maxx, nxbins, miny, maxy, nybins):
        self.minx, self.maxx = minx, maxx
        self.miny, self.maxy = miny, maxy
        self.nxbins, self.nybins = int(nxbins), int(nybins)
        self.dx = (maxx - minx) / float(nxbins)
        self.dy = (maxy - miny) / float(nybins)
        self.data = [[0.0 for j in xrange(self.nybins)] for i in xrange(self.nxbins)]
        self.XE = drange_center(minx, maxx, self.nxbins)  # [minx + self.dx * float(i) for i in range(nxbins)]
        self.YE = drange_center(miny, maxy, self.nybins)  # [miny + self.dy * float(i) for i in range(nybins)]

    def fill(self, x, y, **kwd):
        if x > self.minx and x < self.maxx and y > self.miny and y < self.maxy:
            ix = int((x - self.minx) / self.dx)  # 몫을 구했다.
            iy = int((y - self.miny) / self.dy)

            if kwd.has_key('intensity'):
                self.data[ix][iy] += abs(kwd['intensity'])
            else:
                self.data[ix][iy] += 1.0

    def get_maxcount(self):
        return max([max(yarr) for yarr in self.data])

    def normalize(self, **kwd):
        norm = 1.0
        if kwd.has_key('norm'):
            norm = abs(float(kwd['norm']))
        maxv = self.get_maxcount()
        for i in xrange(self.nxbins):
            for j in xrange(self.nybins):
                self.data[i][j] = float(self.data[i][j]) / float(maxv) * norm

    def serialize(self):
        # XARR = self.XE * self.nybins
        # YARR = list()
        # for yvalue in self.YE:
        #   YARR += [yvalue] * self.nxbins
        try:
            return self.XARR, self.YARR, self.WTARR
        except:
            self.YARR = self.YE * self.nxbins
            self.XARR = list()
            for xvalue in self.XE:
                self.XARR += [xvalue] * self.nybins

            WTARR2D = np.array(self.data)
            self.WTARR = WTARR2D.reshape(self.nxbins * self.nybins)
            return self.XARR, self.YARR, self.WTARR

    def export_as_json(self):
        return {
            'minx': self.minx,
            'maxx': self.maxx,
            'miny': self.miny,
            'maxy': self.maxy,
            'nxbins': self.nxbins,
            'nybins': self.nybins,
            'dx': self.dx,
            'dy': self.dy,
            'dats': self.data
        }

    def import_from_json(self, js):
        self.minx = js['minx']
        self.maxx = js['maxx']
        self.miny = js['miny']
        self.maxy = js['maxy']
        self.nxbins = js['nxbins']
        self.nybins = js['nybins']
        self.dx = js['dx']
        self.dy = js['dy']
        self.data = js['dats']

    # ==========================================================================


# Beam Distributions
# ==========================================================================
def twiss_points(a, b, e, n):
    x_max = sqrt(b * e)
    x_max2 = x_max ** 2
    xp_int = sqrt(e / b)
    x_slop = -a / b
    NP = n
    HNP = NP / 2
    x = [x_max * (-1 + 2.0 * float(i) / float(HNP)) for i in xrange(HNP)]
    xp = []
    for xi in x:
        try:
            xp.append(x_slop * xi + xp_int * sqrt(1. - xi * xi / x_max2))
        except:
            # print xi*xi/x_max2
            xp.append(x_slop * xi)
    # reverse scan
    for i in xrange(HNP, 0, -1):
        xi = x_max * (-1 + 2.0 * float(i) / float(HNP))
        x.append(xi)
        try:
            xp.append(x_slop * xi - xp_int * sqrt(1. - xi * xi / x_max2))
        except:
            # print xi*xi/x_max2
            xp.append(x_slop * xi)
    x.append(x[0])
    xp.append(xp[0])
    return (x, xp)


class TwissPartTr1D:
    def __init__(self):
        self.x = []
        self.xp = []
        self.wt = []
        self.cm_calc_done = False
        self.emit_calc_done = False
        self.cmx = 0.0
        self.cmxp = 0.0

    def add_particle(self, xi, xpi, wi):
        """
        add single particle into this bunch
        """
        self.x.append(xi)
        self.xp.append(xpi)
        self.wt.append(float(wi))
        self.cm_calc_done = False
        self.emit_calc_done = False
        return len(self.x)

    def set_particles(self, xs, xps):
        """
        set particle bunch at once
        """
        self.x = xs
        self.xp = xps
        self.wt = [1.] * len(xs)
        self.cm_calc_done = False
        self.emit_calc_done = False
        return len(self.x)

    def set_particle_bunch(self, xs, xps, wt):
        """
        set particle bunch at once
        """
        self.x = xs
        self.xp = xps
        self.wt = wt
        self.cm_calc_done = False
        self.emit_calc_done = False
        return len(self.x)

    def calc_cm(self):
        """
        calculate center of mass
        """
        self.cmx = 0.0
        self.cmxp = 0.0
        m0 = 0.0
        for i in xrange(len(self.wt)):
            self.cmx += self.x[i] * self.wt[i]
            self.cmxp += self.xp[i] * self.wt[i]
            m0 += self.wt[i]
        self.cmx /= m0
        self.cmxp /= m0
        self.cm_calc_done = True

    def get_cmarray(self):
        """
        returns relative position array of x and xp
        """
        if not self.cm_calc_done:
            self.calc_cm()
        return [e - self.cmx for e in self.x], [e - self.cmxp for e in self.xp]

    def get_cmpoint(self):
        """
        returns center of mass points
        """
        if not self.cm_calc_done:
            self.calc_cm()
        return self.cmx, self.cmxp

    def get_twiss_points(self, N_RMS, N_POINT):
        """
        returns twiss points for plotting
        """
        if not self.emit_calc_done:
            self.get_rms_emittance()

        xarr, xparr = twiss_points(self.alpha_x,
                                   self.beta_x,
                                   N_RMS * self.emit_rms_x,
                                   N_POINT)

        return [e + self.cmx for e in xarr], [e + self.cmxp for e in xparr]

    def get_twiss_points_center(self, N_RMS, N_POINT):
        """
        returns twiss points for plotting
        """
        if not self.emit_calc_done:
            self.get_rms_emittance()

        return twiss_points(self.alpha_x,
                            self.beta_x,
                            N_RMS * self.emit_rms_x,
                            N_POINT)

    def get_rms_emittance(self):
        """
        this returns rms emittance and statistical twiss-parameters

        ( alpha x,      [-]
          beta x,       [m]
          emittance x,  [pi.mm.mrad]
          sigma_x,      [mm]
          sigma_xp )    [mrad]
        """
        if not self.cm_calc_done:
            self.calc_cm()

        N_PARTICLE = len(self.x)

        x_ave = 0.
        x2_ave = 0.
        xp_ave = 0.
        xp2_ave = 0.
        xxp_ave = 0.
        x2xp2_ave = 0.

        for i in xrange(N_PARTICLE):
            wi = self.wt[i]
            xi = self.x[i] - self.cmx
            xi2 = xi ** 2
            xpi = self.xp[i] - self.cmxp
            xpi2 = xpi ** 2
            x_ave += xi * wi
            x2_ave += xi2 * wi
            xp_ave += xpi * wi
            xp2_ave += xpi2 * wi
            xxp_ave += xi * xpi * wi
            x2xp2_ave += xi2 * xpi2 * wi

        np = 0.0
        for e in self.wt:
            np += float(e)
        x_ave /= np
        xp_ave /= np
        x2_ave /= np
        xp2_ave /= np
        xxp_ave /= np
        x2xp2_ave /= np

        x_rms = sqrt(x2_ave - x_ave ** 2)
        xp_rms = sqrt(xp2_ave - xp_ave ** 2)
        xxp_rms = xxp_ave - x_ave * xp_ave

        emit_rms_x = sqrt(x2_ave * xp2_ave - xxp_ave ** 2)
        beta_x = (x_rms * x_rms) / emit_rms_x
        gamma_x = (xp_rms * xp_rms) / emit_rms_x
        alpha_x = - xxp_rms / emit_rms_x

        self.alpha_x = alpha_x
        self.beta_x = beta_x
        self.gamma_x = gamma_x
        self.emit_rms_x = emit_rms_x
        self.x_rms = x_rms
        self.xp_rms = xp_rms

        self.emit_calc_done = True

        return (alpha_x, beta_x, emit_rms_x, x_rms, xp_rms)


ParticleSet1D = TwissPartTr1D


class TrBeam:
    def __init__(self):
        self.ax, self.ay = 0., 0.
        self.bx, self.by = 0., 0.
        self.ex, self.ey = 0., 0.
        self.w0, self.dwow = 0., 0.

    def set_abe(self, ax, bx, ex, ay, by, ey):
        self.ax, self.ay = ax, ay
        self.bx, self.by = bx, by
        self.ex, self.ey = ex, ey

    def set_sizediv(self, x, xp, y, yp):
        self.ax, self.ay = 0., 0.
        self.ex, self.ey = x * xp, y * yp
        self.bx, self.by = x / xp, y / yp

    def set_ke(self, w0, dwow):
        self.w0, self.dwow = w0, dwow

    def gen_gaussdist(self, np):
        if np == 1:
            return ([0.], [0.], [0.], [0.], [self.w0])
        x, xp, y, yp, w = [], [], [], [], []
        xslop, yslop = -self.ax / self.bx, -self.ay / self.by
        xpint, ypint = sqrt(self.ex / self.bx), sqrt(self.ey / self.by)
        xmax, ymax = sqrt(self.ex * self.bx), sqrt(self.ey * self.by)
        for i in xrange(np):
            xi = gauss_shoot(0.0, xmax)
            xpi = gauss_shoot(xslop * xi, xpint)
            yi = gauss_shoot(0.0, ymax)
            ypi = gauss_shoot(yslop * yi, ypint)
            wi = gauss_shoot(self.w0, self.w0 * self.dwow)
            x.append(xi)
            xp.append(xpi)
            y.append(yi)
            yp.append(ypi)
            w.append(wi)
        return (x, xp, y, yp, w)

    def gen_uniformdist(self, np):
        if np == 1:
            return ([0.], [0.], [0.], [0.], [self.w0])
        # 2Dim water-bag distribution
        gx = (1.0 + self.ax ** 2.0) / self.bx
        gy = (1.0 + self.ay ** 2.0) / self.by
        ex, ey = 4.0 * self.ex, 4.0 * self.ey
        x, xp, y, yp, w = [], [], [], [], []
        xslop, yslop = -self.ax / self.bx, -self.ay / self.by
        xpint, ypint = sqrt(ex / self.bx), sqrt(ey / self.by)
        xpmax, ypmax = sqrt(ex * gx), sqrt(ey * gy)
        xmax, ymax = sqrt(ex * self.bx), sqrt(ey * self.by)

        x_ave, x2_ave, xp_ave, xp2_ave = 0., 0., 0., 0.
        y_ave, y2_ave, yp_ave, yp2_ave = 0., 0., 0., 0.
        xxp_ave, xxp2_ave, yyp_ave, yyp2_ave = 0., 0., 0., 0.
        i = 0
        while i != np:
            xi = uniform_shoot(-xmax, xmax)
            xpi = uniform_shoot(-xpmax, xpmax)
            yi = uniform_shoot(-ymax, ymax)
            ypi = uniform_shoot(-ypmax, ypmax)
            # ENERGIES ARE NOT GAUSS DIST
            wi = gauss_shoot(self.w0, self.w0 * self.dwow)
            if gx * xi * xi + 2.0 * self.ax * xpi * xi + self.bx * xpi * xpi < ex and \
                    gy * yi * yi + 2.0 * self.ay * ypi * yi + self.by * ypi * ypi < ey:
                x.append(xi)
                xp.append(xpi)
                y.append(yi)
                yp.append(ypi)
                w.append(wi)
                i += 1
        return (x, xp, y, yp, w)


class fitvalue_ab:
    def __init__(self, ax, bx, ay, by):
        self.ax = ax
        self.bx = bx
        self.gx = (ax ** 2.0 + 1.0) / bx
        self.ay = ay
        self.by = by
        self.gy = (ay ** 2.0 + 1.0) / by
        self.dx_status = False
        self.dx = 0.0
        self.dpx = 0.0

    def initialize(self):
        self.gx = (self.ax ** 2.0 + 1.0) / self.bx
        self.gy = (self.ay ** 2.0 + 1.0) / self.by
