#!/usr/bin/env python
# -*- coding: utf-8 -*-
#===============================================================================
#
#  Linear beam optics
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================


from GlobalUnit import *
from AccBasic import *
from math import ceil, cosh, sinh, tan, sin, cos
from numpy import dot, transpose as tr
from Vector import ThreeVector, RotationMatrix
from random import gauss as gauss_shoot
from wx import ProgressDialog as Progress, PD_ELAPSED_TIME, PD_AUTO_HIDE


def gen_identity_matrix():
  return [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
          [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]


#########################################################################
class BeamSpec:
  def __init__(self, rmi, qi, wi):
    """
    (1) rmi : rest energy (MeV) (rest mass = MeV/c2)
    (2) qi  : charge state (integer number)
    (3) wi  : kinetic energy (MeV) (total energy)
    """
    self.rest_mass = rmi
    self.charge = float(qi)
    self.kinetic_energy = wi
    self.set_particle_vars()
    
    self.ax, self.ay, self.az = 0., 0., 0. # alpha
    self.bx, self.by, self.bz = 0., 0., 0. # beta
    self.gx, self.gy, self.gz = 0., 0., 0. # gamma
    self.ex, self.ey, self.ez = 0., 0., 0. # emittance
    
    self.x, self.xp, self.xs = 0., 0., 0.
    self.y, self.yp, self.ys = 0., 0., 0.
    self.z, self.zp, self.zs = 0., 0., 0.
    
    self.dx, self.dpx = 0., 0.
    self.dy, self.dpy = 0., 0.
    self.dpop = 0.0

  def set_particle_vars(self):
    self.atomic_mass            = self.rest_mass / amu_c2 
    self.kinetic_energy_per_amu = self.kinetic_energy / self.atomic_mass
    
    self.qm     = self.charge / self.atomic_mass
    self.gam    = gamma(self.kinetic_energy_per_amu)
    self.gam2   = self.gam**2.0
    self.bet    = beta(self.kinetic_energy_per_amu)
    self.bet2   = self.bet**2.0
    self.bg     = self.bet * self.gam
    self.bg2    = self.bg**2
    self.br     = brho(self.kinetic_energy_per_amu, self.qm)
    
  def reset(self, rmi, qi, wi):
    self.rest_mass = rmi
    self.charge = float(qi)
    self.kinetic_energy = wi
    self.set_particle_vars()
    
    
  def reset_by_qme(self, qm, e):
    self.qm = qm
    self.kinetic_energy_per_amu  = e
    self.gam    = gamma(self.kinetic_energy_per_amu)
    self.gam2   = self.gam**2.0
    self.bet    = beta(self.kinetic_energy_per_amu)
    self.bet2   = self.bet**2.0
    self.bg     = self.bet * self.gam
    self.bg2    = self.bg**2
    self.br     = brho(self.kinetic_energy_per_amu, qm)
    
    if self.rest_mass != amu:
      self.charge = self.qm * self.atomic_mass
    else:
      self.rest_mass
      self.charge = qm
      self.kinetic_energy = e

  def set_xabe(self, ax, bx, ex):
    """
    ax : alpha
    bx : beta
    ex : emittance
    """
    self.ax = ax
    self.bx = bx
    self.gx = (ax**2.0 + 1.0) / bx
    self.ex = ex

    self.x = sqrt(bx * ex)
    self.xp = sqrt(self.gx * ex)
    self.xs = -ax / bx

  def set_yabe(self, ay, by, ey):
    """
    ay : alpha
    by : beta
    ey : emittance
    """
    self.ay = ay
    self.by = by
    self.gy = (ay**2.0 + 1.0) / by
    self.ey = ey

    self.y = sqrt(by * ey)
    self.yp = sqrt(self.gy * ey)
    self.ys = -ay / by

  def set_xxp(self, x, xp):
    """
    x : x [beam size]
    xp : xp [divergence angle]
    """
    self.x = x
    self.xp = xp
    self.ex = x * xp

    self.bx = x**2.0 / self.ex
    self.gx = xp**2.0 / self.ex

  def set_yyp(self, y, yp):
    """
    y : y [beam size]
    yp : yp [divergence angle]
    """
    self.y = y
    self.yp = yp
    self.ey = y * yp

    self.by = y**2.0 / self.ey
    self.gy = yp**2.0 / self.ey
    
  def set_dpop(self, dpop):
    self.dpop = dpop
  
  def set_zabe(self, a, b, e):
    """
    a : alpha, dimensionless quantity
    b : beta, MeV/deg
    e : emittance, deg*MeV
    """
    self.az = a
    self.bz = b
    self.gz = (a**2.0 + 1.0) / b
    self.ez = e

    self.z = sqrt(b * e)
    self.zp = sqrt(self.gz * e)
    self.zs = -a / b

  def set_dispersion(self, dx, dpx, dy, dpy):
    """
    dx  : Dx
    dpx : D'x
    dy  : Dy
    dpy : D'y
    """
    self.dx = dx
    self.dpx = dpx
    self.dy = dy
    self.dpy = dpy

  def get_sigma_matrix(self):
    return [[ self.bx, -self.ax, 0., 0., 0., 0.],
             [-self.ax,  self.gx, 0., 0., 0., 0.],
             [0., 0.,  self.by, -self.ay, 0., 0.],
             [0., 0., -self.ay,  self.gy, 0., 0.],
             [0., 0., 0., 0., 0., 0.],
             [0., 0., 0., 0., 0., 0.]]

  def get_data(self):
    s = {'rest_mass': self.rest_mass,
         'charge': self.charge,
         'kinetic_energy': self.kinetic_energy,
         'ax': self.ax, 'ay': self.ay,
         'bx': self.bx, 'by': self.by,
         'ex': self.ex, 'ey': self.ey,
         'dpop': self.dpop }
    return s
  
  def gen_gaussian(self, np):
    x, xp, y, yp, zp = [], [], [], [], []
    xslop, yslop = -self.ax/self.bx, -self.ay/self.by
    xpint, ypint = sqrt(self.ex/self.bx), sqrt(self.ey/self.by)
    for i in xrange(np):
      xi = gauss_shoot(0.0, self.x)
      xpi = gauss_shoot(xslop*xi, xpint)
      yi = gauss_shoot(0.0, self.y)
      ypi = gauss_shoot(yslop*yi, ypint)
      zpi = gauss_shoot(0.0, self.dpop)
      x.append(xi)
      xp.append(xpi)
      y.append(yi)
      yp.append(ypi)
      zp.append(zpi)
    return (x, xp, y, yp, zp)
  
  def gen_boundary(self, np):
    """
      x, xp, y, yp are boundary!!
      z, zp is just zero!!
    """
    zeros = [0.0 for x in xrange(np/2)]
    x, xp = twiss_points(self.ax, self.bx, self.ex, np/2)
    x.pop()
    xp.pop()
    x.extend(zeros)
    xp.extend(zeros)
    y, yp = [], []
    y.extend(zeros)
    yp.extend(zeros)
    ye, ype = twiss_points(self.ay, self.by, self.ey, np/2)
    ye.pop()
    ype.pop()
    y.extend(ye)
    yp.extend(ype)
    zp = []
    zp.extend(zeros)
    zp.extend(zeros)
    return (x, xp, y, yp, zp)
  
  def is_inside(self, x, xp, y, yp):
    xres = self.ex > self.gx * x*x + 2.0 * self.ax * x*xp + self.bx * xp*xp
    yres = self.ey > self.gy * y*y + 2.0 * self.ay * y*yp + self.by * yp*yp
    return xres and yres
  
    

#########################################################################
class OptBase(object):
  def __init__(self, t, **prm):
    self.type            = t
    self.arg             = prm
    self.visible         = True
    self.beam            = 0

  #def R(self):
  #  return self.R(1)

  def __repr__(self):
    return self.type + " = " + str(self.arg)

  
# DRIFT
#########################################################################
class Drift(OptBase):
  def __init__(self, **kwd):
    """
    l : length
    r : beam pipe inner radius
    """
    super(Drift, self).__init__("drift", **kwd)
    self.visible = False
    
  
  def R(self, l):
    return [[1.0,   l, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0,   l, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, l/self.beam.gam2],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., self.arg['l'])

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., 0.5 * self.arg['l'])
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False
  
  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']
    
# QUADRUPOLE
#########################################################################
class Quadrupole(OptBase):
  """
  """
  def __init__(self, **kwd):
    """
    G  : focusing gradient  
    l  : effective length
    r  : beam pipe inner radius
    """
    super(Quadrupole, self).__init__("quadrupole", **kwd)
    
  def R(self, l):
    # if there is G: focusing gradient
    G = self.arg['G']
    
    if G > 0.:
      k = sqrt(abs(G / self.beam.br))
      xc = cos(k*l)
      xs = sin(k*l)
      yc = cosh(k*l)
      ys = sinh(k*l)
      kxs = -k*xs
      kys = k*ys
    elif G < 0.:
      k = sqrt(abs(-G / self.beam.br))
      xc = cosh(k*l)
      xs = sinh(k*l)
      yc = cos(k*l)
      ys = sin(k*l)
      kxs = k*xs
      kys = -k*ys
    elif G == 0.:
      return [[1.0,   l, 0.0, 0.0, 0.0, 0.0],
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 1.0,   l, 0.0, 0.0],
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 1.0, l/self.beam.gam2],
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

    return [[ xc, xs/k, 0.0,  0.0, 0.0, 0.0],
            [kxs,   xc, 0.0,  0.0, 0.0, 0.0],
            [0.0,  0.0,  yc, ys/k, 0.0, 0.0],
            [0.0,  0.0, kys,   yc, 0.0, 0.0],
            [0.0,  0.0, 0.0,  0.0, 1.0, l/self.beam.gam2],
            [0.0,  0.0, 0.0,  0.0, 0.0, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., self.arg['l'])

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., 0.5 * self.arg['l'])
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False

  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']

# SOLENOID
#########################################################################
class Solenoid(OptBase):
  """
  """
  def __init__(self, **kwd):
    """
    B  : magnetic field strength  
    l  : effective length
    r  : beam pipe inner radius
    """
    super(Solenoid, self).__init__("solenoid", **kwd)
    
  def R(self, l):
    B = self.arg['B']
    if B == 0.:
      return [[1.0,   l, 0.0, 0.0, 0.0, 0.0],
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 1.0,   l, 0.0, 0.0],
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 1.0, l/self.beam.gam2],
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
    else:
      k = 0.5 * B / self.beam.br
      c = cos(k*l)
      s = sin(k*l)
      cs = c*s
      c2 = c**2
      s2 = s**2
      return [[   c2,  cs/k,    cs, s2/k, 0.0, 0.0],
              [-cs*k,    c2, -s2*k,   cs, 0.0, 0.0],
              [  -cs, -s2/k,    c2, cs/k, 0.0, 0.0],
              [ s2*k,   -cs, -cs*k,   c2, 0.0, 0.0],
              [  0.0,   0.0,   0.0,  0.0, 1.0, l/self.beam.gam2],
              [  0.0,   0.0,   0.0,  0.0, 0.0, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., self.arg['l'])

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., 0.5 * self.arg['l'])
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0:
      return True
    else: 
      return False

  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']



# ThinLens
#########################################################################
class ThinLens(OptBase):
  """
  """
  def __init__(self, **kwd):
    """
    Fx : x focal length  
    Fy : y focal length  
    Fz : z focal length  
    r  : beam pipe inner radius
    """
    super(ThinLens, self).__init__("thinlens", **kwd)
    self.visible = False
    self.arg['l'] = 0.0
    
  def R(self, l):
    # if there is G: focusing gradient
    Fx = self.arg['Fx']
    Fy = self.arg['Fy']
    Fz = self.arg['Fz']
    if Fx < 1.0 * mm:
      Fx = 1.0 * mm
    if Fy < 1.0 * mm:
      Fy = 1.0 * mm
    if Fz < 1.0 * mm:
      Fz = 1.0 * mm

    return [[1.0, 0.0,     0.0, 0.0,     0.0, 0.0],
            [-1.0/Fx, 1.0, 0.0, 0.0,     0.0, 0.0],
            [0.0, 0.0,     1.0, 0.0,     0.0, 0.0],
            [0.0, 0.0,     -1.0/Fy, 1.0, 0.0, 0.0],
            [0.0, 0.0,     0.0, 0.0,     1.0, 0.0],
            [0.0, 0.0,     0.0, 0.0,     -self.beam.gam2/Fz, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False

  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']
    


class QuadrupoleAsym(OptBase):
  """
  """
  def __init__(self, **kwd):
    """
    B0 : magnetic field strength
    r  : aperture size
      or
    G  : focusing gradient  
    l  : effective length
      and common component
    a  : a ratio of Gx/Gy(=G/(1/a * G))
    """
    super(QuadrupoleAsym, self).__init__("quadrupoleasym", **kwd)
    
  def R(self, l):
    # if there is G: focusing gradient
    G = self.arg['G']
    a = self.arg['a']
    Gx = G
    Gy = G/a
    
    if G > 0.:
      kx = sqrt(abs(Gx / self.beam.br))
      ky = sqrt(abs(Gy / self.beam.br))
      xc = cos(kx*l)
      xs = sin(kx*l)
      yc = cosh(ky*l)
      ys = sinh(ky*l)
      kxs = -kx*xs
      kys = ky*ys
    elif G < 0.:
      kx = sqrt(abs(-Gx / self.beam.br))
      ky = sqrt(abs(-Gy / self.beam.br))
      xc = cosh(kx*l)
      xs = sinh(kx*l)
      yc = cos(ky*l)
      ys = sin(ky*l)
      kxs = kx*xs
      kys = -ky*ys
    elif G == 0.:
      return [[1.0,   l, 0.0, 0.0, 0.0, 0.0],
              [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
              [0.0, 0.0, 1.0,   l, 0.0, 0.0],
              [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 1.0, l/self.beam.gam2],
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

    return [[ xc, xs/kx,  0.0,   0.0, 0.0, 0.0],
            [kxs,    xc,  0.0,   0.0, 0.0, 0.0],
            [0.0,   0.0,   yc, ys/kx, 0.0, 0.0],
            [0.0,   0.0,   kys,   yc, 0.0, 0.0],
            [0.0,   0.0,   0.0,  0.0, 1.0, l/self.beam.gam2],
            [0.0,   0.0,   0.0,  0.0, 0.0, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., self.arg['l'])

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., 0.5 * self.arg['l'])
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False

  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']
    


# DIPOLE
#########################################################################
class Dipole(OptBase):
  """
  Bending(without edge angle)
  """
  def __init__(self, **kwd):
    """
    a : angle
    r : curvature of radius
    n : field index (r/By * dBy/dr) r~=x
    ef: edge angle front (beam comes from)
    eb: edge angle back (beam goes to)
    g : full gap width
    pw : full pole width
    """
    super(Dipole, self).__init__("dipole", **kwd)
    self.arg['l'] =  abs(kwd['a']) * kwd['r']
    

  def R(self, l):
    a = self.arg['a']
    r = abs(self.arg['r'])
    n = self.arg['n']

    h  = a / (r * abs(a))
    if n == 0.0:
      kx = abs(h)
      kxl = kx * l
      cx = cos(kxl)
      sx = sin(kxl)
      ky = 0.0     #sqrt(n * h**2.0)
      kyl = 0.0    #ky * l  
      cy = 1.0     #cos(kyl)
      sy = 0.0     #sin(kyl) 
      zz = (sx - kxl * self.beam.bet2) / (kx**3.0 * r**2.0) + l * (1.0 - 1.0 / (r * kx)**2.0) / self.beam.gam2
      xz12 = h * (1.0 - cx) / kx**2.0
      xz22 = h * sx / kx
      zx11 = -xz22
      zx12 = -xz12
      return [[  cx,   sx/kx, 0.0,   0.0, 0.0, xz12],
              [-kx*sx,    cx, 0.0,   0.0, 0.0, xz22],
              [   0.0,   0.0, 1.0,   l,   0.0,  0.0],
              [   0.0,   0.0, 0.0,   1.0, 0.0,  0.0],
              [  zx11,  zx12, 0.0,   0.0, 1.0,   zz],
              [   0.0,   0.0, 0.0,   0.0, 0.0,  1.0]]
    
    # not yet corrected
    elif n == 1.0:
      kx = 0.0    #sqrt((1.0 - n) * h**2.0)
      kxl = 0.0   #kx * l
      cx = 1.0    #cos(kxl)
      sx = 0.0    #sin(kxl)
      ky = abs(h)   #sqrt(n * h**2.0)
      kyl = ky * l  
      cy = cos(kyl) 
      sy = sin(kyl) 
      zz = l / self.beam.gam2 #(sx - kxl * beta(e)**2.0) / (kx**3.0 * r**2.0) + l * (1.0 - 1.0 / (r * kx)**2.0) / gamma(e)**2.0
      xz12 = 0.0 #h * (1.0 - cx) / kx**2.0
      xz22 = h * l # h * sx / kx
      zx11 = -xz22
      zx12 = -xz12
      return [[   1.0,     l,    0.0,   0.0, 0.0, xz12],
              [   0.0,   1.0,    0.0,   0.0, 0.0, xz22],
              [   0.0,   0.0,     cy, sy/ky, 0.0,  0.0],
              [   0.0,   0.0, -ky*sy,    cy, 0.0,  0.0],
              [  zx11,  zx12,    0.0,   0.0, 1.0,   zz],
              [   0.0,   0.0,    0.0,   0.0, 0.0,  1.0]]

    else:
      kx = sqrt((1.0 - n) * h**2.0)
      kxl = kx * l
      cx = cos(kxl)
      sx = sin(kxl)
      ky = sqrt(n * h**2.0)
      kyl = ky * l
      cy = cos(kyl)
      sy = sin(kyl)
      zz = (sx - kxl * self.beam.bet2) / (kx**3.0 * r**2.0) + l * (1.0 - 1.0 / (r * kx)**2.0) / self.beam.gam2
      xz12 = h * (1.0 - cx) / kx**2.0
      xz22 = h * sx / kx
      zx11 = -xz22
      zx12 = -xz12
      return [[  cx,   sx/kx,    0.0,   0.0, 0.0, xz12],
              [-kx*sx,    cx,    0.0,   0.0, 0.0, xz22],
              [   0.0,   0.0,     cy, sy/ky, 0.0,  0.0],
              [   0.0,   0.0, -ky*sy,    cy, 0.0,  0.0],
              [  zx11,  zx12,    0.0,   0.0, 1.0,   zz],
              [   0.0,   0.0,    0.0,   0.0, 0.0,  1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    yhat = rot_m * ThreeVector(0.,1.,0.)
    xhat = rot_m * ThreeVector(1.,0.,0.)
    new_rot = RotationMatrix()
    new_rot.rotate(-self.arg['a'], yhat) # new R
    return pre_vec + self.arg['r'] * (new_rot * xhat - xhat)
    
    #r = RotationMatrix(R = rot_m)
    #r.rotate(0.5 * self.arg['a'], rot_m * ThreeVector(0.,1.,0.))
    #d = self.arg['r'] * sqrt(2.0 * (1. - cos(self.arg['a'])))
    #return pre_vec + r * ThreeVector(0.,0.,d)

  def GetCenterPosition(self, pre_vec, rot_m):
    #pre_rot = RotationMatrix(R = rot_m)
    yhat = rot_m * ThreeVector(0.,1.,0.)
    xhat = rot_m * ThreeVector(1.,0.,0.)
    new_rot = RotationMatrix()
    new_rot.rotate(-0.5 * self.arg['a'], yhat) # new R
    return pre_vec + self.arg['r'] * (new_rot * xhat - xhat)
    
  def GetNextRotMatrix(self, rot_m):
    r = RotationMatrix(R = rot_m)
    r.rotate(-self.arg['a'], rot_m * ThreeVector(0,1.,0))
    return r

  def IsAlive(self, x, y):
    hx = 0.5 * self.arg['pw']
    hy = 0.5 * self.arg['g']
    if (-hx < x) and (x < hx) and (-hy < y) and (y < hy): 
      return True
    else: 
      return False

  def GetXBoundary(self):
    try:
      return 0.5 * self.arg['pw']
    except:
      return 12 * cm
  
  def GetYBoundary(self):
    try:
      return 0.5 * self.arg['g']
    except:
      return 4 * cm



# DIPOLE EDGE
#########################################################################
class DipoleEdge(OptBase):
  """
  Edge angle calculation
  """
  def __init__(self, **kwd):
    """
    b : edge angle
    r : curvature of radius
    g : full gap width
    pw : pole width
    """
    super(DipoleEdge, self).__init__("dipoleEdge", **kwd)
    self.visible = False
    self.arg['l'] = 0.0

  def R(self, l):
    r   = 1./self.arg['r']
    gor = self.arg['g'] * r
    b   = self.arg['b']
    #K1  = 0.45
    #K2  = 2.8
    try:
      K1  = self.arg['k1']
    except:
      self.arg['k1'] = 0.45
      K1 = 0.45
    try:
      K2  = self.arg['k2']
    except:
      self.arg['k2'] = 2.8
      K2 = 2.8
    #print K1, K2
    psi = K1 * gor * (1.0 + sin(b)**2.0) / cos(b) * (1.0 - K1 * K2 * tan(b) * gor)
    return [[     1.0, 0.0,             0.0, 0.0, 0.0, 0.0],
            [tan(b)*r, 1.0,             0.0, 0.0, 0.0, 0.0],
            [     0.0, 0.0,             1.0, 0.0, 0.0, 0.0],
            [     0.0, 0.0, -tan(b - psi)*r, 1.0, 0.0, 0.0],
            [     0.0, 0.0,             0.0, 0.0, 1.0, 0.0],
            [     0.0, 0.0,             0.0, 0.0, 0.0, 1.0]]
        
  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m
  
  def IsAlive(self, x, y):
    return True

  def GetXBoundary(self):
    return 0.5 * self.arg['pw']

  def GetYBoundary(self):
    return 0.5 * self.arg['g']



# RF Gap
#########################################################################
class RFGap(OptBase):
  """
  No emittance growth !! for now !!
  """
  def __init__(self, **kwd):
    """
    ALREADY IMPLEMENTED -------------
    etl : E_0*T*L
    s : synchronous phase
    f : freq
    ---------------------------------
    TODO ----------------------------
    eg : emittance growth
    dW : energy gain flag 0 or 1
    h  : harmonics
    ---------------------------------
    """
    super(RFGap, self).__init__("rfgap", **kwd)
    self.visible = False
    self.arg['l'] = 0.0

  def R(self, l):
    h = 1.0 # harmonics
    lmda = c_light / self.arg['f']
    AA = pi * h * abs(self.beam.charge) * self.arg['etl'] * sin(self.arg['s'])
    BB = self.beam.rest_mass * lmda * self.beam.bet2
    
    x21 = -AA / (BB * self.beam.gam2 * self.beam.bg)
    y21 = x21
    z21 = 2.0 * AA / (BB * self.beam.bg)
        
    return [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [x21, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, y21, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, z21, 1.0]]
        
  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m
  
  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False
  
  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']



# Collimators
#########################################################################
class CollimatorX(OptBase):
  """
  Collimator X
  """
  def __init__(self, **kwd):
    """
    x : x half-gap
    """
    super(CollimatorX, self).__init__("colx", **kwd)
    self.visible = False
    self.arg['l'] = 0.0

  def R(self, l):
    return gen_identity_matrix()
        
  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m
  
  def IsAlive(self, x, y):
    if (-self.arg['x'] < x) and (x < self.arg['x']):
      return True
    else:
      return False

  def GetXBoundary(self):
    return self.arg['x']
  
  def GetYBoundary(self):
    return 0.0
 

class CollimatorY(OptBase):
  """
  Collimator Y
  """
  def __init__(self, **kwd):
    """
    y : y half-gap
    """
    super(CollimatorY, self).__init__("coly", **kwd)
    self.visible = False
    self.arg['l'] = 0.0

  def R(self, l):
    return gen_identity_matrix()
        
  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (-self.arg['y'] < y) and (y < self.arg['y']):
      return True
    else:
      return False

  def GetXBoundary(self):
    return 0.0

  def GetYBoundary(self):
    return self.arg['y']
  

# MONITOR
#########################################################################
class Monitor(OptBase):
  def __init__(self, **kwd):
    """
    l : length
    r : beam pipe inner radius
    """
    super(Monitor, self).__init__("monitor", **kwd)
    self.visible = True
    
  
  def R(self, l):
    return [[1.0,   l, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0,   l, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, l/self.beam.gam2],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

  def GetNextPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., self.arg['l'])

  def GetCenterPosition(self, pre_vec, rot_m):
    return pre_vec + rot_m * ThreeVector(0., 0., 0.5 * self.arg['l'])
    
  def GetNextRotMatrix(self, rot_m):
    return rot_m

  def IsAlive(self, x, y):
    if (x**2.0 + y**2.0) < self.arg['r']**2.0: 
      return True
    else: 
      return False
  
  def GetXBoundary(self):
    return self.arg['r']
  
  def GetYBoundary(self):
    return self.arg['r']
    






# Optic
#########################################################################
class Optic:
  """
  Beam optics calculation class
  """
  dfunc = {'drift'         : Drift,
           'quadrupole'    : Quadrupole,
           'solenoid'      : Solenoid,
           'dipole'        : Dipole,
           'rfgap'         : RFGap,
           'quadrupoleasym': QuadrupoleAsym,
           'thinlens'      : ThinLens,
           'colx'          : CollimatorX,
           'coly'          : CollimatorY,
           'monitor'       : Monitor}
  zleng = 0.0
  bmaxx = 0.0
  bmaxy = 0.0
  def __init__(self, ibeam):
    """
    ibeam : Beam instance, which means initial beam information
    """
    self.beam = ibeam
    self.devs = []
    self.v_init = ThreeVector(0., 0., 0.)

  def get_dzmax(self):
    return self.zleng/200
 
  def reset(self):
    self.devs = []
    self.zleng = 0.0
    
  def add_by_devname(self, devname, pset):
    f = self.dfunc[devname]
    self.add(f(**pset))
    
  def add(self, odif):
    """
    Adding optical devices
      odif : Optical device function
    """
    if odif.type == 'dipole':
        
      devf = DipoleEdge(b = odif.arg['ef'], 
                        r = odif.arg['r'], 
                        g = odif.arg['g'],
                        pw = odif.arg['pw'],
                        k1 = odif.arg['k1'],
                        k2 = odif.arg['k2'])
      devf.beam = self.beam
      self.devs.append(devf)
      
      odif.beam = self.beam
      self.devs.append(odif)
      
      devb = DipoleEdge(b = odif.arg['eb'], 
                        r = odif.arg['r'], 
                        g = odif.arg['g'],
                        pw = odif.arg['pw'],
                        k1 = odif.arg['k1'],
                        k2 = odif.arg['k2'])
      devb.beam = self.beam
      self.devs.append(devb)
      
    else:
      odif.beam = self.beam
      self.devs.append(odif)

    self.zleng += odif.arg['l']
    self.bmaxx = max(self.bmaxx, odif.GetXBoundary())
    self.bmaxy = max(self.bmaxy, odif.GetYBoundary())
    
      
  def get_ellipse_U(self):
    """
    with unit conversion!
    """
    m_1 = 0.001 # 1/m unit of gamma
    
    Dm = [self.beam.dx, self.beam.dpx,
          self.beam.dy, self.beam.dpy,
          0.0, 1.0]
    Sm = self.beam.get_sigma_matrix()
    zz = [0.0]      # z-length array

    ax = [self.beam.ax] # alpha
    ay = [self.beam.ay]
    bx = [self.beam.bx/m] # beta
    by = [self.beam.by/m]
    gx = [self.beam.gx/m_1] # gamma
    gy = [self.beam.gy/m_1]
    dx = [self.beam.dx/m] # dispersion
    dy = [self.beam.dy/m]
    dpx = [self.beam.dpx] # dispersion prime
    dpy = [self.beam.dpy]
    xx = [self.beam.x/mm] # envelope
    yy = [self.beam.y/mm]
    xp = [self.beam.xp/mrad] # divergence
    yp = [self.beam.yp/mrad]
    Rm = gen_identity_matrix()
    for c in self.devs:
      if c.arg['l'] == 0.0:
        if c.type.startswith('dri') or c.type.startswith('mon'): continue
        RMi = c.R(1.0)
        Rm = dot(Rm, RMi)
        Sm = dot(dot(RMi, Sm), tr(RMi))
        Dm = dot(RMi, Dm)
        
        # there are two points...
        ax[-1] = -Sm[0, 1]
        ay[-1] = -Sm[2, 3]
        bx[-1] = Sm[0, 0]/m
        by[-1] = Sm[2, 2]/m
        gx[-1] = Sm[1, 1]/m_1
        gy[-1] = Sm[3, 3]/m_1
        dx[-1] = Dm[0] /m
        dy[-1] = Dm[2] /m
        dpx[-1] = Dm[1]
        dpy[-1] = Dm[3]
        xx[-1] = sqrt(Sm[0, 0]*self.beam.ex) /mm
        yy[-1] = sqrt(Sm[2, 2]*self.beam.ey) /mm
        xp[-1] = sqrt(Sm[1, 1]*self.beam.ex) /mrad
        yp[-1] = sqrt(Sm[3, 3]*self.beam.ey) /mrad

      else:
        np = ceil(c.arg['l'] / self.get_dzmax())
        ds = c.arg['l'] / np
        inp = int(np)
        RMi = c.R(ds)
        for x in xrange(inp):
          Sm = dot(dot(RMi, Sm), tr(RMi))
          Dm = dot(RMi, Dm)
          zz.append(zz[-1] + ds /m)
          ax.append(-Sm[0, 1])
          ay.append(-Sm[2, 3])
          bx.append( Sm[0, 0]/m)
          by.append( Sm[2, 2]/m)
          gx.append( Sm[1, 1]/m_1)
          gy.append( Sm[3, 3]/m_1)
          dx.append(Dm[0] /m)
          dy.append(Dm[2] /m)
          dpx.append(Dm[1])
          dpy.append(Dm[3])
          xx.append(sqrt(Sm[0, 0]*self.beam.ex) /mm)
          yy.append(sqrt(Sm[2, 2]*self.beam.ey) /mm)
          xp.append(sqrt(Sm[1, 1]*self.beam.ex) /mrad)
          yp.append(sqrt(Sm[3, 3]*self.beam.ey) /mrad)
          
    return {'z': zz,
            'ax': ax,   'ay': ay,
            'bx': bx,   'by': by,
            'gx': gx,   'gy': gy,
            'dx': dx,   'dy': dy,
            'dpx': dpx, 'dpy': dpy,
            'x' : xx,   'y' : yy,
            'xp' : xp,   'yp' : yp
            }
  
  def get_tracks(self, tracks):
    """
    tracks[0] = (0.0, x, xp, y, yp, zp, alive_mask) ==> bullet
    tracks[i] = (z_i, x, xp, y, yp, zp, alive_mask) ==> ith
    """
    NPARTICLE = len(tracks[0][1])
    c_thin = None
    zz  = tracks[0][0]
    xo  = tracks[0][1]
    xpo = tracks[0][2]
    yo  = tracks[0][3]
    ypo = tracks[0][4]
    zpo = tracks[0][5]
    amo = tracks[0][6]
    # loop over the device matrixes
    
    NDEVICES = len(self.devs)
    #print NDEVICES
    dialog = Progress("Calculation Status", "Optical Tracking Processing...", NDEVICES, 
                      style=PD_ELAPSED_TIME | PD_AUTO_HIDE)
    for iii, c in enumerate(self.devs):
      if c.arg['l'] == 0.0:
        Ri = c.R(1.0)
        # loop over particles 
        for i in xrange(NPARTICLE):
          is_alive = amo[i]
          if is_alive: is_alive = c.IsAlive(xo[i], yo[i])
                
          if is_alive:
            X0 = (xo[i], xpo[i], yo[i], ypo[i], 0.0, zpo[i])
            X = dot(Ri, X0)
            xo[i] = X[0]
            xpo[i] = X[1]
            yo[i] = X[2]
            ypo[i] = X[3]
            zpo[i] = X[5]
          else:
            amo[i] = False
            xo[i] = X[0]
            xpo[i] = X[1]
            yo[i] = X[2]
            ypo[i] = X[3]
            zpo[i] = X[5]
        #done if l == 0.
      else:
        np = ceil(c.arg['l'] / self.get_dzmax())
        dz = c.arg['l'] / np
        inp = int(np)
        Ri = c.R(dz)
        # loop over each steps
        for i_unused in xrange(inp):
          x, xp, y, yp, zp, alive_mask = [], [], [], [], [], []
          # loop over particles 
          for i in xrange(NPARTICLE):
            is_alive = amo[i]
            if is_alive: is_alive = c.IsAlive(xo[i], yo[i])
                
            if is_alive:
              X0 = (xo[i], xpo[i], yo[i], ypo[i], 0.0, zpo[i])
              X = dot(Ri, X0)
              alive_mask.append(True)
              x.append(X[0])
              xp.append(X[1])
              y.append(X[2])
              yp.append(X[3])
              zp.append(X[5])
            else:
              alive_mask.append(False)
              x.append(0.0)
              xp.append(0.0)
              y.append(0.0)
              yp.append(0.0)
              zp.append(0.0)
            #done.
            
          zz += dz 
          tracks.append([zz, x, [p/mrad for p in xp], y, [p/mrad for p in yp], zp, alive_mask])
          xo, xpo, yo, ypo, zpo, amo = x, xp, y, yp, zp, alive_mask
      dialog.Update(iii+1)
    dialog.Destroy()
  
  def setInitPosition(self, v):
    self.v_init = v
  
  def get_result(self):
    """
    
    """
    m_1 = 0.001 # 1/m unit of gamma
    
    #DmOld = [self.beam.dx, self.beam.dpx,
    #         self.beam.dy, self.beam.dpy,
    #         0.0, 1.0]
    Dm = [self.beam.dx, self.beam.dpx,
          self.beam.dy, self.beam.dpy,
          0.0, 1.0]
    Rm = gen_identity_matrix()
    Sm = self.beam.get_sigma_matrix()
    for c in self.devs:
      #DmOld = [e for e in Dm]
      RMi = c.R(c.arg['l'])
      #Rm = dot(Rm, RMi)
      Rm = dot(RMi, Rm)
      Sm = dot(dot(RMi, Sm), tr(RMi))
      Dm = dot(RMi, Dm)
      
    
    return {'ax': -Sm[0, 1],     'ay': -Sm[2, 3],
            'bx':  Sm[0, 0]/m,   'by':  Sm[2, 2]/m,   # /meter 
            'gx':  Sm[1, 1]/m_1, 'gy':  Sm[3, 3]/m_1, # *meter
            'dx':  Dm[0]/m,      'dy':  Dm[2]/m,      # /meter
            'dpx': Dm[1],        'dpy':  Dm[3],
            'x' :  sqrt(Sm[0, 0]*self.beam.ex), 
            'y' :  sqrt(Sm[2, 2]*self.beam.ey),
            'xp' :  sqrt(Sm[1, 1]*self.beam.ex)/mrad, 
            'yp' :  sqrt(Sm[3, 3]*self.beam.ey)/mrad
            }
    
  def get_result_u(self):
    """
    유닛 문제가 있을듯 하여...
    """
    m_1 = 0.001 # 1/m unit of gamma
    
    #DmOld = [self.beam.dx, self.beam.dpx,
    #         self.beam.dy, self.beam.dpy,
    #         0.0, 1.0]
    Dm = [self.beam.dx, self.beam.dpx,
          self.beam.dy, self.beam.dpy,
          0.0, 1.0]
    Rm = gen_identity_matrix()
    Sm = self.beam.get_sigma_matrix()
    for c in self.devs:
      #DmOld = [e for e in Dm]
      RMi = c.R(c.arg['l'])
      #Rm = dot(Rm, RMi)
      Rm = dot(RMi, Rm)
      Sm = dot(dot(RMi, Sm), tr(RMi))
      Dm = dot(RMi, Dm)
      
    
    return {'ax': -Sm[0, 1], 'ay': -Sm[2, 3],
            'bx':  Sm[0, 0], 'by':  Sm[2, 2], # /meter 
            'gx':  Sm[1, 1], 'gy':  Sm[3, 3], # *meter
            'dx':  Dm[0],    'dy':  Dm[2],    # /meter
            'dpx': Dm[1],    'dpy':  Dm[3],
            'x' :  sqrt(Sm[0, 0]*self.beam.ex), 
            'y' :  sqrt(Sm[2, 2]*self.beam.ey),
            'xp' :  sqrt(Sm[1, 1]*self.beam.ex), 
            'yp' :  sqrt(Sm[3, 3]*self.beam.ey)
            }

  def get_all_result(self):
    """
    
    """
    Dm = [self.beam.dx, self.beam.dpx,
          self.beam.dy, self.beam.dpy,
          0.0, 1.0]
    Rm = gen_identity_matrix()
    Sm = self.beam.get_sigma_matrix()
    s = 0.0
    Sm = dot(Sm, gen_identity_matrix())
    allr = []
    allr.append({'type': 'start',
                 'ax': -Sm[0, 1],     'ay': -Sm[2, 3],
                 'bx':  Sm[0, 0]/m,   'by':  Sm[2, 2]/m,   # /meter 
                 'dx':  Dm[0]/m,      'dy':  Dm[2]/m,      # /meter
                 'x' :  sqrt(Sm[0, 0]*self.beam.ex), 
                 'y' :  sqrt(Sm[2, 2]*self.beam.ey),
                 's' :  0.0}
                )
    
    
    for c in self.devs:
      RMi = c.R(c.arg['l'])
      Rm = dot(RMi, Rm)
      Sm = dot(dot(RMi, Sm), tr(RMi))
      Dm = dot(RMi, Dm)
      s += c.arg['l']
      
      allr.append({'type': c.type,
                   'ax': -Sm[0, 1],     'ay': -Sm[2, 3],
                   'bx':  Sm[0, 0]/m,   'by':  Sm[2, 2]/m,   # /meter 
                   'dx':  Dm[0]/m,      'dy':  Dm[2]/m,      # /meter
                   'x' :  sqrt(Sm[0, 0]*self.beam.ex), 
                   'y' :  sqrt(Sm[2, 2]*self.beam.ey)}
                  )
      allr[-1]['s'] = s/m
    return allr

  #def get_maximum_ellipse(self):
    


