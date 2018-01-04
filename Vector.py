#!/usr/bin/env python
# -*- coding: utf-8 -*-
#===============================================================================
#
#  2D/3D Vectors
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================

from math import sqrt, sin, cos, acos, tan, pi
from exceptions import TypeError

class ThreeVector:
  def __init__(self, x, y, z):
    self.x = x
    self.y = y
    self.z = z
    
  def __add__(self, o):
    return ThreeVector(self.x + o.x, self.y + o.y, self.z + o.z)
  
  def __sub__(self, o):
    return ThreeVector(self.x - o.x, self.y - o.y, self.z - o.z)
    
  def __mul__(self, o):
    try:
      return self.x * o.x + self.y * o.y + self.z * o.z
    except:
      return ThreeVector(self.x * o, self.y * o, self.z * o)
  
  def __div__(self, o):
    return ThreeVector(self.x / o, self.y / o, self.z / o)
  
  def cross(self, o):
    return ThreeVector(self.y * o.z - self.z * o.y,
                       self.z * o.x - self.x * o.z,
                       self.x * o.y - self.y * o.x)
  def dot(self, o):
    return self.x * o.x + self.y * o.y + self.z * o.z
  
  def __pow__(self, o):
    return self.x ** o + self.y ** o + self.z ** o
  
  def __rmul__(self, o):
    return ThreeVector(self.x * o, self.y * o, self.z * o)
  
  def __iadd__(self, o):
    self.x += o.x
    self.y += o.y
    self.z += o.z
  
  def __isub__(self, o):
    self.x -= o.x
    self.y -= o.y
    self.z -= o.z
  
  def __neg__(self):
    return ThreeVector(self.x * -1.0,
                       self.y * -1.0,
                       self.z * -1.0)
  
  def __pos__(self):
    pass
  
  def __abs__(self):
    return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
      
  def getTheta(self, o):
    return acos((self * o)/(abs(self) * abs(o)))
  
  def getBase(self):
    return self / abs(self)
  
  def rotate(self, o, ax, th):
    """
    o  : origin vector
    ax : rotate axis vector
    th : rotation angle
    """
    d = ax.getBase() # direction vector
    cc = cos(th)
    ss = sin(th)
    x, y, z = self.x, self.y, self.z
    a, b, c = o.x, o.y, o.z
    u, v, w = d.x, d.y, d.z
    uxvywz = u*x + v*y + w*z
    return ThreeVector(
      (a * (v**2 + w**2) - u * (b*v + c*w - uxvywz)) * (1 - cc) 
       + x * cc + (-c*v +b*w -w*y + v*z) * ss,
      (b * (u**2 + w**2) - v * (a*u + c*w - uxvywz)) * (1 - cc) 
       + y * cc + (c*u -a*w +w*x - u*z) * ss,
      (c * (u**2 + v**2) - w * (a*u + b*v - uxvywz)) * (1 - cc) 
       + z * cc + (-b*u +a*v -v*x + u*y) * ss)
  
  def rot_origin(self, ax, th):
    """
    ax : rotate axis vector along the origin
    th : rotation angle
    """
    return self.rotate(ThreeVector(0., 0., 0.), ax, th)
  
  
  def print_out(self):
    print "(%.3f, %.3f, %.3f)" % (self.x, self.y, self.z)
  
  def getTuple(self):
    return (self.x, self.y, self.z)
  
  def getXY(self):
    return TwoVector(self.x, self.y)
  
  def getYZ(self):
    return TwoVector(self.y, self.z)
  
  def getZX(self):
    return TwoVector(self.z, self.x)
  
  
  def __repr__(self):
    return '(%.6f, %.6f, %.6f)' % (self.x,  self.y, self.z)
  
  def new(self):
    return ThreeVector(self.x, self.y, self.z)

  def unit( self ):
    return self / abs(self)


xhat = ThreeVector(1., 0., 0.)
yhat = ThreeVector(0., 1., 0.)
zhat = ThreeVector(0., 0., 1.)


#==============================================================================
class TwoVector:
  def __init__(self, x, y):
    self.x = float(x)
    self.y = float(y)
    
  def __add__(self, o):
    return TwoVector(self.x + o.x, self.y + o.y)
  
  def __sub__(self, o):
    return TwoVector(self.x - o.x, self.y - o.y)
    
  def __mul__(self, o):
    return TwoVector(o * self.x, o * self.y)
  
  def __div__(self, o):
    return TwoVector(self.x / o, self.y / o)
  
  def __pow__(self, o):
    return self.x ** o + self.y ** o
  
  def __rmul__(self, o):
    return TwoVector(self.x * o, self.y * o)
  
  def __iadd__(self, o):
    self.x += o.x
    self.y += o.y
  
  def __isub__(self, o):
    self.x -= o.x
    self.y -= o.y
  
  def __neg__(self):
    return TwoVector(self.x * -1.0,
                       self.y * -1.0)
  
  def __pos__(self):
    pass
  
  def __abs__(self):
    return sqrt(self.x ** 2 + self.y ** 2)
      
  def dot(self, o):
    return self.x * o.x + self.y * o.y
  
  def getTheta(self, o):
    return acos(self.dot(o)/(abs(self) * abs(o)))
  
  def getBase(self):
    return self / abs(self)
  
  def get_rotation(self, th, **kwd):
    """
    th : rotation angle
    center = two vector center point
    
   """
    try:
      dv = TwoVector( self.x - kwd['center'].x, self.y - kwd['center'].y )
      cc = cos(th)
      ss = sin(th)
      dv.x = cc*dv.x - ss*dv.y
      dv.y = ss*dv.x + cc*dv.y
      dv = dv + kwd['center']
      return dv
  
    except:
      cc = cos(th)
      ss = sin(th)
      return TwoVector(cc*self.x - ss*self.y, ss*self.x + cc*self.y)
  
  def rotate(self, th, **kwd):
    """
    th : rotation angle
    center = two vector center point
    
    ex)
      v.rotate(pi, v2center)  
    
    """
    try:
      dv = TwoVector( self.x - kwd['center'].x, self.y - kwd['center'].y )
      cc = cos(th)
      ss = sin(th)
      dv.x = cc*dv.x - ss*dv.y
      dv.y = ss*dv.x + cc*dv.y
      self = dv + kwd['center']
  
    except:
      cc = cos(th)
      ss = sin(th)
      self.x, self.y = cc*self.x - ss*self.y, ss*self.x + cc*self.y
      
    return self
  
  def print_out(self):
    print "(%.3f, %.3f)" % (self.x, self.y)
  
  def getTuple(self):
    return (self.x, self.y)
      
  def __repr__(self):
    return '(%.6f, %.6f)' % (self.x,  self.y)
  
  def new(self):
    return TwoVector(self.x, self.y)
  
  def convXY3D(self, **kwd):
    try:
      return ThreeVector(self.x, self.y, kwd['z'])
    except:
      return ThreeVector(self.x, self.y, 0.0)
      
  def convYZ3D(self, **kwd):
    try:
      return ThreeVector(kwd['x'], self.x, self.y)
    except:
      return ThreeVector(0.0, self.x, self.y)

  def convZX3D(self, **kwd):
    try:
      return ThreeVector(self.y, kwd['y'], self.x)
    except:
      return ThreeVector(self.y, 0.0, self.x)
  
  def unit( self ):
    return self / abs(self)

  def get_mirror_yaxis( self ):
    return TwoVector(-self.x, self.y)

  def get_mirror_xaxis( self ):
    return TwoVector(self.x, -self.y)

xhat2D = TwoVector(1., 0.)
yhat2D = TwoVector(0., 1.)




def ToTwoVector( iterable ):
  if len(iterable) >= 2:
    return ThreeVector(iterable[0], iterable[1])
  else:
    return 0

def ToThreeVector( iterable ):
  if len(iterable) >= 3:
    return ThreeVector(iterable[0], iterable[1], iterable[2])
  else:
    return 0

def StrToTwoVector( str ):
  """
  :param str:  comma seperated pair numbers
  :return: TwoVector
  """
  str = str.replace('(', '')
  str = str.replace(')', '')
  slst = str.split( ',' )
  try:
    return TwoVector( float(slst[0].strip()) , float(slst[1].strip()) )
  except:
    return 0


def StrToThreeVector(str):
  """
  :param str:  comma seperated pair numbers
  :return: TwoVector
  """
  str = str.replace('(', '')
  str = str.replace(')', '')
  slst = str.split(',')
  try:
    return ThreeVector(float(slst[0].strip()),
                       float(slst[1].strip()),
                       float(slst[2].strip()))
  except:
    return 0


#==============================================================================
class RotationMatrix:
  def __init__(self, **kwd):
    """
    RotationMatrix() or
    RotationMatrix(R = R0)
    """
    try:
      self.xx = kwd['R'].xx
      self.xy = kwd['R'].xy
      self.xz = kwd['R'].xz
      self.yx = kwd['R'].yx
      self.yy = kwd['R'].yy
      self.yz = kwd['R'].yz
      self.zx = kwd['R'].zx
      self.zy = kwd['R'].zy
      self.zz = kwd['R'].zz
    except:
      self.xx = 1.
      self.xy = 0
      self.xz = 0
      self.yx = 0
      self.yy = 1.
      self.yz = 0
      self.zx = 0
      self.zy = 0
      self.zz = 1.
  
  def set(self, *args):
    self.xx = args[0]
    self.xy = args[1]
    self.xz = args[2]
    self.yx = args[3]
    self.yy = args[4]
    self.yz = args[5]
    self.zx = args[6]
    self.zy = args[7]
    self.zz = args[8]
    
  def __mul__(self, r):
    try: # rr = self * V
      return ThreeVector(
        self.xx * r.x + self.xy * r.y + self.xz * r.z,
        self.yx * r.x + self.yy * r.y + self.yz * r.z,
        self.zx * r.x + self.zy * r.y + self.zz * r.z)
    except:
      try: # rr = self * R
        rr = RotationMatrix()
        rr.set(
          self.xx * r.xx + self.xy * r.yx + self.xz * r.zx,
          self.xx * r.xy + self.xy * r.yy + self.xz * r.zy,
          self.xx * r.xz + self.xy * r.yz + self.xz * r.zz,
          self.yx * r.xx + self.yy * r.yx + self.yz * r.zx,
          self.yx * r.xy + self.yy * r.yy + self.yz * r.zy,
          self.yx * r.xz + self.yy * r.yz + self.yz * r.zz,
          self.zx * r.xx + self.zy * r.yx + self.zz * r.zx,
          self.zx * r.xy + self.zy * r.yy + self.zz * r.zy,
          self.zx * r.xz + self.zy * r.yz + self.zz * r.zz)
        return rr
      except: # rr = self * Const.
        rr = RotationMatrix()
        rr.set(
          self.xx * r, self.xy * r, self.xz * r,
          self.yx * r, self.yy * r, self.yz * r,
          self.zx * r, self.zy * r, self.zz * r)
        return rr
    
  def __rmul__(self, l):
    try: # rr = Const * self
      rr = RotationMatrix()
      rr.set(
        self.xx * l, self.xy * l, self.xz * l,
        self.yx * l, self.yy * l, self.yz * l,
        self.zx * l, self.zy * l, self.zz * l)
      return rr
    except:
      raise TypeError
  
  def __neg__(self):
    try: # rr = Const * self
      rr = RotationMatrix()
      rr.set(
        self.xx * -1.0, self.xy * -1.0, self.xz * -1.0,
        self.yx * -1.0, self.yy * -1.0, self.yz * -1.0,
        self.zx * -1.0, self.zy * -1.0, self.zz * -1.0)
      return rr
    except:
      raise TypeError

  def __pos__(self):
    pass
    
  def __repr__(self):
    form = '((%.2f, %.2f, %.2f),\n (%.2f, %.2f, %.2f),\n (%.2f, %.2f, %.2f))\n'
    return form % (self.xx, self.xy, self.xz,
                   self.yx, self.yy, self.yz,
                   self.zx, self.zy, self.zz)
                   
  def rotate(self, a, ax):
    """
    a : rotation angle
    ax : rotate axis vector
    """
    if (a != 0.0):
      ll = abs(ax)
      if (ll == 0.0):
        print "error"
      else:
        sa = sin(a)
        ca = cos(a)
        dx = ax.x/ll
        dy = ax.y/ll
        dz = ax.z/ll
        m = RotationMatrix()
        m.set(
          ca+(1-ca)*dx*dx,    (1-ca)*dx*dy-sa*dz, (1-ca)*dx*dz+sa*dy,
          (1-ca)*dy*dx+sa*dz, ca+(1-ca)*dy*dy,    (1-ca)*dy*dz-sa*dx,
          (1-ca)*dz*dx-sa*dy, (1-ca)*dz*dy+sa*dx,  ca+(1-ca)*dz*dz)
        self.set(
          self.xx * m.xx + self.xy * m.yx + self.xz * m.zx,
          self.xx * m.xy + self.xy * m.yy + self.xz * m.zy,
          self.xx * m.xz + self.xy * m.yz + self.xz * m.zz,
          self.yx * m.xx + self.yy * m.yx + self.yz * m.zx,
          self.yx * m.xy + self.yy * m.yy + self.yz * m.zy,
          self.yx * m.xz + self.yy * m.yz + self.yz * m.zz,
          self.zx * m.xx + self.zy * m.yx + self.zz * m.zx,
          self.zx * m.xy + self.zy * m.yy + self.zz * m.zy,
          self.zx * m.xz + self.zy * m.yz + self.zz * m.zz)
  
  def get_rotation(self, a, ax):
    """
    returns new instance
    """
    rr = RotationMatrix(R = self)
    rr.rotate(a, ax)
    return rr
  

  


