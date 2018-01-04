#!/usr/bin/env python
# -*- coding: utf-8 -*-
#===============================================================================
#
#  Physics Units
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================

millimeter  = 1.                        
millimeter2 = millimeter*millimeter
millimeter3 = millimeter*millimeter*millimeter

centimeter  = 10.*millimeter   
centimeter2 = centimeter*centimeter
centimeter3 = centimeter*centimeter*centimeter

meter  = 1000.*millimeter                  
meter2 = meter*meter
meter3 = meter*meter*meter

kilometer = 1000.*meter                   
kilometer2 = kilometer*kilometer
kilometer3 = kilometer*kilometer*kilometer

parsec = 3.0856775807e+16*meter

micrometer = 1.e-6 *meter             
nanometer  = 1.e-9 *meter
angstrom   = 1.e-10*meter
fermi      = 1.e-15*meter

barn      = 1.e-28*meter2
millibarn = 1.e-3 *barn
microbarn = 1.e-6 *barn
nanobarn  = 1.e-9 *barn
picobarn  = 1.e-12*barn

# symbols
nm  = nanometer                        
um  = micrometer                        

mm  = millimeter                        
mm2 = millimeter2
mm3 = millimeter3

cm  = centimeter   
cm2 = centimeter2
cm3 = centimeter3

m  = meter                  
m2 = meter2
m3 = meter3

km  = kilometer                   
km2 = kilometer2
km3 = kilometer3

pc = parsec


# Angle

radian      = 1.                  
milliradian = 1.e-3*radian
degree      = (3.14159265358979323846/180.0)*radian

steradian   = 1.
  
# symbols
rad  = radian
mrad = milliradian
sr   = steradian
deg  = degree

#
# Time [T]
#
nanosecond  = 1.
second      = 1.e+9 *nanosecond
millisecond = 1.e-3 *second
microsecond = 1.e-6 *second
picosecond  = 1.e-12*second
minute      = 60.   *second

hertz     = 1./second
Hz        = hertz
kilohertz = 1.e+3*hertz
megahertz = 1.e+6*hertz
gigahertz = 1.e+9*hertz

# symbols

ns = nanosecond
s  = second
ms = millisecond
us = microsecond
#ps = picosecond

kHz = kilohertz
MHz = megahertz
GHz = gigahertz

#
# Electric charge [Q]
#
eplus = 1. # positron charge
e_SI  = 1.602176487e-19 # positron charge in coulomb
coulomb = eplus/e_SI # coulomb = 6.24150 e+18 * eplus

#
# Energy [E]
#
megaelectronvolt = 1. 
electronvolt     = 1.e-6*megaelectronvolt
kiloelectronvolt = 1.e-3*megaelectronvolt
gigaelectronvolt = 1.e+3*megaelectronvolt
teraelectronvolt = 1.e+6*megaelectronvolt
petaelectronvolt = 1.e+9*megaelectronvolt

joule = electronvolt/e_SI  # joule = 6.24150 e+12 * MeV

calorie = 0.239005736 * joule  # 1 joule = 0.239005736 calorie
cal = calorie
kcal = calorie * 1000.0

# symbols
MeV = megaelectronvolt
eV = electronvolt
keV = kiloelectronvolt
GeV = gigaelectronvolt
TeV = teraelectronvolt
PeV = petaelectronvolt


#
# Potentials [E/Q]
#
# symbols
MV = 1.
V = 1.e-6
kV = 1.e-3

#
# Mass [E][T^2][L^-2]
#
kilogram = joule*second*second/(meter*meter)   
gram = 1.e-3*kilogram
milligram = 1.e-3*gram

# symbols
kg = kilogram
g = gram
mg = milligram

#
# Power [E][T^-1]
#
watt = joule/second# watt = 6.24150 e+3 * MeV/ns
kW = 1.e3 * watt
MW = 1.e6 * watt


#
# Force [E][L^-1]
#
newton = joule/meter# newton = 6.24150 e+9 * MeV/mm

#
# Pressure [E][L^-3]
#
#define pascal hep_pascal # a trick to avoid warnings 
hep_pascal = newton/m2 #pascal = 6.24150 e+3 * MeV/mm3
pascal = hep_pascal
bar        = 100000*pascal # bar    = 6.24150 e+8 * MeV/mm3
atmosphere = 101325*pascal # atm    = 6.32420 e+8 * MeV/mm3

#
# Electric current [Q][T^-1]
#
ampere      = coulomb/second # ampere = 6.24150 e+9 * eplus/ns
milliampere = 1.e-3*ampere
microampere = 1.e-6*ampere
nanoampere = 1.e-9*ampere

mA = milliampere
uA = microampere
nA = nanoampere

#
# Electric potential [E][Q^-1]
#
megavolt = megaelectronvolt/eplus
kilovolt = 1.e-3*megavolt
volt     = 1.e-6*megavolt

#
# Electric resistance [E][T][Q^-2]
#
ohm = volt/ampere # ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)

#
# Electric capacitance [Q^2][E^-1]
#
farad = coulomb/volt # farad = 6.24150e+24 * eplus/Megavolt
millifarad = 1.e-3*farad
microfarad = 1.e-6*farad
nanofarad = 1.e-9*farad
picofarad = 1.e-12*farad

#
# Magnetic Flux [T][E][Q^-1]
#
weber = volt*second # weber = 1000*megavolt*ns

#
# Magnetic Field [T][E][Q^-1][L^-2]
#
tesla     = volt*second/meter2 #tesla = 0.001*megavolt*ns/mm2

gauss     = 1.e-4*tesla
kilogauss = 1.e-1*tesla

#
# Inductance [T^2][E][Q^-2]
#
henry = weber/ampere# henry = 1.60217e-7*MeV*(ns/eplus)**2

#
# Temperature
#
kelvin = 1.

#
# Amount of substance
#
mole = 1.

#
# Activity [T^-1]
#
becquerel = 1./second 
curie = 3.7e+10 * becquerel

#
# Absorbed dose [L^2][T^-2]
#
gray = joule/kilogram 
kilogray = 1.e+3*gray
milligray = 1.e-3*gray
microgray = 1.e-6*gray

#
# Luminous intensity [I]
#
candela = 1.

#
# Luminous flux [I]
#
lumen = candela*steradian

#
# Illuminance [I][L^-2]
#
lux = lumen/meter2

#
# Miscellaneous
#
perCent     = 0.01 
perThousand = 0.001
perMillion  = 0.000001




## constants
pi  = 3.14159265358979323846
twopi  = 2*pi
halfpi  = pi/2
pi2 = pi*pi

#
# 
#
Avogadro = 6.02214179e+23/mole

#
# c   = 299.792458 mm/ns
# c^2 = 898.7404 (mm/ns)^2 
#
c_light   = 2.99792458e+8 * m/s
c_squared = c_light * c_light

#
# h     = 4.13566e-12 MeV*ns
# hbar  = 6.58212e-13 MeV*ns
# hbarc = 197.32705e-12 MeV*mm
#
h_Planck      = 6.62606896e-34 * joule*s
hbar_Planck   = h_Planck/twopi
hbarc         = hbar_Planck * c_light
hbarc_squared = hbarc * hbarc

#
#
#
electron_charge = - eplus # see SystemOfUnits.h
e_squared = eplus * eplus

#
# amu_c2 - atomic equivalent mass unit
#        - AKA, unified atomic mass unit (u)
# amu    - atomic mass unit
#
electron_mass_c2 = 0.510998910 * MeV
proton_mass_c2 = 938.272013 * MeV
neutron_mass_c2 = 939.56536 * MeV
amu_c2 = 931.494028 * MeV
amu = amu_c2/c_squared

#
# permeability of free space mu0    = 2.01334e-16 Mev*(ns*eplus)^2/mm
# permittivity of free space epsil0 = 5.52636e+10 eplus^2/(MeV*mm)
#
mu0      = 4*pi*1.e-7 * henry/m
epsilon0 = 1./(c_squared*mu0)

#
# electromagnetic coupling = 1.43996e-12 MeV*mm/(eplus^2)
#
elm_coupling           = e_squared/(4*pi*epsilon0)
fine_structure_const   = elm_coupling/hbarc
classic_electr_radius  = elm_coupling/electron_mass_c2
electron_Compton_length = hbarc/electron_mass_c2
Bohr_radius = electron_Compton_length/fine_structure_const

alpha_rcl2 = fine_structure_const *classic_electr_radius *classic_electr_radius

twopi_mc2_rcl2 = twopi*electron_mass_c2 *classic_electr_radius *classic_electr_radius
#
#
#
k_Boltzmann = 8.617343e-11 * MeV/kelvin

#
#
#
STP_Temperature = 273.15*kelvin
STP_Pressure    = 1.*atmosphere
kGasThreshold   = 10.*mg/cm3

#
#
#
universe_mean_density = 1.e-25*g/cm3
