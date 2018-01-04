#!/usr/bin/env python
# -*- coding: utf-8 -*-

#===============================================================================
#
#  OpticsExpert entry point
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================

optics_expert_title = 'Optics Expert v0.1.1'

# platform check
from platform import python_version_tuple as pyversion
PYVER = pyversion()[0]

if PYVER is not '2':
  print 'please install python 2.7'
  exit()

from Widgets import *

app = wx.App(False)
gm = GMixer(optics_expert_title)
app.MainLoop()


