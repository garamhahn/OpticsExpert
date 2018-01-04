#!/usr/bin/env python
# -*- coding: utf-8 -*-

#===============================================================================
#
#  Windows widgets
#
#                                                 2017. 11. 18. by Garam Hahn
#
#===============================================================================



# ===================================================================================================
# Imports
# ===================================================================================================

# 1. wxWidget python module for window widget instead of tcl/tk
import wx

# 2. Matplotlib for plotting
import matplotlib
matplotlib.use('WXAgg')

import matplotlib.lines as mlines

from matplotlib.patches import Rectangle as rect
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_wxagg import \
                                       FigureCanvasWxAgg as FigCanvas, \
                                       NavigationToolbar2WxAgg as NavigationToolbar

# 3. python basic modules
from sys import maxint as MAXINT, argv
from os import getcwd, path as ospath
from json import dump, load
from random import uniform
from scipy.optimize import minimize, fmin

# 4. Linear optics modules
from LinearOptics import *

# 5. others
from datetime import datetime
from codecs import open as copen
from platform import system as plsys



# ===================================================================================================
# Global variables
# ===================================================================================================
nTICKS = 1000
nfTICKS = float(nTICKS)
fTICKS = 1./nfTICKS
f2s = lambda f : '%.6f' % f

ICONPATH = ospath.abspath(ospath.dirname(__file__))+'/icons/'

#color maps
CLX = ['#ff0000', '#ea5300', '#f49b00', '#fdcc00', '#b3d465']
CLY = ['#0000ff', '#0067b3', '#00a0e9', '#7ecef4', '#8957a1']

# donot change ITC, it is linked to histo1d class
ITC = ['#ccccff', '#00ff00', '#44ff00', '#88ff00', '#ccff00', '#ffff00', '#ffcc00', '#ff8800', '#ff4400', '#ff0000']
#ITC = ['#ff0000', '#ff4400', '#ff8800', '#ffcc00', '#ffff00', '#ccff00', '#88ff00', '#44ff00', '#00ff00', '#ccccff']

c_table = {'dipole'    : '#333333',
           'D'         : '#333333',
           'quadrupole': '#FF0000',
           'QF'        : '#FF0000',
           'QD'        : '#0000ff',
           'solenoid'  : '#00AA00',
           'S'         : '#00AA00',
           'curr_pos'  : '#ffcc00',
           'monitor'   : '#d0d0ff',
           'M'         : '#d0d0ff'}

OSNAME = plsys()



# ===================================================================================================
class GMixer(wx.Frame):
    optviewer = None
    trkviewer = None
    filepath = ''
    filedir = ''
    filename = ''
    beam = []
    cutbullet = None

    def __init__(self, title):
        wx.Frame.__init__(self, None, 1, title)
        self.title_name = title

        # dummy fit value
        self.fitv = fitvalue_ab(0., 1. * m, 0., 1. * m)

        # dummy beam
        b = BeamSpec(11174.67 * MeV,
                     6.0,
                     400.0 * 12.0 * MeV)
        b.set_xxp(1 * mm, 1 * mrad)
        b.set_yyp(1 * mm, 1 * mrad)
        b.set_dpop(0.125 * perCent)
        self.beam.append(b)

        self.initStatusBar()
        self.createLists()  # self.reflc, self.comlc
        self.createPlot()

        self.menuData = [
            ("&File",
             (("&New", "New Optics", self.OnNew),
              ("&Open", "Open Optics Set-up File", self.OnOpen),
              ("&Save", "Save Optics Result File", self.OnSave),
              ("Save As", "Save As A New Result File", self.OnSaveAs),
              ("", "", ""),
              ("Export",
               (("TRACE3D", "Export to trace3d input", self.OnExportT3D),
                ("MAD", "Export to mad input", self.OnExportMAD),
                ("", "", ""),
                ("LaTeX", "Export to LaTeX document", self.OnExportTEX),
                ("", "", ""),
                ("Excel", "Export to EXCEL CSV", self.OnExportCSV)
                )),
              ("", "", ""),
              ("Append", "Append Other Optics Here", self.OnAppend),
              ("", "", ""),
              ("About...", "About This Program", self.OnAbout),
              ("", "", ""),
              ("&Quit", "Quit", self.OnQuit)
              )),
            ("&Beam",
             (("&Twiss Ellipse Boundary", "Initial Twiss Beam Information", self.OnSetTwissBeam),
              ("&Arbitrary Distribution(for tracking only)", "Arbitrary Input Beam Generator", self.OnSetArbitraryBeam),
              ("", "", ""),
              ("&Chagne particle with holding optics", "Optical component scaling", self.OnParticleChange),
              ("", "", ""),
              ("&Inverse tracing", "Inverse the beam line", self.OnInverse)
              )),
            ("&Result View",
             (("&Optic", "Optic Functions", self.OnShowOpticsWin),
              ("&Tracks", "Optic Tracking", self.OnShowTrackWin),
              ("&Layout", "Layout Viewer", self.OnShowLayoutWin)
              ))]
        self.createMenuBar()

        self.toolbarData = (
            ("", "", "", ""),
            ("New", "new.png", "Create new sketch", self.OnNew),
            ("", "", "", ""),
            ("Open", "open.png", "Open existing sketch", self.OnOpen),
            ("Save", "save.png", "Save existing sketch", self.OnSave),
            ("", "", "", ""),
            ("Add Reference", "addref.png", "Add reference item", self.reflc.OnAddItem),
            ("Del Reference", "delref.png", "Delete reference item", self.reflc.OnRemoveItem),
            ("Modify Reference", "modifyref.png", "Modify reference item", self.reflc.OnModifyItem),
            ("", "", "", ""),
            ("Push Up Item", "pushupitem.png", "Push up selected item", self.comlc.OnPushupItem),
            ("Push Down Item", "pushdownitem.png", "Push down selected item", self.comlc.OnPushdownItem),
            ("", "", "", ""),
            ("Fits", "fit.png", "Fits", self.OnFit))
        self.createToolBar()

        # print argv

        if len(argv) > 2:
            self.filepath = argv[1]

            # read file
            try:
                f = open(self.filepath)
                data = load(f)
                f.close()
            except:
                print 'data reading error'
                return

            try:
                self.beam = []
                for b in data['beam']:
                    ibeam = BeamSpec(b['rest_mass'], b['charge'], b['kinetic_energy'])
                    ibeam.set_xabe(b['ax'], b['bx'], b['ex'])
                    ibeam.set_yabe(b['ay'], b['by'], b['ey'])
                    ibeam.dpop = b['dpop']
                    self.beam.append(ibeam)
            except:
                print 'beam loading error'
                return

            try:
                self.comlc.clist = data['clist']
                self.reflc.rlist = data['rlist']
                self.comlc.lastuid = int(data['clastuid'])
                self.reflc.lastuid = int(data['rlastuid'])
            except:
                print 'reference and component list loading error'
                return

            try:
                self.fitv = fitvalue_ab(data['fitv']['ax'],
                                        data['fitv']['bx'],
                                        data['fitv']['ay'],
                                        data['fitv']['by'])

            except:
                print 'target fit value loading error'

            self.reflc.ReGenerateList()
            self.comlc.ReGenerateList()

            self.Plot()
            self.SetTitle(self.filename)

        self.Show()

    def initStatusBar(self):
        self.statusbar = self.CreateStatusBar()
        self.statusbar.SetFieldsCount(3)
        self.statusbar.SetStatusWidths([-1, -1, -1])

    def createPlot(self):
        self.optviewer = GOpticsView()
        self.optviewer.Show()
        self.UpdateStatusbar('Application Started')
        self.optview_status = True

    def createLists(self):
        """ """
        self.reflc = GReferenceList(self)
        self.comlc = GComponentList(self, self.reflc)  # runtime problem, assigned dually
        self.reflc.comlc = self.comlc

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.reflc, 0, wx.EXPAND)
        hbox.Add(self.comlc, 1, wx.EXPAND)
        self.SetSizer(hbox)

    def createToolBar(self):
        toolbar = self.CreateToolBar()
        for each in self.toolbarData:
            self.createSimpleTool(toolbar, *each)
        toolbar.Realize()

    def createSimpleTool(self, toolbar, label, filename, hlp, handler):
        if not label:
            toolbar.AddSeparator()
            return
        bmp = wx.Bitmap(ICONPATH + filename, wx.BITMAP_TYPE_PNG)
        # tool = toolbar.AddTool(-1, bmp, label, hlp)  #
        tool = toolbar.AddTool(-1, label, bmp, hlp)
        self.Bind(wx.EVT_MENU, handler, tool)

    def createMenuBar(self):
        menuBar = wx.MenuBar()
        for eachMenuData in self.menuData:
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1]
            menuBar.Append(self.createMenu(menuItems), menuLabel)
        self.SetMenuBar(menuBar)

    def createMenu(self, menuData):
        menu = wx.Menu()
        for eachItem in menuData:
            if len(eachItem) == 2:
                label = eachItem[0]
                subMenu = self.createMenu(eachItem[1])
                menu.Append(wx.NewId(), label, subMenu)
            else:
                self.createMenuItem(menu, *eachItem)
        return menu

    def createMenuItem(self, menu, label, status, handler, kind=wx.ITEM_NORMAL):
        if not label:
            menu.AppendSeparator()
            return
        menuItem = menu.Append(-1, label, status, kind)
        self.Bind(wx.EVT_MENU, handler, menuItem)

    def Show(self):
        self.Fit()
        super(wx.Frame, self).Show()

    def GetOptList(self):
        ol = []
        for b in self.beam:
            o = Optic(b)
            for ci in self.comlc.clist:
                ri = self.reflc.ru(ci['ruid'])
                rfo = self.reflc.rform[ri['type']]
                pset = dict(**ri['pset'])
                pset[rfo['cval']['name']] = ci['cval']
                o.add_by_devname(ri['type'], pset)
            ol.append(o)
        return ol

    def GetFirstOptics(self):
        o = Optic(self.beam[0])
        # print self.comlc.clist
        for ci in self.comlc.clist:
            ri = self.reflc.ru(ci['ruid'])
            rfo = self.reflc.rform[ri['type']]
            pset = dict(**ri['pset'])
            pset[rfo['cval']['name']] = ci['cval']
            o.add_by_devname(ri['type'], pset)
            # print ri['type'], pset
        return o

    def FitFunc(self, x):
        """
    TRACE-3D Method
    """
        # print "try -> ", x
        o = Optic(self.beam[0])
        j = 0
        for ci in self.comlc.clist:
            if ci['tune']:
                ci['cval'] = x[j]
                j = j + 1
            ri = self.reflc.ru(ci['ruid'])
            rfo = self.reflc.rform[ri['type']]
            pset = dict(**ri['pset'])
            pset[rfo['cval']['name']] = ci['cval']

            o.add_by_devname(ri['type'], pset)
            if ci['mp']: break

        res = o.get_result_u()
        # res = o.get_result()

        # mismatching factor calculation
        rx = self.fitv.bx * res['gx'] \
             + self.fitv.gx * res['bx'] \
             - 2.0 * self.fitv.ax * res['ax']
        ry = self.fitv.by * res['gy'] \
             + self.fitv.gy * res['by'] \
             - 2.0 * self.fitv.ay * res['ay']
        # print rx, ry

        # 귀찮아서 에러 없에려고 음수 제거
        Qx1 = sqrt(abs(rx ** 2 - 4.0))
        Qy1 = sqrt(abs(ry ** 2 - 4.0))
        Qx2 = sqrt(0.5 * abs(rx + sqrt(Qx1)))
        Qy2 = sqrt(0.5 * abs(ry + sqrt(Qy1)))

        Mx = Qx2 - 1.0
        My = Qy2 - 1.0

        if self.fitv.dx_status:
            dDx = res['dx'] - self.fitv.dx
            dDpx = res['dpx'] - self.fitv.dpx
            return Mx ** 2 + My ** 2 + dDx ** 2 + dDpx ** 2
        else:
            return Mx ** 2 + My ** 2

    def FitResult(self):
        o = Optic(self.beam[0])
        i = 0
        for ci in self.comlc.clist:
            if ci['tune']:
                ci['cval'] = self.init_arr[i]
                i = i + 1
            ri = self.reflc.ru(ci['ruid'])
            rfo = self.reflc.rform[ri['type']]
            pset = dict(**ri['pset'])
            pset[rfo['cval']['name']] = ci['cval']

            o.add_by_devname(ri['type'], pset)
            if ci['mp']: break

        return o.get_result_u()

    def GetCenterPositions(self):
        self.center_positions = []
        self.start_positions = []
        self.center_svalues = []
        self.start_svalues = []
        r0 = RotationMatrix()
        v0 = ThreeVector(0.0, 0.0, 0.0)
        slength_last = 0.0
        for c in self.GetFirstOptics().devs:
            if c.type.startswith('dipoleEdge'): continue
            self.start_positions.append(v0)

            r_next = c.GetNextRotMatrix(r0)
            v_next = c.GetNextPosition(v0, r0)
            v_cent = c.GetCenterPosition(v0, r0)
            v0 = v_next
            r0 = r_next

            # scalar
            self.start_svalues.append(slength_last)
            self.center_svalues.append(slength_last + c.arg['l'] * .5)
            slength_last = self.center_svalues[-1] + c.arg['l'] * .5

            # vector
            self.center_positions.append(v_cent)

        return self.center_positions

    def OnFit(self, e):
        self.init_arr = []
        for ci in self.comlc.clist:
            if ci['tune']:
                self.init_arr.append(ci['cval'])
            if ci['mp']: break

        self.gfit = GFitter()

    def OnFitTry(self):
        """
    callback by GFitter()
    """
        self.init_arr = []
        for ci in self.comlc.clist:
            if ci['tune']:
                self.init_arr.append(ci['cval'])
            if ci['mp']: break

        # res = minimize(self.FitFunc, self.init_arr) sdfasdf
        res = fmin(self.FitFunc, self.init_arr, ftol=0.0000001)
        # print res

        for i, ri in enumerate(res):
            self.init_arr[i] = ri

    def OnFitApply(self):
        """
    callback by GFitter()

    """
        i = 0
        for ci in self.comlc.clist:
            if ci['tune']:
                ci['cval'] = self.init_arr[i]
                i = i + 1
        self.comlc.ReGenerateList()
        self.Plot()

    def OnNew(self, e):
        self.comlc.clist = []
        self.reflc.rlist = []
        self.comlc.lastuid = 0
        self.reflc.lastuid = 0
        self.reflc.ReGenerateList()
        self.comlc.ReGenerateList()
        # temp save file must be here!

    def OnOpen(self, e):
        dialog = wx.FileDialog(None,  # parent
                               message="Choose an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="Optics Jason Data Format (*.opt)|*.opt|All files (*.*)|*.*",
                               style=wx.FD_OPEN | wx.FD_CHANGE_DIR)

        if dialog.ShowModal() == wx.ID_OK:
            self.filepath = dialog.GetPath()
            self.filename = dialog.GetFilename()
            self.filedir = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # read file
        try:
            f = open(self.filepath)
            data = load(f)
            f.close()
        except:
            print 'data reading error'
            return

        try:
            self.beam = []
            for b in data['beam']:
                ibeam = BeamSpec(b['rest_mass'], b['charge'], b['kinetic_energy'])
                ibeam.set_xabe(b['ax'], b['bx'], b['ex'])
                ibeam.set_yabe(b['ay'], b['by'], b['ey'])
                ibeam.dpop = b['dpop']
                self.beam.append(ibeam)
        except:
            print 'beam loading error'
            return

        try:
            self.comlc.clist = data['clist']
            self.reflc.rlist = data['rlist']
            self.comlc.lastuid = int(data['clastuid'])
            self.reflc.lastuid = int(data['rlastuid'])
        except:
            print 'reference and component list loading error'
            return

        try:
            self.fitv = fitvalue_ab(data['fitv']['ax'],
                                    data['fitv']['bx'],
                                    data['fitv']['ay'],
                                    data['fitv']['by'])

        except:
            print 'target fit value loading error'

        #     try:
        #       self.reflc.ReGenerateList()
        #       self.comlc.ReGenerateList()
        #     except:
        #       print 'refresh table error'
        self.reflc.ReGenerateList()
        self.comlc.ReGenerateList()

        self.Plot()
        self.SetTitle(self.filename)

    def OnAppend(self, e):
        dialog = wx.FileDialog(None,
                               message="Choose an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="Optics Jason Data Format (*.opt)|*.opt|All files (*.*)|*.*",
                               style=wx.FD_OPEN | wx.FD_CHANGE_DIR)

        if dialog.ShowModal() == wx.ID_OK:
            appfilepath = dialog.GetPath()
            appfilename = dialog.GetFilename()
            appfiledir = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # read file
        try:
            f = open(appfilepath)
            data = load(f)
            f.close()
        except:
            print 'data reading error'

        this_rid = self.reflc.lastuid
        this_cid = self.comlc.lastuid
        try:
            # reference uid correction
            ci_modified_index = []
            for ri in data['rlist']:
                rid_org = ri['id']
                ri['id'] = this_rid
                for i, ci in enumerate(data['clist']):
                    if ci_modified_index.count(i):
                        continue
                    if ci['ruid'] == rid_org:
                        ci['ruid'] = this_rid
                        ci_modified_index.append(i)

                this_rid += 1

            # component uid correction
            for ci in data['clist']:
                cid_org = ci['id']
                ci['id'] = this_cid
                for uid in ci['mirroredby']:
                    for cci in data['clist']:
                        if cci['mirror'] == uid:
                            cci['mirror'] = this_cid
                this_cid += 1

        except:
            print 'converting error'

        for ri in data['rlist']:
            self.reflc.rlist.append(ri)
        for ci in data['clist']:
            self.comlc.clist.append(ci)

        self.reflc.lastuid = this_rid
        self.comlc.lastuid = this_cid

        # print self.reflc.rlist
        # print self.comlc.clist
        try:
            self.reflc.ReGenerateList()
            self.comlc.ReGenerateList()
        except:
            print 'refresh table error'
            return

        self.Plot()

    def OnSave(self, event):
        if self.filepath == '':
            self.OnSaveAs(event)
            return
        else:
            self.SaveToFile()

    def OnSaveAs(self, event):
        dialog = wx.FileDialog(None,
                               message="Name an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="Optics Jason Data Format (*.opt)|*.opt|All files (*.*)|*.*",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            self.filepath = dialog.GetPath()
            self.filename = dialog.GetFilename()
            self.filedir = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()
        self.SaveToFile()
        self.SetTitle(self.filename)

    def OnExportT3D(self, event):
        dialog = wx.FileDialog(None,
                               message="Name an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="Trace3D  (*.t3d)|*.t3d",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            # filedir  = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # Trace 3D
        # ==========================================================================
        if filename.endswith('t3d') or filename.endswith('T3D'):
            beam_text = """
 &DATA
 ER={rest_mass:f} Q={charge:f} W={kinetic_energy:f} XI=0.000
 EMITI={emitx:f} {emity:f} 1944.0
 BEAMI={alpx:f} {betx:f} {alpy:f} {bety:f} {alpz:f} {betz:f}
 FREQ=216.82   ICHROM=0  IBS=0  XC=0.0000
      """.format(rest_mass=self.beam[0].rest_mass,
                 charge=self.beam[0].charge,
                 kinetic_energy=self.beam[0].kinetic_energy,
                 emitx=self.beam[0].ex / (mm * mrad) * 4.0,
                 emity=self.beam[0].ey / (mm * mrad) * 4.0,  # 4rms
                 alpx=self.beam[0].ax,
                 betx=self.beam[0].bx / m,
                 alpy=self.beam[0].ay,
                 bety=self.beam[0].by / m,
                 alpz=0.,
                 betz=1.)

            frame_text = """
 XM=40  XPM={xpm:f}  YM=40  DPM={dpm:f}  DWM=100  DPP=180
 XMI=20  XPMI=20  XMF=20  XPMF=20
 DPMI=180 DWMI=500 DPMF=180  DWMF=500
      """.format(xm=1.2 * max(self.optviewer.last_elps[0]['x']),
                 xpm=1.2 * max(self.optviewer.last_elps[0]['xp']),
                 ym=1.2 * max(self.optviewer.last_elps[0]['y']),
                 dpm=100.,  # longitudinal
                 dwm=100.,  # longitudinal
                 dpp=30.,
                 dpmi=30.,
                 dpmf=30.,
                 dwmi=100.,
                 dwmf=100.)

            comp_text = ''
            count = 1
            for ci in self.comlc.clist:
                ri = self.reflc.ru(ci['ruid'])
                if ri['type'] == 'drift' or ri['type'] == 'monitor':

                    if ci['cval'] == 0.0:
                        continue

                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=1 A(1,{index:03d})={length:f}""".format(
                        index=count,
                        name=ri['name'],
                        length=ci['cval'])
                    count += 1

                elif ri['type'] == 'quadrupole':
                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=3 A(1,{index:03d})={fgrad:f} {leng:f} 0 0 0""".format(
                        index=count,
                        name=ri['name'],
                        fgrad=ci['cval'] / (tesla / m),
                        leng=ri['pset']['l'])
                    count += 1

                elif ri['type'] == 'solenoid':
                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=5 A(1,{index:03d})={fgrad:f} {leng:f}""".format(
                        index=count,
                        name=ri['name'],
                        fgrad=ci['cval'] / gauss,
                        leng=ri['pset']['l'])
                    count += 1

                elif ri['type'] == 'rfgap':
                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=10 A(1,{index:03d})={E0TL:f} {phase:f} 0 0 1""".format(
                        index=count,
                        name=ri['name'],
                        E0TL=ci['cval'] / MV,
                        phase=ri['pset']['s'] / deg)
                    count += 1

                elif ri['type'] == 'dipole':
                    fint_k2 = ri['pset']['k2']
                    if fint_k2 == 0.0:
                        fint_k2 = 0.000001

                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=9 A(1,{index:03d})={eb_angle:f} {rho:f} {gap:f} {fint1:f} {fint2:f}""".format(
                        index=count,
                        name=ri['name'],
                        rho=ci['cval'],
                        gap=ri['pset']['g'],
                        eb_angle=ri['pset']['ef'] / degree,
                        fint1=ri['pset']['k1'],
                        fint2=fint_k2)
                    count += 1

                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=8 A(1,{index:03d})={bang:f} {rho:f} {field_index:f} 0""".format(
                        index=count,
                        name=ri['name'],
                        bang=ri['pset']['a'] / degree,
                        rho=ci['cval'],
                        field_index=ri['pset']['n'])
                    count += 1

                    comp_text += """
CMT({index:03d})='{name:s}' NT({index:03d})=9 A(1,{index:03d})={eb_angle:f} {rho:f} {gap:f} {fint1:f} {fint2:f}""".format(
                        index=count,
                        name=ri['name'],
                        rho=ci['cval'],
                        gap=ri['pset']['g'],
                        eb_angle=ri['pset']['eb'] / degree,
                        fint1=ri['pset']['k1'],
                        fint2=fint_k2)
                    count += 1
            comp_text += "\n&END"

            frame_text += 'N1=1 N2={mx:d} SMAX=5.0 PQSMAX=2.0 NEL1=1 NEL2={mx:d} NP1=1 NP2={mx:d}\n'.format(
                mx=count - 1)

            # print comp_text
            f = open(filepath, 'w')
            f.write(beam_text)
            f.write(frame_text)
            f.write(comp_text)
            f.close()

    def OnExportMAD(self, event):
        dialog = wx.FileDialog(None,
                               message="Name an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="MadX (*.mad)|*.mad",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            # filedir  = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # MAD
        # ==========================================================================
        if filename.endswith('mad') or filename.endswith('MAD'):
            beam_text = """
//=============================================================================
// MadX input file, automatically generated by "Optics Expert"
//=============================================================================
//
// original filename : {fname:s}
//
// particle information
// rest mass, m      = {rest_mass:f} (MeV/c2)
// charge state, q   = {charge:f} (eplus) 
// kinetic energy, w = {kinetic_energy:f} (MeV)
// emittance ex, ey  = {emitx:f}, {emity:f} (pi.mm.mrad)
// alpha ax,  ay     = {alpx:f} {alpy:f}
// beta bx, by       = {betx:f} {bety:f}

OPTION, RBARC=FALSE;
    \n""".format(rest_mass=self.beam[0].rest_mass,
                 charge=self.beam[0].charge,
                 kinetic_energy=self.beam[0].kinetic_energy,
                 emitx=self.beam[0].ex / (m * mrad) * 4.0,
                 emity=self.beam[0].ey / (m * mrad) * 4.0,  # 4rms
                 alpx=self.beam[0].ax,
                 betx=self.beam[0].bx / m,
                 alpy=self.beam[0].ay,
                 bety=self.beam[0].by / m,
                 fname=self.filename)

            ref_text = "// class definition\n"
            for i, ri in enumerate(self.reflc.rlist):
                if ri['type'] == 'quadrupole':
                    reftext = "%s: QUADRUPOLE, L:=%.10f;\n"
                    vararr = (ri['name'].replace('-', ''),  # label
                              ri['pset']['l'] / m)
                    tmptxt = reftext % vararr


                elif ri['type'] == 'solenoid':
                    reftext = "%s: SOLENOID, L:=%.10f;\n"
                    vararr = (ri['name'].replace('-', ''),  # label
                              ri['pset']['l'] / m)  # eff_leng
                    tmptxt = reftext % vararr


                elif ri['type'] == 'dipole':
                    K1 = ri['pset']['k1']
                    K2 = ri['pset']['k2']

                    COMTERM = (1. + sin(ri['pset']['ef']) ** 2) / cos(ri['pset']['ef'])
                    I2 = K1 * (1.0 - K1 * K2 * tan(ri['pset']['ef']) / COMTERM)
                    edgefront = "%sE1: DIPEDGE, H:=%.10f, E1:=%.10f, FINT:=%.10f, HGAP=%.10f;\n" % (
                        ri['name'].replace('-', ''),  # label
                        1.0 / ri['cmax'] * m,  # 1/rho
                        ri['pset']['ef'] / rad,  # E1
                        I2,  # FINT
                        0.5 * ri['pset']['g'] / m)  # HGAP
                    COMTERM = (1. + sin(ri['pset']['eb']) ** 2) / cos(ri['pset']['eb'])
                    I2 = K1 * (1.0 - K1 * K2 * tan(ri['pset']['eb']) / COMTERM)
                    edgeback = "%sE2: DIPEDGE, H:=%.10f, E1:=%.10f, FINT:=%.10f, HGAP=%.10f;\n" % (
                        ri['name'].replace('-', ''),  # label
                        1.0 / ri['cmax'] * m,  # 1/rho
                        ri['pset']['eb'] / rad,  # E1
                        I2,  # FINT
                        0.5 * ri['pset']['g'] / m)  # HGAP

                    reftext = edgefront + "%s: SBEND, L:=%.10f, ANGLE:=%.10f;\n" + edgeback
                    vararr = (ri['name'].replace('-', ''),  # label
                              ri['cmax'] * ri['pset']['a'] / m,  # L
                              ri['pset']['a'] / rad)

                    tmptxt = reftext % vararr
                else:
                    tmptxt = ''

                ref_text += tmptxt
            ref_text += '\n'

            # component array=========================================================
            total_length = self.optviewer.last_elps[0]['z'][-1]
            comp_text = '//component array\nBEAMLINE: SEQUENCE, REFER=CENTRE, L=%s;\n' % total_length
            add_opt = ''
            for ii, ci in enumerate(self.comlc.clist):
                ri = self.reflc.ru(ci['ruid'])
                if ri['type'] == 'drift' or ri['type'] == 'monitor': continue

                if ri['type'] == 'quadrupole':
                    add_opt = 'K1:=%.10f, ' % (ci['cval'] / self.beam[0].br * m * m)
                elif ri['type'] == 'dipole':
                    add_opt = 'K0:=%.10f, ' % (1.0 / ci['cval'] * m)
                elif ri['type'] == 'solenoid':
                    add_opt = 'KS:=%.10f, ' % (ci['cval'] / self.beam[0].br / (rad / m))
                else:
                    add_opt = ''

                if ri['type'] == 'dipole':
                    comp_text += '  %sE1: %sE1, AT=%.10f;\n' % (
                        ri['name'].replace('-', '') + "_%d" % ii,
                        ri['name'].replace('-', ''),
                        self.start_svalues[ii] / m)

                comp_text += '  %s: %s, %s AT=%.10f;\n' % (
                    ri['name'].replace('-', '') + "_%d" % ii,
                    ri['name'].replace('-', ''),
                    add_opt,
                    self.center_svalues[ii] / m)

                if ri['type'] == 'dipole':
                    comp_text += '  %sE2: %sE2, AT=%.10f;\n' % (
                        ri['name'].replace('-', '') + "_%d" % ii,
                        ri['name'].replace('-', ''),
                        self.start_svalues[ii + 1] / m)

            comp_text += 'endsequence;\n'

            last_text = """
// Define Beam
// ======================================================================
BEAM, PARTICLE=ion,
      CHARGE={charge:d}, 
      GAMMA={gamm:f}, 
      MASS={rest_mass:f},
      EX={emitx:f}, 
      EY={emity:f}, 
      SEQUENCE=BEAMLINE;
USE, SEQUENCE=BEAMLINE;

SELECT, flag=sectormap, clear;
SELECT, flag=twiss, column=name,s,betx,bety,dx,dpx,x,y;
TWISS, save, betx={betx:f}, bety={bety:f}, alfx={alpx:f}, alfy={alpy:f}, 
       file=optics.dat, sectormap;
PLOT, haxis=s, vaxis=betx,bety, interpolate=true, colour=100;
PLOT, haxis=s, vaxis=dx,dy,dpx, interpolate=true, colour=100;
//PLOT, haxis=s, vaxis=x,y, interpolate=true, colour=100;
stop;
"""
            last_text = last_text.format(
                rest_mass=self.beam[0].rest_mass / GeV,
                charge=int(self.beam[0].charge),
                gamm=gamma(self.beam[0].kinetic_energy_per_amu),
                emitx=self.beam[0].ex / (m * mrad) * 4.0,
                emity=self.beam[0].ey / (m * mrad) * 4.0,  # 4rms
                alpx=self.beam[0].ax,
                betx=self.beam[0].bx / m,
                alpy=self.beam[0].ay,
                bety=self.beam[0].by / m)

            f = open(filepath, 'w')
            f.write(beam_text)
            f.write(ref_text)
            f.write(comp_text)
            f.write(last_text)
            f.close()

    def OnExportAGILE4(self, event):
        dialog = wx.FileDialog(None,
                               message="Name an optics configuration file",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="WinAgile4 (*.lat)|*.lat",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            # filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            # filedir  = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # AGILE
        # ==========================================================================
        if self.filename.endswith('lat') or self.filename.endswith('LAT'):
            header_text = """{date:s}
{file_name:s}
transfer line
{title:s}
403
{n_element:s}
""".format(date=datetime.now().strftime("%Y-%m-%d"),
           file_name=filename[:-3],
           title=filename[:-3],
           n_element=len(self.comlc.clist))

            # component array=========================================================
            component = """{name:s}

{type:s}
{length_m:f}
{horbend_rad:f}
{vertbend_rad:f}
{edge1_rad:f}
{edge2_rad:f}
{k1_m2_inv:f}
{k2_m3_inv:f}
-{aperr_m:f}
{aperr_m:f}
-{aperr_m:f}
{aperr_m:f}
RRRR
{halfgap_m:f}
{fint:f}
0
{rotation_rad:f}
"""
            comp_text = ''
            for ii, ci in enumerate(self.comlc.clist):
                ri = self.reflc.ru(ci['ruid'])
                if ri['type'] == 'drift' or ri['type'] == 'monitor':
                    comp_text += component.format(type='DRIFT',
                                                  length_m=ci['cval'] / m,
                                                  horbend_rad=0.,
                                                  vertbend_rad=0.,
                                                  edge1_rad=0.,
                                                  edge2_rad=0.,
                                                  k1_m2_inv=0.,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)
                elif ri['type'] == 'quadrupole':
                    comp_text += component.format(type='QUADR',
                                                  length_m=ri['pset']['l'] / m,
                                                  horbend_rad=0.,
                                                  vertbend_rad=0.,
                                                  edge1_rad=0.,
                                                  edge2_rad=0.,
                                                  k1_m2_inv=-ci['cval'] / self.beam[0].br * m * m,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)
                elif ri['type'] == 'dipole':
                    comp_text += component.format(type='SBEND',
                                                  length_m=ci['cval'] * ri['pset']['a'] / m,
                                                  horbend_rad=ri['pset']['a'] / rad,
                                                  vertbend_rad=0.,
                                                  edge1_rad=ri['pset']['ef'] / rad,
                                                  edge2_rad=ri['pset']['eb'] / rad,
                                                  k1_m2_inv=0.,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)
                elif ri['type'] == 'solenoid':
                    comp_text += component.format(type='SOL',
                                                  length_m=ri['pset']['l'] / m,
                                                  horbend_rad=0.,
                                                  vertbend_rad=0.,
                                                  edge1_rad=0.,
                                                  edge2_rad=0.,
                                                  k1_m2_inv=0.,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)
                elif ri['type'] == 'marker':
                    comp_text += component.format(type='MARK',
                                                  length_m=ci['cval'] / m,
                                                  horbend_rad=0.,
                                                  vertbend_rad=0.,
                                                  edge1_rad=0.,
                                                  edge2_rad=0.,
                                                  k1_m2_inv=0.,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)
                else:
                    comp_text += component.format(type='MARK',
                                                  length_m=ci['cval'] / m,
                                                  horbend_rad=0.,
                                                  vertbend_rad=0.,
                                                  edge1_rad=0.,
                                                  edge2_rad=0.,
                                                  k1_m2_inv=0.,
                                                  k2_m3_inv=0.,
                                                  aperr_m=ri['pset']['r'] / m,
                                                  halfgap_m=0.,
                                                  rotation_rad=0.)

            last_text = """{start_pos:d}
0
2
 {betx:f}
 {alpx:f}
 {dx:f}
 {ddx:f}
 {bety:f}
 {alpy:f}
 {dy:f}
 {ddy:f}
 {betz:f}
 {alpz:f}
 4.28377856600000E+0001
 3.60568566000000E+0000
 2.75367404400000E+0001
 0.00000000000000E+0000
 0.00000000000000E+0000
 0.00000000000000E+0000
 5.00000000000000E-0001 
 6.00000000000000E-0001
 4.00000000000000E-0001
 7.50000000000000E-0001
 3.00000000000000E-0001
 2.00000000000000E-0001
 1.00000000000000E-0001
 3.00000000000000E-0001
 2.00000000000000E-0001
 {emitx:f}
 {sigx:f}
 {emity:f}
 {sigy:f}
 {emitz:f}
 {sigz:f}
 0.00000000000000E+0000
 1.50000000000000E-0003
0
  12C 6+ 파티클 이름
16
 1.11746700000000E+0001
6
6
 2.00000000000000E-0003
 1.00000000000000E+0000
 7.00000000000000E-0003
User 1
 9.38255000000000E-0001
1
0
1
User 2
 9.38255000000000E-0001
1
0
1
1
0
0
0
0
12
1
0
0
0
0
0
"""
            last_text = last_text.format(
                start_pos=0,
                emitx=self.beam[0].ex / (m * mrad),
                emity=self.beam[0].ey / (m * mrad),
                emitz=1.0,
                alpx=self.beam[0].ax,
                betx=self.beam[0].bx / m,
                alpy=self.beam[0].ay,
                bety=self.beam[0].by / m,
                dx=0.,
                ddx=0.,
                dy=0.,
                ddy=0.,
                betz=20.,
                alpz=0.,
                sigx=4.0,
                sigy=4.0,
                sigz=1.0)

            f = open(self.filepath, 'w')
            f.write(header_text)
            f.write(comp_text)
            f.write(last_text)
            f.close()

    def OnExportCSV(self, e):

        dialog = wx.FileDialog(None,
                               message="",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="Comma Saperated Version (*.csv)|*.csv|All files (*.*)|*.*",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            self.filedir = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        o = self.GetFirstOptics()
        ar = o.get_all_result()
        text_ = "{name:s}, {type:s}, {s:.2f}, {bx:.2f}, {ax:.2f}, {dx:.2f}, {x:.2f}, {by:.2f}, {ay:.2f}, {dy:.2f}, {y:.2f}\n"
        f = copen(filepath, 'w', 'utf-8')
        f.write('name, type, s[m], bx, ax, dx, x*2.24, by, ay, dy, y*2.24\n')
        j = -1
        for opt_inf in ar:
            opt_inf['x'] *= 2.24
            opt_inf['y'] *= 2.24
            if opt_inf['type'] != 'dipoleEdge' and j != -1:
                ri = self.reflc.ru(self.comlc.clist[j]['ruid'])
                opt_inf['name'] = ri['name']
                j += 1
                f.write(text_.format(**opt_inf))
            if j == -1:
                j += 1
        f.close()

    def OnExportTEX(self, event):

        dialog = wx.FileDialog(None,
                               message=" ",
                               defaultDir=getcwd(),
                               defaultFile="",
                               wildcard="LaTeX  (*.tex)|*.tex|All files (*.*)|*.*",
                               style=wx.FD_SAVE)

        if dialog.ShowModal() == wx.ID_OK:
            filepath = dialog.GetPath()
            filename = dialog.GetFilename()
            self.filedir = dialog.GetDirectory()
        else:
            dialog.Destroy()
            return 0
        dialog.Destroy()

        # Latex
        # ==========================================================================
        if not (filename.endswith('tex') or filename.endswith('TEX')):
            return

        title = self.filename[:-4].replace("_", " ")
        tex = """
\\documentclass[a4paper, 11pt, subfigure]{oblivoir}
\\usepackage{subfigure}
\\usepackage{graphicx}
\\usepackage{amssymb}
\\usepackage{amsmath}
\\usepackage{bm}
\\title{%s Optics Summary}
""" % title

        # author skip

        tex += "\date{%s}\n" % datetime.now().strftime("%Y-%m-%d %H:%M")

        tex += """
\\begin{document}

\\maketitle
"""

        # ===========================================================================
        # Opt Table, reference array
        rtable_header = """
\\begin{table}
\\begin{center}
\\caption{Reference table %d}
\\label{tab:ref_table_%d}

\\begin{tabular}{ c | c | c c }
\\hline
\\hline
No. &
Name &
Var. &
Value \\\\
\\hline 
"""
        rtable_footer = """
\\hline
\\end{tabular}

\\end{center}
\\end{table}
"""

        opt_table = rtable_header % (1, 1)
        for i, ri in enumerate(self.reflc.rlist):
            # check whather exceeding row or not ?
            # n_next_row = 0
            # if ri['type'] == 'quadrupole': n_next_row = 4
            # elif ri['type'] == 'solenoid': n_next_row = 4
            # elif ri['type'] == 'dipole': n_next_row = 8
            # elif ri['type'] == 'drift': n_next_row = 1
            # elif ri['type'] == 'monitor': n_next_row = 1

            # print ri['type']
            if ri['type'] == 'quadrupole':
                # get max f grad
                max_fg = 0.0
                for ci in self.comlc.clist:
                    if ri['id'] == ci['ruid']:
                        max_fg = max(max_fg, abs(ci['cval']))

                tmptxt = """
{no:d} & {name:s}(quadrupole) & N of used    & {n_used:d} \\\\
       &          & Max. F.Grad & {max_fgrad:.4f} (T/m) \\\\
       &          & Eff. Length & {l:.2f} (cm) \\\\
       &          & Bore radius & {r:.2f} (cm) \\\\
       &          & Max. F. at Surf. & {max_f:.4f} (T) \\\\
        """.format(no=i + 1,
                   name=ri['name'],
                   n_used=ri['citation'],
                   max_fgrad=max_fg / (tesla / m),
                   l=ri['pset']['l'] / cm,
                   r=ri['pset']['r'] / cm,
                   max_f=max_fg * ri['pset']['r'] / tesla)


            elif ri['type'] == 'solenoid':
                # get max f grad
                max_f = 0.0
                for ci in self.comlc.clist:
                    if ri['id'] == ci['ruid']:
                        max_f = max(max_f, abs(ci['cval']))

                tmptxt = """
{no:d} & {name:s}(solenoid) & N of used   & {n_used:d} \\\\
       &          & Max. field  & {max_f:.4f} (T) \\\\
       &          & Eff. Length & {l:.2f} (cm) \\\\
       &          & Bore radius & {r:.2f} (cm) \\\\
        """.format(no=i + 1,
                   name=ri['name'],
                   n_used=ri['citation'],
                   max_f=max_f / tesla,
                   l=ri['pset']['l'] / cm,
                   r=ri['pset']['r'] / cm)


            elif ri['type'] == 'drift':

                tmptxt = """
{no:d} & {name:s}(drift) & Bore radius & {r:.2f} (cm) \\\\
        """.format(no=i + 1,
                   name=ri['name'],
                   r=ri['pset']['r'] / cm)

            elif ri['type'] == 'monitor':

                tmptxt = """
{no:d} & {name:s}(monitor) & Bore radius & {r:.2f} (cm) \\\\
        """.format(no=i + 1,
                   name=ri['name'],
                   r=ri['pset']['r'] / cm)

            elif ri['type'] == 'dipole':
                tmptxt = """
{no:d} & {name:s}(dipole) & N of used    & {n_used:d} \\\\
       &          & Radius of Curvature & {rho:.4f} (m) \\\\
       &          & Bending angle & {a:.2f} (deg) \\\\
       &          & Edge angle(start) & {eb:.2f} (deg) \\\\
       &          & Edge angle(end) & {ef:.2f} (deg) \\\\
       &          & Field index & {n:.4f} (--) \\\\
       &          & Full pole-gap & {g:.1f} (cm) \\\\
       &          & Good field width & {pw:.1f} (cm) \\\\
        """.format(no=i + 1,
                   name=ri['name'],
                   n_used=ri['citation'],
                   rho=ri['cmax'] / m,
                   a=ri['pset']['a'] / deg,
                   pw=ri['pset']['pw'] / cm,
                   g=ri['pset']['g'] / cm,
                   ef=ri['pset']['ef'] / deg,
                   eb=ri['pset']['eb'] / deg,
                   n=ri['pset']['n'])
            opt_table += tmptxt + '\n\\hline'

        opt_table += rtable_footer

        tex += opt_table

        # ===========================================================================
        # Opt Table, component array
        table_header = """
\\begin{table}
\\begin{center}
\\caption{Beam line optics configuration %d}
\\label{tab:optics_table_%d}

\\begin{tabular}{ c | c | c   c | c }
\\hline
\\hline
No. &
Name &
Type &
Length (mm) &
Position (mm)\\\\
\\hline 
"""
        table_footer = """
\\hline
\\hline
\\end{tabular}

\\end{center}
\\end{table}
"""

        opt_table = ""
        max_row = 30
        table_no = 1
        # cp = self.center_positions
        cp = self.start_positions
        for i, ci in enumerate(self.comlc.clist):
            new_table = i % max_row
            if new_table == 0:
                opt_table += table_header % (table_no, table_no)
                table_no += 1
            ri = self.reflc.ru(ci['ruid'])

            tmp = i + 1
            opt_table += " %d & " % tmp
            opt_table += " %s & " % ri['name']
            opt_table += " %s & " % ri['type']

            if ri['type'] == 'dipole':
                tmp = ri['pset']['a'] * ci['cval']
                opt_table += " %.1f &" % tmp

            elif ri['type'] == 'quadrupole' or ri['type'] == 'solenoid':
                opt_table += " %.1f & " % ri['pset']['l']
            else:
                opt_table += " %.1f & " % ci['cval']

            tmp0 = (cp[i].x, cp[i].y, cp[i].z)
            tmp = u"(%.1f, %.1f, %.1f) " % tmp0
            opt_table += " %s \\\\ \n" % tmp

            # tmp = ci['cval'] / m
            # opt_table += " & & $\\rho$ & %.3f (m)& \\\\ \n" % tmp

            if new_table == max_row - 1:
                opt_table += table_footer

        if new_table != max_row - 1:
            opt_table += table_footer

        tex += opt_table

        texplot = """
\\begin{figure}
\\begin{center}
  \\includegraphics[width=12cm]{%s.png}
  \\caption{%s}
  \\label{fig:%s}
\\end{center}
\\end{figure}
"""

        # tex += texplot % ('betadisp', 'Beta and dispersion function', 'betadisp')
        tex += texplot % ('env', 'Envelope function', 'env')
        # tex += texplot % ('t3d', 'Trace 3D envelope function', 't3d')

        tex += "\\end{document}\n"
        # f.write(tex)

        f = copen(filepath, 'w', 'utf-8')
        f.write(tex)
        # f.write(opt_table)

        f.close()

    def SaveToFile(self):
        f = open(self.filepath, 'wb')
        beams = [b.get_data() for b in self.beam]
        data = {'beam': beams,
                'rlist': self.reflc.rlist,
                'clist': self.comlc.clist,
                'clastuid': self.comlc.lastuid,
                'rlastuid': self.reflc.lastuid,
                'fitv': {'ax': self.fitv.ax,
                         'bx': self.fitv.bx,
                         'ay': self.fitv.ay,
                         'by': self.fitv.by}
                }
        dump(data, f)
        f.close()

    def OnAbout(self, event):
        pass

    def OnParticleChange(self, event):
        """
    """
        m_org = self.beam[0].rest_mass
        q_org = self.beam[0].charge
        w_org = self.beam[0].kinetic_energy

        # rest mass
        dialog = wx.TextEntryDialog(None,
                                    "Rest mass? (MeV/c2)",
                                    "Target rest mass", str(m_org), style=wx.OK | wx.CANCEL)
        if dialog.ShowModal() != wx.ID_OK:
            dialog.Destroy()
            return
        try:
            m1 = abs(float(dialog.GetValue()))
        except:
            dialog.Destroy()
            return
        if m1 == 0.:
            dialog.Destroy()
            return
        dialog.Destroy()

        # charge state
        dialog = wx.TextEntryDialog(None,
                                    "Charge state? (eplus)",
                                    "Target charge state", str(q_org), style=wx.OK | wx.CANCEL)
        if dialog.ShowModal() != wx.ID_OK:
            dialog.Destroy()
            return
        try:
            q1 = abs(float(dialog.GetValue()))
        except:
            dialog.Destroy()
            return
        if q1 <= 0.:
            dialog.Destroy()
            return
        dialog.Destroy()

        # kinetic energy
        dialog = wx.TextEntryDialog(None,
                                    "Kinetic energy? (MeV)",
                                    "Target kinetic energy", str(w_org), style=wx.OK | wx.CANCEL)
        if dialog.ShowModal() != wx.ID_OK:
            dialog.Destroy()
            return
        try:
            w1 = abs(float(dialog.GetValue()))
        except:
            dialog.Destroy()
            return
        if w1 <= 0.:
            dialog.Destroy()
            return
        dialog.Destroy()

        # scale factor
        scale_f = (q_org / q1) * (m1 / m_org) * (betagammaT(w1, m1) / betagammaT(w_org, m_org))

        # reference item limitation check
        has_limit = False
        for ci in self.comlc.clist:
            for ri in self.reflc.rlist:
                if (ci['ruid'] == ri['id']) and \
                        (ri['type'] == 'quadrupole' or ri['type'] == 'solenoid'):
                    # min max check
                    if ci['cval'] * scale_f < ri['cmin'] or \
                            ci['cval'] * scale_f > ri['cmax']:
                        wx.MessageBox('component item %d has an upper/lower limit!' % ci['id'], 'Warning',
                                      wx.OK | wx.ICON_WARNING)
                        has_limit = True

        if has_limit:
            mdlg = wx.MessageDialog(None,
                                    "Reference limitation was detected. Do you want to change min/max value AND perform energy scaling?" % scale_f,
                                    'Ref. item modification confirmation',
                                    wx.YES_NO | wx.ICON_QUESTION)

            if mdlg.ShowModal() != wx.ID_YES:
                mdlg.Destroy()
                return
        for ci in self.comlc.clist:
            for ri in self.reflc.rlist:
                # print ri['type']
                if (ci['ruid'] == ri['id']) and \
                        (ri['type'] == 'quadrupole' or \
                         ri['type'] == 'solenoid'):
                    # min max check
                    # print ri['type']
                    if ci['cval'] * scale_f < ri['cmin'] or ci['cval'] * scale_f > ri['cmax']:
                        ri['cmin'] = ci['cval'] * scale_f
                    ci['cval'] = ci['cval'] * scale_f

        for b in self.beam:
            b.reset(m1, q1, w1)

        # apply to widget
        try:
            self.reflc.ReGenerateList()
            self.comlc.ReGenerateList()
        except:
            print 'refresh table error'

        self.Plot()

    def OnQMScaling(self, event):
        """
    qovera_0 = self.beam[0].qm

    qovera_1 = qovera_0
    qovera_1s = str(qovera_1)
    dialog = wx.TextEntryDialog(None,
                                "Currently %.5f, To what charge to mass ratio? (Q/A)" % qovera_0,
                                "Target Q/A", qovera_1s, style=wx.OK|wx.CANCEL)
    if dialog.ShowModal() != wx.ID_OK:
      dialog.Destroy()
      return
    try:
      qovera_1 = abs(float( dialog.GetValue() ))
    except:
      dialog.Destroy()
      return
    if qovera_1 == 0.:
      dialog.Destroy()
      return

    dialog.Destroy()

    alpha10 = qovera_1 / qovera_0

    mdlg = wx.MessageDialog(None, "Apply scale factor %f to all the magnetic field?" % alpha10,
                           'Scaling confirmation',
                            wx.YES_NO | wx.ICON_QUESTION)

    if mdlg.ShowModal() != wx.ID_YES:
      mdlg.Destroy()
      return
    mdlg.Destroy()

    # reference item limitation check
    has_limit = False
    for ci in self.comlc.clist:
      for ri in self.reflc.rlist:
        if (ci['ruid'] == ri['id']) and \
           ( ri['type'] == 'quadrupole' or \
             ri['type'] == 'solenoid'):
          # min max check
          if ci['cval']/alpha10 < ri['cmin'] or \
             ci['cval']/alpha10 > ri['cmax']:
            wx.MessageBox('component item %d has an upper/lower limit!' % ci['id'], 'Warning', wx.OK | wx.ICON_WARNING )
            has_limit = True

    if has_limit:
      mdlg = wx.MessageDialog(None, "Reference limitation was detected. Do you want to change min/max value AND perform energy scaling?" % alpha10,
                             'Ref. item modification confirmation',
                              wx.YES_NO | wx.ICON_QUESTION)

      if mdlg.ShowModal() != wx.ID_YES:
        mdlg.Destroy()
        return
    for ci in self.comlc.clist:
      for ri in self.reflc.rlist:
        #print ri['type']
        if (ci['ruid'] == ri['id']) and \
           ( ri['type'] == 'quadrupole' or \
             ri['type'] == 'solenoid'):
          # min max check
          #print ri['type']
          if alpha10*ci['cval'] < ri['cmin']:
            ri['cmin'] = ci['cval']/alpha10
          elif alpha10*ci['cval'] > ri['cmax']:
            ri['cmax'] = ci['cval']/alpha10
          ci['cval'] /= alpha10

    for b in self.beam:
      b.reset(qovera_1, b.e)

    # apply to widget
    try:
      self.reflc.ReGenerateList()
      self.comlc.ReGenerateList()
    except:
      print 'refresh table error'

    self.Plot()"""
        pass

    def OnInverse(self, event):
        """
    Inverse optics configurations
    """
        # beam upside-down
        for ii, bi in enumerate(self.beam):
            bi.set_xabe(-self.optviewer.last_elps[ii]['ax'][-1],
                        self.optviewer.last_elps[ii]['bx'][-1] * m,
                        bi.ex)
            bi.set_yabe(-self.optviewer.last_elps[ii]['ay'][-1],
                        self.optviewer.last_elps[ii]['by'][-1] * m,
                        bi.ey)

        # reverse components
        self.comlc.clist.reverse()

        try:
            self.reflc.ReGenerateList()
            self.comlc.ReGenerateList()
        except:
            print 'refresh table error'
            return

        self.Plot()

    def OnQuit(self, event):
        self.Destroy()
        # self.Close(True)

    def OnSetTwissBeam(self, e):
        """
    """
        GTwissBeamSetter(self.beam, updatef=self.Plot)

    def OnSetArbitraryBeam(self, e):
        GBeamDistSetter(self.beam)

    def OnShowOpticsWin(self, e):
        if self.optview_status:
            self.optviewer.OnVisibleSwitch('dummy')
        else:
            self.createPlot()

    def OnShowTrackWin(self, e):
        """
    Finally it must be selected from the MainMenu.
    """
        if not (self.trkviewer is None):
            self.trkviewer.OnQuit('dummy')
        self.trkviewer = GTracksView()
        self.trkviewer.Show()

    def OnShowLayoutWin(self, e):
        pass

    def Plot(self):
        if self.optviewer:
            self.optviewer.DrawNew(self.GetOptList())
            self.UpdateStatusbar('RePlotted')

    def UpdateStatusbar(self, p0):
        """
    p0 : log msg
    p1 : line length
    p2 : last vector calc
    """
        self.statusbar.SetStatusText(p0, 0)
        self.statusbar.SetStatusText("Length(m) : %.6f (m)" %
                                     self.optviewer.last_elps[0]['z'][-1], 1)
        self.statusbar.SetStatusText("Dim.(m) : %s" %
                                     repr(self.optviewer.last_vector), 2)

    def OnSetTargetBeam(self, e):
        pass


class GOpticsView(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, wx.ID_ANY,
                          'Optics View', size=(700, 450))
        self.Bind(wx.EVT_CLOSE, self.OnQuit)

        self.menuData = [
            ("&Canvas", (("&Beta", "Beta function", self.OnBetaView),
                         ("&Gamma", "Gamma Function", self.OnGammaView),
                         ("&Alpha", "Alpha Function", self.OnAlphaView),
                         ("", "", ""),
                         ("&Envelope", "Half Size Function", self.OnEnvelopeView),
                         ("&Dispersion", "Dispersion Function", self.OnDispersionView),
                         ("", "", ""),
                         ("&Beta + Dispersion", "Beta with Dispersion Function", self.OnBetaDispView),
                         ("&Envelope + Dispersion", "Envelope with Dispersion Function", self.OnEnvDispView)
                         )),
            ("&View", (("&Legend", "Legend On/Off", self.OnSwitchLegend),
                       ("", "", "")
                       ))]

        gen = lambda t, yt, x, y: dict(title=t, name=yt, v1=x, v2=y)
        self.p = \
            {'envelope': gen('Envelope Function', 'Half Size [mm]', 'x', 'y'),
             'beta': gen('Beta Function', 'Length [m]', 'bx', 'by'),
             'gamma': gen('Gamma Function', 'Divergence [1/m]', 'gx', 'gy'),
             'alpha': gen('Alpha Function', 'Arb. Unit', 'ax', 'ay'),
             'dispersion': gen('Dispersion Function', '100% dp/p Distance [m]', 'dx', 'dy'),
             'betadisp': gen('Beta Function', 'Length [m]', 'bx', 'by'),
             'envdisp': gen('Envelope Function', 'Half Size [mm]', 'x', 'y')}

        self.last_tx = None

        # create some sizers
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.createMenuBar()

        # canvas
        self.canvas_fig = Figure((5., 4.), dpi=65)
        self.canvas = FigCanvas(self, -1, self.canvas_fig)
        self.canvas_ax = self.canvas_fig.add_subplot(111)
        self.canvas_toolbar = NavigationToolbar(self.canvas)
        self.canvas_grid = False
        self.canvas_legend = True
        self.viewset = 'beta'

        # slider
        self.canvas_slider = wx.Slider(self,
                                       id=-1,
                                       value=0,
                                       minValue=0,
                                       maxValue=nTICKS,
                                       style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        if OSNAME == 'Windows':
            self.canvas_slider.SetTickFreq(nTICKS / 10)
        self.canvas_slider.Bind(wx.EVT_SCROLL, self.OnScrolledCS)
        self.canvas_slider.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.TextUpdate)
        self.UpdateCSBoundary()
        self.canvas_slider.SetValue(0)

        # phase space
        self.phase_fig = Figure((2., 2.), dpi=65)
        self.phase = FigCanvas(self, -1, self.phase_fig)
        self.phase_ax = self.phase_fig.add_subplot(111)
        self.phase_toolbar = NavigationToolbar(self.phase)
        self.phase_z_pos = 0.0
        self.phase_grid = True
        self.phase_legend = True

        # info pad
        self.info = wx.Panel(self, -1)
        self.info_textc = wx.StaticText(self.info, -1, '', (5, 5))

        # layout
        l_vbox = wx.BoxSizer(wx.VERTICAL)
        l_vbox.Add(self.canvas_toolbar)
        l_vbox.Add(self.canvas, 90, wx.LEFT | wx.TOP | wx.GROW)
        l_vbox.Add(self.canvas_slider, 5, wx.EXPAND)
        hbox.Add(l_vbox, 2, wx.EXPAND)

        r_vbox = wx.BoxSizer(wx.VERTICAL)
        r_vbox.Add(self.phase_toolbar)
        r_vbox.Add(self.phase, 50, wx.LEFT | wx.TOP | wx.GROW)
        r_vbox.Add(self.info, 45, wx.EXPAND)
        hbox.Add(r_vbox, 1, wx.EXPAND)

        self.SetSizer(hbox)

        self.updating = True
        self.DrawNew(self.Parent.GetOptList())

    def OnScrolledCS(self, event):
        try:
            value = float(self.canvas_slider.GetValue())
        except:
            value = 0.0
        self.phase_z_pos = self.canvas_slider.minv \
                           + value * fTICKS * \
                           (self.canvas_slider.maxv - self.canvas_slider.minv)
        self.DrawCanvas()
        self.DrawPhase()

    def TextUpdate(self, e):
        bi = self.GenBeamInfo(**self.iv)
        self.info_textc.Destroy()
        self.info_textc = wx.StaticText(self.info, -1, bi, (5, 5))
        self.info_textc.SetFont(wx.Font(8, wx.FONTFAMILY_MODERN, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
        self.info.Refresh()

    def UpdateCSBoundary(self):
        try:
            self.canvas_slider.maxv = max(self.last_elps[0]['z'])
            self.canvas_slider.minv = min(self.last_elps[0]['z'])
        except:
            self.canvas_slider.maxv = 0.0
            self.canvas_slider.minv = 0.0

    def createMenuBar(self):
        menuBar = wx.MenuBar()
        for eachMenuData in self.menuData:
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1]
            menuBar.Append(self.createMenu(menuItems), menuLabel)
        self.SetMenuBar(menuBar)

    def createMenu(self, menuData):
        menu = wx.Menu()
        for eachItem in menuData:
            if len(eachItem) == 2:
                label = eachItem[0]
                subMenu = self.createMenu(eachItem[1])
                menu.AppendMenu(wx.NewId(), label, subMenu)
            else:
                self.createMenuItem(menu, *eachItem)
        return menu

    def createMenuItem(self, menu, label, status, handler, kind=wx.ITEM_NORMAL):
        if not label:
            menu.AppendSeparator()
            return
        menuItem = menu.Append(-1, label, status, kind)
        self.Bind(wx.EVT_MENU, handler, menuItem)

    def OnSwitchGrid(self, e):
        """ """
        self.canvas_grid = not self.canvas_grid
        self.canvas.SetEnableGrid(self.canvas_grid)

    def OnSwitchLegend(self, e):
        """ """
        self.canvas_legend = not self.canvas_legend
        # help(self.canvas)
        # self.canvas.SetEnableLegend(self.canvas_legend)
        self.DrawAllAgain()

    def OnVisibleSwitch(self, e):
        if self.updating:
            self.updating = False
        else:
            self.updating = True
        self.Show(self.updating)

    def OnQuit(self, e):
        self.Parent.optview_status = False
        self.canvas.Destroy()
        self.canvas_toolbar.Destroy()
        self.phase.Destroy()
        self.phase_toolbar.Destroy()
        self.Destroy()

    def OnBetaView(self, e):
        self.viewset = 'beta'
        self.DrawCanvas()

    def OnGammaView(self, e):
        self.viewset = 'gamma'
        self.DrawCanvas()

    def OnAlphaView(self, e):
        self.viewset = 'alpha'
        self.DrawCanvas()

    def OnDispersionView(self, e):
        self.viewset = 'dispersion'
        self.DrawCanvas()

    def OnEnvelopeView(self, e):
        self.viewset = 'envelope'
        self.DrawCanvas()

    def OnBetaDispView(self, e):
        self.viewset = 'betadisp'
        self.DrawCanvas()

    def OnEnvDispView(self, e):
        self.viewset = 'envdisp'
        self.DrawCanvas()

    def DrawEmptyCanvas(self):
        v = self.p[self.viewset]
        title, name = v['title'], v['name']

        ax = self.canvas_ax
        ax.clear()
        ax.grid(True)
        ax.set_title(title)
        ax.set_xlabel('Length [m]')
        ax.set_ylabel(name)

        self.canvas.draw()

    ### @property
    def DrawCanvas(self):
        rp = self.last_elps[0]
        v = self.p[self.viewset]
        title, name, v1, v2 = v['title'], v['name'], v['v1'], v['v2']

        if self.last_tx != None:
            self.last_tx.clear()
            # 뭐 이렇게 해도 지워지지는 않는다. fig 전체를 지우면, 바인딩 되어 있는 그림이며 툴바가 망가지기 때문에,
            # FIXME 버그 정도로 해 두자.
        ax = self.canvas_ax
        ax.clear()
        ax.grid(True)
        ax.set_title(title)
        ax.set_xlabel('Length [m]')
        ax.set_ylabel(name)
        ax.set_xlim(0., self.opt.zleng / m)
        if len(self.Parent.comlc.clist) == 0:
            self.canvas.draw()
            return False

        ax.plot(rp['z'], rp[v1], color='#ff0000')
        ax.plot(rp['z'], rp[v2], color='#0000ff')

        # about legend
        if v1 == 'bx':
            lab1 = "$\\beta_{x}$"
            lab2 = "$\\beta_{y}$"
        elif v1 == 'x':
            lab1 = "$x$"
            lab2 = "$y$"
        elif v1 == 'dx':
            lab1 = "$D_{x}$"
            lab2 = "$D_{y}$"
        elif v1 == 'ax':
            lab1 = "$\\alpha_{x}$"
            lab2 = "$\\alpha_{y}$"
        elif v1 == 'gx':
            lab1 = "$\\gamma_{x}$"
            lab2 = "$\\gamma_{y}$"

        if self.viewset == 'dispersion':
            self.last_tx = None
            dz = rp['z'][1] - rp['z'][0]
            ax.plot(rp['z'], rp['dpx'], '--', color='#00ff00')
            ax.plot(rp['z'], rp['dpy'], '--', color='#ffcccc')
            vmax = max(max(rp[v1]), max(rp[v2]), max(rp['dpx']), max(rp['dpy']))
            vmin = min(min(rp[v1]), min(rp[v2]), min(rp['dpx']), min(rp['dpy']))
            if self.canvas_legend:
                ax.legend((lab1, lab2, "$D'_{x}$", "$D'_{y}$"), loc='best', shadow=True)

        elif self.viewset == 'betadisp' or self.viewset == 'envdisp':
            tx = ax.twinx()
            self.last_tx = tx
            tx.set_ylabel('Dispersion [m]')
            tx.plot(rp['z'], rp['dx'], '--', color='#00ff00')
            tx.plot(rp['z'], rp['dy'], '--', color='#ffcccc')
            vmax = max(max(rp[v1]), max(rp[v2]), max(rp['dx']), max(rp['dy']))
            vmin = min(min(rp[v1]), min(rp[v2]), min(rp['dx']), min(rp['dy']))
            if self.canvas_legend:
                ax.legend((lab1, lab2), loc='left', shadow=True)
                tx.legend(('$D_{x}$', '$D_{y}$'), loc='right', shadow=True)

        else:  # beta or envelope only
            self.last_tx = None
            vmax = max(max(rp[v1]), max(rp[v2]))
            # vmin = min(min(rp[v1]), min(rp[v2]))
            vmin = 0.0
            if self.canvas_legend:
                ax.legend((lab1, lab2), loc='best', shadow=True)
        vv = absmaxabs(vmax, vmin) * 0.05

        # optical components
        ax.set_ylim(vmin - vv, vmax + vv)
        z_pos = 0.0
        for c in self.opt.devs:
            if not c.visible:
                z_pos = z_pos + c.arg['l'] / m
                continue
            if c.type == 'dipole':
                ax.add_patch(rect((z_pos, vmin - vv),
                                  c.arg['l'] / m, vv,
                                  facecolor=c_table['D']))
                ax.add_patch(rect((z_pos, vmax),
                                  c.arg['l'] / m, vv,
                                  facecolor=c_table['D']))
            elif c.type == 'quadrupole':
                if c.arg['G'] > 0.:
                    ax.add_patch(rect((z_pos, vmin - 0.5 * vv),
                                      c.arg['l'] / m, 0.5 * vv,
                                      facecolor=c_table['QF']))
                    ax.add_patch(rect((z_pos, vmax + 0.5 * vv),
                                      c.arg['l'] / m, 0.5 * vv,
                                      facecolor=c_table['QF']))
                else:
                    ax.add_patch(rect((z_pos, vmin - vv),
                                      c.arg['l'] / m, 0.5 * vv,
                                      facecolor=c_table['QD']))
                    ax.add_patch(rect((z_pos, vmax),
                                      c.arg['l'] / m, 0.5 * vv,
                                      facecolor=c_table['QD']))
            elif c.type == 'solenoid':
                ax.add_patch(rect((z_pos, vmin - 0.75 * vv),
                                  c.arg['l'] / m, 0.5 * vv,
                                  facecolor=c_table['S']))
                ax.add_patch(rect((z_pos, vmax + 0.25 * vv),
                                  c.arg['l'] / m, 0.5 * vv,
                                  facecolor=c_table['S']))
            elif c.type == 'monitor':
                ax.add_patch(rect((z_pos, vmin - 0.75 * vv),
                                  c.arg['l'] / m, 0.5 * vv,
                                  facecolor=c_table['M']))
                ax.add_patch(rect((z_pos, vmax + 0.25 * vv),
                                  c.arg['l'] / m, 0.5 * vv,
                                  facecolor=c_table['M']))
            z_pos = z_pos + c.arg['l'] / m
        ax.add_line(mlines.Line2D([self.phase_z_pos, self.phase_z_pos],
                                  [vmin - vv, vmax + vv],
                                  linestyle='--', lw=1.0, color=c_table['curr_pos']))

        # other components
        self.DrawCanvasOthers()

        self.canvas.draw()

    def DrawCanvasOthers(self):
        v = self.p[self.viewset]
        v1, v2 = v['v1'], v['v2']
        ax = self.canvas_ax
        for i, rp in enumerate(self.last_elps[1:]):
            if i == 4: return
            ax.plot(rp['z'], rp[v1], color=CLX[i + 1])
            ax.plot(rp['z'], rp[v2], color=CLY[i + 1])

    def DrawPhase(self):
        data = self.last_elps[0]
        z = self.last_elps[0]['z']
        idx = 1
        if self.phase_z_pos == self.canvas_slider.maxv:
            idx = len(z) - 1
        elif self.phase_z_pos != 0.0:
            for i, zi in enumerate(z):
                if self.phase_z_pos < zi:
                    idx = i
                    break

        if len(self.Parent.comlc.clist) == 0:
            ax, bx, gx = data['ax'][0], data['bx'][0], data['gx'][0]
            ay, by, gy = data['ay'][0], data['by'][0], data['gy'][0]
            dx, dpx, dy, dpy = data['dx'][0], 0.0, data['dy'][0], 0.0
        else:
            dz = z[idx] - z[idx - 1]
            ratio = (self.phase_z_pos - z[idx - 1]) / dz
            iv = {}
            for k in data.keys():
                iv[k] = data[k][idx - 1] + ratio * (data[k][idx] - data[k][idx - 1])
        iv['ex'], iv['ey'] = self.beam.ex / (mm * mrad), self.beam.ey / (mm * mrad)
        iv['xint'], iv['yint'] = sqrt(iv['ex'] / iv['gx']), sqrt(iv['ey'] / iv['gy'])
        iv['xpint'], iv['ypint'] = sqrt(iv['ex'] / iv['bx']), sqrt(iv['ey'] / iv['by'])
        iv['sx'], iv['sy'] = -iv['ax'] / iv['bx'], -iv['ay'] / iv['by']
        self.iv = iv

        axs = self.phase_ax
        axs.clear()
        axs.grid(self.phase_grid)
        axs.set_title("Phase Space")
        axs.set_xlabel('Size [mm]')
        axs.set_ylabel('Divergence [mrad]')
        x, xp = twiss_points(iv['ax'], iv['bx'], iv['ex'], 100)
        y, yp = twiss_points(iv['ay'], iv['by'], iv['ey'], 100)
        xymax = 1.2 * max(iv['x'], iv['y'])
        xypmax = 1.2 * max(iv['xp'], iv['yp'])
        axs.set_xlim(-xymax, xymax)
        axs.set_ylim(-xypmax, xypmax)
        axs.add_line(mlines.Line2D(x, xp, lw=1.0, color=CLX[0]))
        axs.add_line(mlines.Line2D(y, yp, lw=1.0, color=CLY[0]))
        if self.phase_legend:
            axs.legend((r'$x$', r'$y$'), loc='best', shadow=True)

        # other plots!!
        self.DrawPhaseOthers(idx)
        self.phase.draw()

    def DrawPhaseOthers(self, idx):
        if len(self.Parent.comlc.clist) == 0: return

        for i, data in enumerate(self.last_elps[1:]):
            if i == 4: return
            z = data['z']
            dz = z[idx] - z[idx - 1]
            ratio = (self.phase_z_pos - z[idx - 1]) / dz
            ax = data['ax'][idx - 1] + ratio * (data['ax'][idx] - data['ax'][idx - 1])
            bx = data['bx'][idx - 1] + ratio * (data['bx'][idx] - data['bx'][idx - 1])
            ay = data['ay'][idx - 1] + ratio * (data['ay'][idx] - data['ay'][idx - 1])
            by = data['by'][idx - 1] + ratio * (data['by'][idx] - data['by'][idx - 1])
            ex = self.optlist[i + 1].beam.ex / (mm * mrad)
            ey = self.optlist[i + 1].beam.ey / (mm * mrad)
            axs = self.phase_ax
            x, xp = twiss_points(ax, bx, ex, 100)
            y, yp = twiss_points(ay, by, ey, 100)
            axs.add_line(mlines.Line2D(x, xp, lw=1.0, color=CLX[i + 1]))
            axs.add_line(mlines.Line2D(y, yp, lw=1.0, color=CLY[i + 1]))

    def DrawAllAgain(self):
        try:
            self.DrawCanvas()
            self.DrawPhase()
            self.TextUpdate('dummy')
        except:
            print "DrawAllAgain ERROR"

    def DrawNew(self, optlist):
        self.optlist = optlist
        self.last_elps = [o.get_ellipse_U() for o in optlist]
        self.beam = optlist[0].beam
        self.opt = optlist[0]
        # position calculation
        r0 = RotationMatrix()
        v0 = ThreeVector(0.0, 0.0, 0.0)
        for c in self.opt.devs:
            r_next = c.GetNextRotMatrix(r0)
            v_next = c.GetNextPosition(v0, r0)
            # v_cent = c.GetCenterPosition(v0, r0)
            v0 = v_next
            r0 = r_next
        self.last_vector = v0 * 0.001  # 1/m
        # self.canvas_slider.UpdateBoundary()
        self.UpdateCSBoundary()
        self.canvas_slider.SetValue(nTICKS)
        self.phase_z_pos = self.canvas_slider.maxv

        self.DrawAllAgain()

    def GenBeamInfo(self, **kwd):
        s = """  z (m)      : {z:.3f}
  beta(m)    : ({bx:.3f}, {by:.3f})
  alpha      : ({ax:.3f}, {ay:.3f})
  gamma(1/m) : ({gx:.3f}, {gy:.3f})
  slope      : ({sx:.3f}, {sy:.3f})
  D (m)      : ({dx:.5f}, {dy:.3f})
  D'         : ({dpx:.5f}, {dpy:.5f})
  X (mm)     : ({x:.2f}, {y:.2f})
  X' (mrad)  : ({xp:.2f}, {yp:.2f})
  Xi (mm)    : ({xint:.2f}, {yint:.2f})
  X'i (mrad) : ({xpint:.2f}, {ypint:.2f})"""
        return s.format(**kwd)


class GTracksView(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, wx.ID_ANY,
                          'Tracks View', size=(1100, 700))
        self.tracks = []
        self.NINTENSITY = len(ITC)
        self.o = self.Parent.GetOptList()[0]

        self.Bind(wx.EVT_CLOSE, self.OnQuit)

        # canvas
        self.canvas_fig = Figure((5., 4.), dpi=75)
        self.canvas = FigCanvas(self, -1, self.canvas_fig)
        self.canvas_grid = gridspec.GridSpec(3, 1, height_ratios=[2.5, 2.5, 1])
        self.canvas_axx = self.canvas_fig.add_subplot(self.canvas_grid[0])
        self.canvas_axy = self.canvas_fig.add_subplot(self.canvas_grid[1])
        self.canvas_axs = self.canvas_fig.add_subplot(self.canvas_grid[2])
        # self.canvas_axx = self.canvas_fig.add_subplot(311)
        # self.canvas_axy = self.canvas_fig.add_subplot(312)
        # self.canvas_axs = self.canvas_fig.add_subplot(313)
        self.canvas_toolbar = NavigationToolbar(self.canvas)
        self.canvas_fig.subplots_adjust(hspace=0)

        # slider
        self.surv_slider = wx.Slider(self,
                                     id=-1,
                                     value=0,
                                     minValue=0,
                                     maxValue=nTICKS,
                                     style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        if OSNAME == 'Windows':
            self.surv_slider.SetTickFreq(nTICKS / 10)
        self.surv_slider.Bind(wx.EVT_SCROLL, self.OnScrolledCS)
        self.surv_slider.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.DrawPhasePlot)
        self.surv_slider.maxv = self.o.zleng / m
        self.surv_slider.SetValue(0)

        # surv
        # self.survcan_fig = Figure((5., 0.6), dpi=75)
        # self.survcan = FigCanvas(self, -1, self.survcan_fig)
        # self.survcan_ax = self.survcan_fig.add_subplot(111)
        # self.survcan_toolbar = NavigationToolbar(self.survcan)
        self.surv = None

        # phase space
        self.real_fig = Figure((2., 4.), dpi=75)
        self.real = FigCanvas(self, -1, self.real_fig)
        self.real_axx = self.real_fig.add_subplot(211)
        self.real_axy = self.real_fig.add_subplot(212)
        self.real_toolbar = NavigationToolbar(self.real)
        self.real_z_pos = 0.0
        self.real_grid = True
        self.real_legend = True
        self.real_fig.subplots_adjust(hspace=0.41, left=0.18)

        # layout
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        l_vbox = wx.BoxSizer(wx.VERTICAL)
        l_vbox.Add(self.canvas_toolbar)
        # l_vtoolbox = wx.BoxSizer(wx.HORIZONTAL)
        # l_vtoolbox.Add(self.canvas_toolbar)
        # l_vtoolbox.Add(self.survcan_toolbar)
        # l_vbox.Add(l_vtoolbox)

        l_vbox.Add(self.canvas, 70, wx.LEFT | wx.TOP | wx.GROW)
        # l_vbox.Add(self.survcan, 20, wx.EXPAND)
        l_vbox.Add(self.surv_slider, 5, wx.EXPAND)
        hbox.Add(l_vbox, 3, wx.EXPAND)

        r_vbox = wx.BoxSizer(wx.VERTICAL)
        r_vbox.Add(self.real_toolbar)
        r_vbox.Add(self.real, 90, wx.EXPAND)
        hbox.Add(r_vbox, 1, wx.EXPAND)

        self.SetSizer(hbox)

        self.CalcAndDraw()

    def OnScrolledCS(self, event):
        try:
            value = float(self.surv_slider.GetValue())
        except:
            value = 0.0
        self.real_z_pos = value * fTICKS * self.surv_slider.maxv
        self.DrawSurvPlot()

    def InitCanvas(self):
        self.canvas_axx.clear()
        self.canvas_axy.clear()
        # self.canvas_axs.clear()
        self.canvas_axx.grid()
        self.canvas_axy.grid()
        # self.canvas_axs.grid()
        self.canvas_axx.set_ylabel('X [mm]')
        self.canvas_axy.set_ylabel('Y [mm]')
        self.canvas_axx.xaxis.set_visible(False)
        self.canvas_axy.xaxis.set_visible(False)
        # self.canvas_axy.set_xlabel('Length [m]')

        # beam pipe boundary and optics element
        bmax = max(self.o.bmaxx, self.o.bmaxy)
        self.bmax = bmax
        vv = bmax * 0.2
        self.canvas_axx.set_ylim(-1.2 * bmax, 1.2 * bmax)
        self.canvas_axx.set_xlim(0.0, self.o.zleng / m)
        self.canvas_axy.set_ylim(-1.2 * bmax, 1.2 * bmax)
        self.canvas_axy.set_xlim(0.0, self.o.zleng / m)
        # self.canvas_axs.set_ylim(0.0, 1.05)

        z_pos = 0.0
        axx = self.canvas_axx
        axy = self.canvas_axy
        last_x_boundary, last_y_boundary = 0.0, 0.0
        for c in self.o.devs:
            z_next = z_pos + c.arg['l'] / m
            if c.type.startswith("col"):
                if c.type[-1] == 'x':
                    axx.add_line(mlines.Line2D([z_pos, z_next], [c.GetXBoundary(), bmax],
                                               color='#666666', linewidth=1.0))
                    axx.add_line(mlines.Line2D([z_pos, z_next], [-c.GetXBoundary(), -bmax],
                                               color='#666666', linewidth=2.0))
                else:  # y
                    axy.add_line(mlines.Line2D([z_pos, z_next], [c.GetYBoundary(), bmax],
                                               color='#666666', linewidth=2.0))
                    axy.add_line(mlines.Line2D([z_pos, z_next], [-c.GetYBoundary(), -bmax],
                                               color='#666666', linewidth=2.0))

            else:

                if last_x_boundary != c.GetXBoundary():
                    axx.add_line(mlines.Line2D([z_pos, z_pos],
                                               [last_x_boundary, c.GetXBoundary()],
                                               color='#666666'))
                    axx.add_line(mlines.Line2D([z_pos, z_pos],
                                               [-last_x_boundary, -c.GetXBoundary()],
                                               color='#666666'))

                if last_y_boundary != c.GetYBoundary():
                    axy.add_line(mlines.Line2D([z_pos, z_pos],
                                               [last_y_boundary, c.GetYBoundary()],
                                               color='#666666'))
                    axy.add_line(mlines.Line2D([z_pos, z_pos],
                                               [-last_y_boundary, -c.GetYBoundary()],
                                               color='#666666'))

                axx.add_line(mlines.Line2D([z_pos, z_next],
                                           [c.GetXBoundary(), c.GetXBoundary()],
                                           color='#666666'))
                axy.add_line(mlines.Line2D([z_pos, z_next],
                                           [c.GetYBoundary(), c.GetYBoundary()],
                                           color='#666666'))
                axx.add_line(mlines.Line2D([z_pos, z_next],
                                           [-c.GetXBoundary(), -c.GetXBoundary()],
                                           color='#666666'))
                axy.add_line(mlines.Line2D([z_pos, z_next],
                                           [-c.GetYBoundary(), -c.GetYBoundary()],
                                           color='#666666'))

                last_x_boundary, last_y_boundary = c.GetXBoundary(), c.GetYBoundary()

                if c.type == 'dipole':
                    axx.add_patch(rect((z_pos, bmax), c.arg['l'] / m, vv, facecolor=c_table['D']))
                    # axy.add_patch(rect((z_pos, bmax), c.arg['l']/m, vv, facecolor=c_table['D']))
                    # axx.add_patch(rect((z_pos, -bmax), c.arg['l']/m, -vv, facecolor=c_table['D']))
                    axy.add_patch(rect((z_pos, -bmax), c.arg['l'] / m, -vv, facecolor=c_table['D']))
                    # axs.add_patch(rect((z_pos, 1.0), c.arg['l']/m, 0.05, facecolor='#666666'))
                elif c.type == 'quadrupole':
                    rr = 1.1
                    color_foc = c_table['QF']
                    if c.arg['G'] < 0.:
                        rr = 1.0
                        color_foc = c_table['QD']
                    axx.add_patch(rect((z_pos, rr * bmax), c.arg['l'] / m, 0.5 * vv, facecolor=color_foc))
                    # axy.add_patch(rect((z_pos, rr*bmax), c.arg['l']/m, 0.5 * vv, facecolor=color_foc))
                    # axx.add_patch(rect((z_pos, -rr*bmax), c.arg['l']/m, 0.5 * -vv, facecolor=color_foc))
                    axy.add_patch(rect((z_pos, -rr * bmax), c.arg['l'] / m, 0.5 * -vv, facecolor=color_foc))
                    # axs.add_patch(rect((z_pos, 1.0), c.arg['l']/m, 0.05, facecolor='#666666'))
                elif c.type == 'solenoid':
                    axx.add_patch(rect((z_pos, 1.05 * bmax), c.arg['l'] / m, 0.5 * vv, facecolor=c_table['S']))
                    # axy.add_patch(rect((z_pos, 1.05*bmax), c.arg['l']/m, 0.5*vv, facecolor=c_table['S']))
                    # axx.add_patch(rect((z_pos, -1.05*bmax), c.arg['l']/m, -0.5*vv, facecolor=c_table['S']))
                    axy.add_patch(rect((z_pos, -1.05 * bmax), c.arg['l'] / m, -0.5 * vv, facecolor=c_table['S']))
                    # axs.add_patch(rect((z_pos, 1.0), c.arg['l']/m, 0.05, facecolor='#666666'))
                elif c.type == 'monitor':
                    axx.add_patch(rect((z_pos, 1.05 * bmax), c.arg['l'] / m, 0.5 * vv, facecolor=c_table['M']))
                    # axy.add_patch(rect((z_pos, 1.05*bmax), c.arg['l']/m, 0.5 * vv, facecolor=c_table['M']))
                    # axx.add_patch(rect((z_pos, -1.05*bmax), c.arg['l']/m, -0.5*vv, facecolor=c_table['M']))
                    axy.add_patch(rect((z_pos, -1.05 * bmax), c.arg['l'] / m, -0.5 * vv, facecolor=c_table['M']))
            z_pos = z_next

    def OnDispersion(self, e):
        pass

    def OnQuit(self, e):
        self.Parent.trkviewer = None
        self.Destroy()

    def CalcAndDraw(self):
        self.InitCanvas()

        # input beam profile and particle numbers
        beams = ['beam%d' % i for i in xrange(len(self.Parent.beam))]

        if not self.Parent.cutbullet == None:
            beams.append("Cut Dist.")

        dialog = wx.SingleChoiceDialog(self, "Pick A Beam", "Choices", beams)
        if dialog.ShowModal() == wx.ID_OK:
            i = dialog.GetSelection()
            dialog.Destroy()
        else:
            dialog.Destroy()
            self.OnQuit('dummy')

        if i == len(self.Parent.beam):  # last component is user-cut beam
            self.tracks.append(self.Parent.cutbullet)
            self.np = len(self.Parent.cutbullet[1])
            # this is very dirty!!
            self.beam = BeamSpec(self.Parent.beam[0].rest_mass, self.Parent.beam[0].charge,
                                 self.Parent.beam[0].kinetic_energy)
            self.beam.dpop = self.Parent.beam[0].dpop
            self.beam.set_xxp(max(self.Parent.cutbullet[1]), max(self.Parent.cutbullet[2]))
            self.beam.set_yyp(max(self.Parent.cutbullet[3]), max(self.Parent.cutbullet[4]))

        else:  # otherwise gaussian distribution beam that was generated at "beam gen" section
            self.beam = self.Parent.beam[i]
            # input number of particles
            while True:
                dialog = wx.TextEntryDialog(self,
                                            "How many particles do you want to put",
                                            "N of particle : ", "500", style=wx.OK | wx.CANCEL)
                if dialog.ShowModal() == wx.ID_OK:
                    try:
                        self.np = int(dialog.GetValue())
                        dialog.Destroy()
                        if self.np > 0:
                            break
                    except:
                        dialog.Destroy()

                    dialog = wx.MessageDialog(self, "Insert a possitive integer number",
                                              'ERROR',
                                              wx.OK | wx.ICON_ERROR)
                    dialog.ShowModal()
                    dialog.Destroy()

            # particle generations
            x, xp, y, yp, zp = self.beam.gen_gaussian(self.np)
            alive_mask = [self.o.devs[0].IsAlive(x[i], y[i]) for i in xrange(self.np)]
            self.tracks.append((0.0, x, xp, y, yp, zp, alive_mask))

        self.o.get_tracks(self.tracks)  # code exception!!

        # ith-tick = (z_i, x, xp, y, yp, zp, alive_mask)
        dialog = wx.ProgressDialog("Post Processing Status",
                                   "Post-processing...",
                                   3,
                                   style=wx.PD_ELAPSED_TIME | wx.PD_AUTO_HIDE)
        xmin = -self.bmax
        xmax = self.bmax
        ymin = -self.bmax
        ymax = self.bmax
        nbins = 40
        # maxcountx = 0
        # maxcounty = 0

        # print self.tracks

        # histogram binning and filling up!
        # ith-tick = (z_i, x, xp, y, yp, zp, alive_mask)
        zt, hxs, hys = [], [], []
        for tick in self.tracks:
            zt.append(tick[0] / m)

            hx = histo1d(xmin, xmax, nbins)
            for i, xi in enumerate(tick[1]):
                if tick[6][i] == 1:
                    hx.fill(xi)
            # maxcountx = max((maxcountx, hx.get_maxcount()))
            hxs.append(hx)

            hy = histo1d(ymin, ymax, nbins)
            for i, yi in enumerate(tick[3]):
                if tick[6][i] == 1:
                    hy.fill(yi)
            # maxcounty = max((maxcounty, hy.get_maxcount()))
            hys.append(hy)
        dialog.Update(1)

        # normalization and histomaking
        zx = [[] for k in xrange(self.NINTENSITY)]
        xx = [[] for k in xrange(self.NINTENSITY)]
        for i in xrange(len(ITC)):
            for idt, hx in enumerate(hxs):
                for xi in hx.get_bincenter_by_itclevel(i):
                    zx[i].append(zt[idt])
                    xx[i].append(xi)
        dialog.Update(2)

        zy = [[] for k in xrange(self.NINTENSITY)]
        yy = [[] for k in xrange(self.NINTENSITY)]
        for i in xrange(len(ITC)):
            for idt, hy in enumerate(hys):
                for yi in hy.get_bincenter_by_itclevel(i):
                    zy[i].append(zt[idt])
                    yy[i].append(yi)
        dialog.Update(3)

        self.zx = zx
        self.xx = xx
        self.zy = zy
        self.yy = yy

        for i in xrange(self.NINTENSITY):
            self.canvas_axx.scatter(self.zx[i], self.xx[i], c=ITC[i], s=15, linewidths=0)
            self.canvas_axy.scatter(self.zy[i], self.yy[i], c=ITC[i], s=15, linewidths=0)

        dialog.Destroy()
        self.DrawSurvPlot()
        self.canvas.draw()

    def DrawSurvPlot(self):
        if len(self.tracks) < 2: return
        # ith-tick = (z_i, x, xp, y, yp, zp, alive_mask)

        # make surviv ratio plot data
        if self.surv == None:
            surv = 1. / float(self.np)
            zz, ss = [], []
            for tick in self.tracks:
                n_alive = 0
                for i in xrange(self.np):
                    if tick[6][i]: n_alive += 1
                zz.append(tick[0] / m)
                ss.append(float(n_alive) * surv)
            self.surv = (zz, ss)

        self.canvas_axs.clear()
        self.canvas_axs.grid()
        self.canvas_axs.set_ylabel('Surv. Ratio [-]')
        self.canvas_axs.set_xlabel('Length [m]')

        # beam pipe boundary and optics element
        vv = 0.2
        self.canvas_axs.set_ylim(0.0, 1.05)
        self.canvas_axs.set_xlim(0.0, self.o.zleng / m)

        z_pos = 0.0
        axs = self.canvas_axs
        # for c in self.o.devs:
        #  z_next = z_pos + c.arg['l']/m
        #  if c.type == 'dipole':
        #    axs.add_patch(rect((z_pos, -0.2), c.arg['l']/m, vv, facecolor=c_table['D']))
        #  elif c.type == 'quadrupole':
        #    rr = -0.1
        #    color_foc = c_table['QF']
        #    if c.arg['G'] < 0.:
        #      rr = -0.2
        #      color_foc = c_table['QD']
        #    axs.add_patch(rect((z_pos, rr), c.arg['l']/m, 0.5 * vv, facecolor=color_foc))
        #  elif c.type == 'solenoid':
        #    axs.add_patch(rect((z_pos, -0.15), c.arg['l']/m, 0.5 * vv, facecolor=c_table['S']))
        #  z_pos = z_next
        axs.plot(self.surv[0], self.surv[1], c='b')
        axs.add_line(mlines.Line2D([self.real_z_pos, self.real_z_pos],
                                   [-0.2, 1.05],
                                   linestyle='--', color=c_table['curr_pos']))
        self.canvas.draw()

    def DrawPhasePlot(self, e):
        if len(self.tracks) < 2: return

        data = self.Parent.optviewer.last_elps[0]
        z = data['z']
        idx = 1
        if self.real_z_pos == self.surv_slider.maxv:
            idx = len(z) - 1
        elif self.real_z_pos != 0.0:
            for i, zi in enumerate(z):
                if self.real_z_pos < zi:
                    idx = i
                    break

        T = self.tracks[idx]
        i_track = [[], [], [], []]
        for i in xrange(self.np):
            if T[6][i]:
                i_track[0].append(T[1][i])
                i_track[1].append(T[2][i])
                i_track[2].append(T[3][i])
                i_track[3].append(T[4][i])

        # ax, bx = data['ax'][idx], data['bx'][idx]
        # ay, by = data['ay'][idx], data['by'][idx]

        # ex, ey = self.beam.ex / (mm*mrad), self.beam.ey / (mm*mrad)

        emittance_unit = mm * mrad
        xphase = TwissPartTr1D()
        xphase.set_particles(i_track[0], [e * mrad for e in i_track[1]])
        erms_x = xphase.get_rms_emittance()[2] / emittance_unit
        print 'emittance_x rms geom= %.2f' % erms_x

        yphase = TwissPartTr1D()
        yphase.set_particles(i_track[2], [e * mrad for e in i_track[3]])
        erms_y = yphase.get_rms_emittance()[2] / emittance_unit
        print 'emittance_y rms geom= %.2f' % erms_y

        axs = self.real_axx
        axs.clear()
        axs.grid()
        # axs.set_title("Phase Space")
        axs.set_xlabel('Size [mm]')
        axs.set_ylabel('Divergence [mrad]')
        # x, xp = twiss_points(ax, bx, ex, 100)
        # y, yp = twiss_points(ay, by, ey, 100)

        t_xmax = max(i_track[0])
        t_ymax = max(i_track[2])
        xmax = max(t_xmax, t_ymax)
        t_xpmax = max(i_track[1])
        t_ypmax = max(i_track[3])
        xpmax = max(t_xpmax, t_ypmax)
        axs.set_xlim(-1.2 * xmax, 1.2 * xmax)
        axs.set_ylim(-1.2 * xpmax, 1.2 * xpmax)

        # now scatter plots
        axs.scatter(i_track[0], i_track[1], c='r', s=5, alpha=0.5, linewidths=0)
        axs.scatter(i_track[2], i_track[3], c='b', s=5, alpha=0.5, linewidths=0)

        # phase envelope
        # axs.add_line(mlines.Line2D(x, xp, lw=1.0, color=CLX[0]))
        # axs.add_line(mlines.Line2D(y, yp, lw=1.0, color=CLY[0]))

        axs = self.real_axy
        axs.clear()
        axs.grid()
        # axs.set_title("Real Space")
        axs.set_xlabel('X-Size [mm]')
        axs.set_ylabel('Y-Size [mm]')
        axs.set_xlim(-1.2 * xmax, 1.2 * xmax)
        axs.set_ylim(-1.2 * xmax, 1.2 * xmax)
        axs.scatter(i_track[0], i_track[2], c='black', s=5, linewidths=0)
        svr = 100.0 * float(len(i_track[0])) / float(self.np)
        # print svr, len(i_track[1]), self.np
        axs.text(.6 * xmax, .9 * xmax, "%.2f %%" % svr)
        self.real.draw()


class GTwissBeamSetter(wx.Frame):
    def __init__(self, beam, **kwd):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, 1111, 'Beam Set')
        self.beam = beam
        if kwd.has_key('updatef'):
            self.has_apply = True
            self.updatefunc = kwd['updatef']
        else:
            self.has_apply = False

        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, 1111)

        # beam phase canvas
        self.beamcanv_fig = Figure((5., 5.), dpi=60)
        self.beamcanv = FigCanvas(panel, -1, self.beamcanv_fig)
        self.beamcanv_ax = self.beamcanv_fig.add_subplot(111)

        # beam listbox
        self.lbox = wx.ListBox(panel, -1, (20, 20), (80, 240),
                               self.GetBeamList(), wx.LB_SINGLE)
        if OSNAME == 'Windows':
            self.lbox.SetSelection(0)
        self.lbox.Bind(wx.EVT_LISTBOX, self.OnChange)

        # beam informations
        txtSizer = wx.GridBagSizer(hgap=5, vgap=5)
        plist = [["Rest mass (MeV) :", MeV, beam[0].rest_mass, 0],  # 0
                 ["Kinetic energy (MeV) :", MeV, beam[0].kinetic_energy, 0],  # 1
                 ["dp/p (percent) :", perCent, beam[0].dpop, 0],  # 2
                 ["Charge state (eplus) :", 1.0, beam[0].charge, 0],  # 3
                 ["alpx :", 1.0, beam[0].ax, 0],  # 4
                 ["betx (m):", m, beam[0].bx, 0],  # 5
                 ["emitx (pi.mm.mrad, rms):", mm * mrad, beam[0].ex, 0],  # 6
                 ["alpy :", 1.0, beam[0].ay, 0],  # 7
                 ["bety (m):", m, beam[0].by, 0],  # 8
                 ["emity (pi.mm.mrad, rms):", mm * mrad, beam[0].ey, 0]]  # 9
        i = 0
        for p in plist:
            txtSizer.Add(wx.StaticText(panel, -1, p[0]),
                         pos=(i, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
            p[3] = wx.TextCtrl(panel, -1, f2s(p[2] / p[1]))
            # p[3].Bind(wx.EVT_KILL_FOCUS, self.OnChange)
            txtSizer.Add(p[3],
                         pos=(i, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
            i = i + 1
        self.plist = plist

        btnSizer = wx.BoxSizer(wx.VERTICAL)
        bAdd = wx.Button(panel, -1, "Add New")
        bAdd.Bind(wx.EVT_BUTTON, self.OnAdd)
        btnSizer.Add(bAdd)
        bMod = wx.Button(panel, -1, "Modify")
        bMod.Bind(wx.EVT_BUTTON, self.OnModify)
        btnSizer.Add(bMod)
        bDel = wx.Button(panel, -1, "Delete")
        bDel.Bind(wx.EVT_BUTTON, self.OnDelete)
        btnSizer.Add(bDel)
        bApply = wx.Button(panel, -1, "Apply")
        bApply.Bind(wx.EVT_BUTTON, self.OnApply)
        btnSizer.Add(bApply)
        bDone = wx.Button(panel, -1, "Done")
        bDone.Bind(wx.EVT_BUTTON, self.OnDone)
        btnSizer.Add(bDone)

        # layout the widgets
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        mainSizer.Add(self.beamcanv)
        mainSizer.Add(self.lbox)
        mainSizer.Add(txtSizer)
        mainSizer.Add(btnSizer, 1, wx.EXPAND)
        panel.SetSizer(mainSizer)
        mainSizer.Fit(self)

        self.DrawBeamsOnLeft()
        self.Show()

    def GetBeamList(self):
        return ['beam%0d' % i for i in xrange(len(self.beam))]

    def OnChange(self, e):
        i = self.lbox.GetSelection()

        self.plist[0][3].SetValue(str(self.beam[i].rest_mass / self.plist[0][1]))
        self.plist[1][3].SetValue(str(self.beam[i].kinetic_energy / self.plist[1][1]))
        self.plist[2][3].SetValue(str(self.beam[i].dpop / self.plist[2][1]))
        self.plist[3][3].SetValue(str(self.beam[i].charge / self.plist[3][1]))
        self.plist[4][3].SetValue(str(self.beam[i].ax / self.plist[4][1]))
        self.plist[5][3].SetValue(str(self.beam[i].bx / self.plist[5][1]))
        self.plist[6][3].SetValue(str(self.beam[i].ex / self.plist[6][1]))
        self.plist[7][3].SetValue(str(self.beam[i].ay / self.plist[7][1]))
        self.plist[8][3].SetValue(str(self.beam[i].by / self.plist[8][1]))
        self.plist[9][3].SetValue(str(self.beam[i].ey / self.plist[9][1]))
        self.DrawBeamsOnLeft()

    def DrawBeamsOnLeft(self):
        ii = self.lbox.GetSelection()
        self.beamcanv_ax.clear()
        xmin, xpmin = 0.0, 0.0
        for i, b in enumerate(self.beam):
            x, xp = twiss_points(b.ax, b.bx * 0.001, b.ex * 1000.0, 100)  # /1000.0 == /(mm*mrad)
            y, yp = twiss_points(b.ay, b.by * 0.001, b.ey * 1000.0, 100)
            if i == ii:
                wid = 3.0
            else:
                wid = 1.0
            xmin = min(min(xmin, x[0]), y[0])
            xpmin = min(xpmin, min(min(xp), min(yp)))
            self.beamcanv_ax.add_line(mlines.Line2D(x, xp, lw=wid, color=CLX[i]))
            self.beamcanv_ax.add_line(mlines.Line2D(y, yp, lw=wid, color=CLY[i]))
        self.beamcanv_ax.grid(True)
        self.beamcanv_ax.set_xlabel('size [mm]')
        self.beamcanv_ax.set_ylabel('div [mrad]')
        self.beamcanv_ax.set_xlim(1.1 * xmin, -1.1 * xmin)
        self.beamcanv_ax.set_ylim(1.1 * xpmin, -1.1 * xpmin)

        self.beamcanv.draw()

    def OnAdd(self, e):
        b = BeamSpec(float(self.plist[0][3].GetValue()) * self.plist[0][1],
                     float(self.plist[3][3].GetValue()) * self.plist[3][1],
                     float(self.plist[1][3].GetValue()) * self.plist[1][1])
        b.set_dpop(float(self.plist[2][3].GetValue()) * self.plist[2][1])
        b.set_xabe(float(self.plist[4][3].GetValue()) * self.plist[4][1],
                   float(self.plist[5][3].GetValue()) * self.plist[5][1],
                   float(self.plist[6][3].GetValue()) * self.plist[6][1])
        b.set_yabe(float(self.plist[7][3].GetValue()) * self.plist[7][1],
                   float(self.plist[8][3].GetValue()) * self.plist[8][1],
                   float(self.plist[9][3].GetValue()) * self.plist[9][1])
        i = len(self.beam) - 1
        if i == 4:
            d = wx.MessageDialog(self, "Exceed maximum number of beam (5)",
                                 'Twiss Ellipse Adding Error',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return
        i += 1
        self.beam.append(b)
        self.lbox.Append('beam%0d' % i)
        self.lbox.Refresh()
        self.DrawBeamsOnLeft()

    def OnModify(self, e):
        i = self.lbox.GetSelection()
        if i < 0: return
        self.beam[i].reset(float(self.plist[0][3].GetValue()) * self.plist[0][1],
                           float(self.plist[3][3].GetValue()) * self.plist[3][1],
                           float(self.plist[1][3].GetValue()) * self.plist[1][1])
        self.beam[i].set_dpop(float(self.plist[2][3].GetValue()) * self.plist[2][1])
        self.beam[i].set_xabe(float(self.plist[4][3].GetValue()) * self.plist[4][1],
                              float(self.plist[5][3].GetValue()) * self.plist[5][1],
                              float(self.plist[6][3].GetValue()) * self.plist[6][1])
        self.beam[i].set_yabe(float(self.plist[7][3].GetValue()) * self.plist[7][1],
                              float(self.plist[8][3].GetValue()) * self.plist[8][1],
                              float(self.plist[9][3].GetValue()) * self.plist[9][1])
        self.DrawBeamsOnLeft()

    def OnDelete(self, e):
        i = self.lbox.GetSelection()
        if len(self.beam) > 1:
            self.beam.remove(self.beam[i])
            self.lbox.Delete(i)
            self.lbox.Refresh()
            self.DrawBeamsOnLeft()

    def OnDone(self, e):
        self.OnModify('dummy')
        self.Destroy()

    def OnApply(self, e):
        if self.has_apply:
            self.OnModify('dummy')
            self.updatefunc()


class GBeamDistSetter(wx.Frame):
    def __init__(self, beam, **kwd):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, 1111, 'Beam Dist Set')
        self.beam = beam

        self.twiss_apply = False
        self.div_apply = False
        self.col_apply = False
        self.done = False

        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, 1111)

        # beam phase canvas
        self.beamcanv_fig = Figure((12., 5.), dpi=60)
        self.beamcanv = FigCanvas(panel, -1, self.beamcanv_fig)
        self.beamcanv_ax = self.beamcanv_fig.add_subplot(121)
        self.beamcanv_ay = self.beamcanv_fig.add_subplot(122)
        # self.beamcanv_toolbar = NavigationToolbar(self.beamcanv)

        # beam cut selection
        beamCutSizer = wx.GridBagSizer(hgap=5, vgap=5)
        # twiss cut info ___________________________________________________________
        beamCutSizer.Add(wx.StaticText(panel, -1, "Select Twiss :"),
                         pos=(0, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        twissbeams = ["beam%d" % i for i in xrange(len(beam))]
        self.twiss_choice = wx.Choice(panel, -1, (100, 18), choices=twissbeams)
        self.twiss_choice.Bind(wx.EVT_CHOICE, self.OnTwissSelect)
        beamCutSizer.Add(self.twiss_choice,
                         pos=(0, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(wx.StaticLine(panel), pos=(1, 0), span=(1, 2), flag=wx.EXPAND)

        # divergence cut info ______________________________________________________
        line = 2
        self.xpint_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        self.xslop_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        self.ypint_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        self.yslop_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        beamCutSizer.Add(wx.StaticText(panel, -1, "Xp_interior_max :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.xpint_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        beamCutSizer.Add(wx.StaticText(panel, -1, "-Alpx/Betx :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.xslop_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        beamCutSizer.Add(wx.StaticText(panel, -1, "Yp_interior_max :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.ypint_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        beamCutSizer.Add(wx.StaticText(panel, -1, "-Alpy/Bety :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.yslop_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        self.div_btn = wx.Button(panel, -1, "Divergence Cut Apply")
        self.div_btn.Bind(wx.EVT_BUTTON, self.OnDivApply)
        beamCutSizer.Add(self.div_btn, pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1
        beamCutSizer.Add(wx.StaticLine(panel), pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1

        # collimator cut info ______________________________________________________
        self.xcol_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        self.ycol_txtctrl = wx.TextCtrl(panel, -1, "0.0")
        beamCutSizer.Add(wx.StaticText(panel, -1, "X_max :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.xcol_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        beamCutSizer.Add(wx.StaticText(panel, -1, "Y_max :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.ycol_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        self.col_btn = wx.Button(panel, -1, "Collimator Cut Apply")
        self.col_btn.Bind(wx.EVT_BUTTON, self.OnColApply)
        beamCutSizer.Add(self.col_btn, pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1
        beamCutSizer.Add(wx.StaticLine(panel), pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1

        # number of particles ______________________________________________________
        self.np_txtctrl = wx.TextCtrl(panel, -1, "500")
        beamCutSizer.Add(wx.StaticText(panel, -1, "Particle Sample :"),
                         pos=(line, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        beamCutSizer.Add(self.np_txtctrl,
                         pos=(line, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        line += 1
        beamCutSizer.Add(wx.StaticLine(panel), pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1

        # gaussian & uniform fill-up dist button ___________________________________
        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.gaussian_btn = wx.Button(panel, -1, "Gaussian Gen.")
        self.gaussian_btn.Bind(wx.EVT_BUTTON, self.OnGaussDist)
        self.uniform_btn = wx.Button(panel, -1, "Uniform Gen.")
        self.uniform_btn.Bind(wx.EVT_BUTTON, self.OnUniformDist)
        btnSizer.Add(self.gaussian_btn)
        btnSizer.Add(self.uniform_btn)
        beamCutSizer.Add(btnSizer, pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1
        beamCutSizer.Add(wx.StaticLine(panel), pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1
        self.runOpTrack = wx.Button(panel, -1, "Run OpTrack")
        self.runOpTrack.Bind(wx.EVT_BUTTON, self.OnRunOpTrack)
        beamCutSizer.Add(self.runOpTrack, pos=(line, 0), span=(1, 2), flag=wx.EXPAND)
        line += 1

        # layout the widgets _______________________________________________________
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        mainSizer.Add(self.beamcanv)
        mainSizer.Add(beamCutSizer, 1, wx.EXPAND)
        panel.SetSizer(mainSizer)
        mainSizer.Fit(self)

        self.DrawFirst()
        self.beamcanv.draw()
        self.Show()

    def DrawFirst(self):
        """
    draw all twiss beams on canvas as a background
    """
        self.beamcanv_ax.clear()
        self.beamcanv_ay.clear()
        xmin, xpmin = 0.0, 0.0
        for b in self.beam:
            x, xp = twiss_points(b.ax, b.bx * 0.001, b.ex * 1000.0, 100)  # /1000.0 == /(mm*mrad)
            y, yp = twiss_points(b.ay, b.by * 0.001, b.ey * 1000.0, 100)
            xmin = min(min(xmin, x[0]), y[0])
            xpmin = min(xpmin, min(min(xp), min(yp)))
            self.beamcanv_ax.add_line(mlines.Line2D(x, xp, lw=1.0, alpha=0.1, color=CLX[0]))
            self.beamcanv_ay.add_line(mlines.Line2D(y, yp, lw=1.0, alpha=0.1, color=CLY[0]))
        self.beamcanv_ax.grid(True)
        self.beamcanv_ax.set_xlabel('size [mm]')
        self.beamcanv_ax.set_ylabel('div [mrad]')
        self.beamcanv_ay.grid(True)
        self.beamcanv_ay.set_xlabel('size [mm]')
        self.beamcanv_ay.set_ylabel('div [mrad]')
        self.xmax = -1.2 * xmin
        self.xpmax = -1.2 * xpmin

        # update frame boundaries
        if not self.twiss_apply and self.div_apply:
            xpmax = max(-self.xslop * -self.xmax + self.xpint, -self.xslop * self.xmax + self.xpint)
            ypmax = max(-self.yslop * -self.xmax + self.ypint, -self.yslop * self.xmax + self.ypint)
            self.xpmax = max(self.xpmax, max(xpmax, ypmax))

        if self.col_apply:
            self.xmax = max(self.xmax, self.xcol)
            self.xmax = max(self.xmax, self.ycol)

        self.beamcanv_ax.set_xlim(-self.xmax, self.xmax)
        self.beamcanv_ax.set_ylim(-self.xpmax, self.xpmax)
        self.beamcanv_ay.set_xlim(-self.xmax, self.xmax)
        self.beamcanv_ay.set_ylim(-self.xpmax, self.xpmax)

        # plot it
        if self.twiss_apply:
            b = self.beam[self.twiss_selected]
            x, xp = twiss_points(b.ax, b.bx * 0.001, b.ex * 1000.0, 100)  # /1000.0 == /(mm*mrad)
            y, yp = twiss_points(b.ay, b.by * 0.001, b.ey * 1000.0, 100)
            self.beamcanv_ax.add_line(mlines.Line2D(x, xp, lw=2.0, color=CLX[0]))
            self.beamcanv_ay.add_line(mlines.Line2D(y, yp, lw=2.0, color=CLY[0]))

        if self.div_apply:
            xy = [-self.xmax, self.xmax]
            xu = [self.xslop * -self.xmax + self.xpint, self.xslop * self.xmax + self.xpint]
            xl = [self.xslop * -self.xmax - self.xpint, self.xslop * self.xmax - self.xpint]
            yu = [self.yslop * -self.xmax + self.ypint, self.yslop * self.xmax + self.ypint]
            yl = [self.yslop * -self.xmax - self.ypint, self.yslop * self.xmax - self.ypint]
            self.beamcanv_ax.add_line(mlines.Line2D(xy, xu, lw=2.0, color=CLX[1]))
            self.beamcanv_ax.add_line(mlines.Line2D(xy, xl, lw=2.0, color=CLX[1]))
            self.beamcanv_ay.add_line(mlines.Line2D(xy, yu, lw=2.0, color=CLY[1]))
            self.beamcanv_ay.add_line(mlines.Line2D(xy, yl, lw=2.0, color=CLY[1]))

        if self.col_apply:
            xl = [-self.xcol, -self.xcol]
            xr = [self.xcol, self.xcol]
            yl = [-self.ycol, -self.ycol]
            yr = [self.ycol, self.ycol]
            c = [-self.xpmax, self.xpmax]
            self.beamcanv_ax.add_line(mlines.Line2D(xl, c, lw=2.0, color=CLX[1]))
            self.beamcanv_ax.add_line(mlines.Line2D(xr, c, lw=2.0, color=CLX[1]))
            self.beamcanv_ay.add_line(mlines.Line2D(yl, c, lw=2.0, color=CLY[1]))
            self.beamcanv_ay.add_line(mlines.Line2D(yr, c, lw=2.0, color=CLY[1]))

    def OnTwissSelect(self, e):
        """
    variable arrange and enables cuts
    """
        # print "called"
        i = self.twiss_choice.GetCurrentSelection()
        if -1 < i < len(self.beam):
            # print i
            self.twiss_apply = True
            self.twiss_selected = i
            self.DrawFirst()
            self.beamcanv.draw()
        else:
            self.twiss_apply = False

    def OnDivApply(self, e):
        """
    variable arrange and enables cuts
    """
        try:
            xpint = abs(float(self.xpint_txtctrl.GetValue()))
            xslop = float(self.xslop_txtctrl.GetValue())
            ypint = abs(float(self.ypint_txtctrl.GetValue()))
            yslop = float(self.yslop_txtctrl.GetValue())
        except:
            self.div_apply = False
            print 'errrrr'
        self.xpint = xpint
        self.xslop = xslop
        self.ypint = ypint
        self.yslop = yslop
        self.div_apply = True

        self.DrawFirst()
        self.beamcanv.draw()

    def OnColApply(self, e):
        """
    variable arrange and enables cuts
    """
        try:
            xcol = abs(float(self.xcol_txtctrl.GetValue()))
            ycol = abs(float(self.ycol_txtctrl.GetValue()))
        except:
            self.col_apply = False
            print 'errrrr'

        self.xcol = xcol
        self.ycol = ycol
        self.col_apply = True

        self.DrawFirst()
        self.beamcanv.draw()

    def OnGaussDist(self, e):

        np = int(self.np_txtctrl.GetValue())
        if np < 10:
            d = wx.MessageDialog(self, "At least 10 particles required",
                                 'N of particle error',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if self.twiss_apply and not self.col_apply and not self.div_apply:
            d = wx.MessageDialog(self, "Why dont you use normal gaussian beam tracking",
                                 'Only gaussian distribution ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if not self.twiss_apply:
            d = wx.MessageDialog(self, "Try uniform distribution rather this without gaussian selection",
                                 'Only gaussian distribution ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if not self.twiss_apply and not self.col_apply and not self.div_apply:
            d = wx.MessageDialog(self, "Select at least one twiss condition",
                                 'Only gaussian distribution ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        count = 0
        alive_mask = [True for i in xrange(np)]
        bx, bxp, by, byp, bzp = [], [], [], [], []
        if self.div_apply:
            xs = self.xslop * 0.001
            ys = self.yslop * 0.001
            xpint = self.xpint * 0.001
            ypint = self.ypint * 0.001
        while count != np:
            x, xp, y, yp, zp = self.beam[self.twiss_selected].gen_gaussian(1000)
            for i in xrange(1000):
                # collimator test
                if self.col_apply and not (-self.xcol < x[i] < self.xcol and -self.ycol < y[i] < self.ycol):
                    continue
                if self.div_apply:
                    if not (xs * x[i] - xpint < xp[i] < xs * x[i] + xpint and ys * y[i] - ypint < yp[i] < ys * y[
                        i] + ypint):
                        continue
                # print count
                bx.append(x[i])
                bxp.append(xp[i])
                by.append(y[i])
                byp.append(yp[i])
                bzp.append(zp[i])
                count += 1
                if count == np:
                    break

        self.Parent.cutbullet = (0.0, bx, bxp, by, byp, bzp, alive_mask)
        xxp = [1000.0 * xpi for xpi in bxp]
        yyp = [1000.0 * ypi for ypi in byp]

        self.DrawFirst()
        self.beamcanv_ax.scatter(bx, xxp, c='r', s=3, alpha=0.5, linewidth=0)
        self.beamcanv_ay.scatter(by, yyp, c='b', s=3, alpha=0.5, linewidth=0)

        # draw scatter
        self.beamcanv.draw()
        self.done = True

    def OnUniformDist(self, e):

        np = int(self.np_txtctrl.GetValue())
        if np < 10:
            d = wx.MessageDialog(self, "At least 10 particles required",
                                 'N of particle error',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if not self.twiss_apply and self.col_apply and not self.div_apply:
            d = wx.MessageDialog(self, "Only Collimators?",
                                 'boundary define ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if not self.twiss_apply and not self.col_apply and self.div_apply:
            d = wx.MessageDialog(self, "Only Divergence?",
                                 'boundary define ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        if not self.twiss_apply and not self.col_apply and not self.div_apply:
            d = wx.MessageDialog(self, "Select at least one twiss condition",
                                 'Only gaussian distribution ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return

        count = 0
        alive_mask = [True for i in xrange(np)]
        bx, bxp, by, byp, bzp = [], [], [], [], []
        if self.div_apply:
            xs = self.xslop
            ys = self.yslop
            xpint = self.xpint
            ypint = self.ypint
        while count != np:
            xi = uniform(-self.xmax, self.xmax)
            xpi = uniform(-self.xpmax, self.xpmax)
            yi = uniform(-self.xmax, self.xmax)
            ypi = uniform(-self.xpmax, self.xpmax)
            # collimator test
            if self.col_apply and not (-self.xcol < xi < self.xcol and -self.ycol < yi < self.ycol):
                continue
            if self.div_apply:
                if not (xs * xi - xpint < xpi < xs * xi + xpint and ys * yi - ypint < ypi < ys * yi + ypint):
                    continue
            if self.twiss_apply:
                if not self.beam[self.twiss_selected].is_inside(xi, xpi * 0.001, yi, ypi * 0.001):
                    continue
            bx.append(xi)
            bxp.append(xpi * 0.001)
            by.append(yi)
            byp.append(ypi * 0.001)
            bzp.append(uniform(-self.beam[self.twiss_selected].dpop, self.beam[self.twiss_selected].dpop))
            count += 1

        self.Parent.cutbullet = (0.0, bx, bxp, by, byp, bzp, alive_mask)
        xxp = [1000.0 * xpi for xpi in bxp]
        yyp = [1000.0 * ypi for ypi in byp]

        self.DrawFirst()
        self.beamcanv_ax.scatter(bx, xxp, c='r', s=3, alpha=0.5, linewidth=0)
        self.beamcanv_ay.scatter(by, yyp, c='b', s=3, alpha=0.5, linewidth=0)

        # draw scatter
        self.beamcanv.draw()
        self.done = True

    def OnRunOpTrack(self, e):
        if not self.done:
            d = wx.MessageDialog(self, "Generation First!",
                                 'No particle bunch ERROR',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return
        self.Parent.trkviewer = GTracksView()
        self.Parent.trkviewer.Show()


class GReferenceList(wx.ListCtrl):
    pgen = lambda n, ni, u, ui: dict(name=n, ninfo=ni, unit=u, uinfo=ui)
    cols = ('No.', 'Name', 'Type', 'Used')
    rform = {'drift':
                 {'pset': [pgen('r', 'beam pipe inner radius', cm, 'cm')],
                  'cval': {'name': 'l', 'info': 'L length'},
                  'cunit': {'name': 'm', 'value': m}},

             'dipole':
                 {'pset': [pgen('n', 'field index', 1.0, ' '),
                           pgen('a', 'bending angle', deg, 'degree'),
                           pgen('g', 'full gap size', cm, 'cm'),
                           pgen('pw', 'full pole width', cm, 'cm'),
                           pgen('ef', 'front edge angle', deg, 'degree'),
                           pgen('eb', 'back edge angle', deg, 'degree'),
                           pgen('k1', 'FINT(k1)', 1.0, ' '),
                           pgen('k2', 'Higher order FINT(k2)', 1.0, ' ')],
                  'cval': {'name': 'r', 'info': 'R curv radii'},
                  'cunit': {'name': 'm', 'value': m}},

             'quadrupole':
                 {'pset': [pgen('l', 'effective length', cm, 'cm'),
                           pgen('r', 'beam pipe inner radius', cm, 'cm')],
                  'cval': {'name': 'G', 'info': 'G foc. grad.'},
                  'cunit': {'name': 'T/m', 'value': tesla / m}},

             'solenoid':
                 {'pset': [pgen('l', 'effective length', cm, 'cm'),
                           pgen('r', 'beam pipe inner radius', cm, 'cm')],
                  'cval': {'name': 'B', 'info': 'B field'},
                  'cunit': {'name': 'T', 'value': tesla}},

             # 'quadrupoleasym':
             # {'pset' : [pgen('l', 'effective length', cm, 'cm'),
             #            pgen('r', 'beam pipe inner radius', cm, 'cm'),
             #            pgen('a', 'Gx/Gy ratio', 1.0, ' ')],
             #  'cval' : {'name': 'G', 'info': 'Gx foc. grad.'},
             #  'cunit': {'name': 'T/m', 'value': tesla/m}},

             'rfgap':
                 {'pset': [pgen('s', 'synchronous phase', deg, 'deg'),
                           pgen('f', 'frequency', MHz, 'MHz'),
                           pgen('r', 'tube radius', cm, 'cm')],
                  'cval': {'name': 'etl', 'info': 'Eff. gap voltage'},
                  'cunit': {'name': 'MV', 'value': MV}},

             'monitor':
                 {'pset': [pgen('r', 'beam pipe inner radius', cm, 'cm')],
                  'cval': {'name': 'l', 'info': 'L length'},
                  'cunit': {'name': 'cm', 'value': cm}},

             'thinlens':
                 {'pset': [pgen('Fy', 'focal length y', cm, 'cm'),
                           pgen('Fz', 'focal length z', cm, 'cm'),
                           pgen('r', 'beam pipe inner radius', cm, 'cm')],
                  'cval': {'name': 'Fx', 'info': 'focal length x'},
                  'cunit': {'name': 'cm', 'value': cm}},

             #
             # Collimators
             #
             'colx':
                 {'pset': [],
                  'cval': {'name': 'x', 'info': 'x half gap'},
                  'cunit': {'name': 'cm', 'value': cm}},

             'coly':
                 {'pset': [],
                  'cval': {'name': 'y', 'info': 'y half gap'},
                  'cunit': {'name': 'cm', 'value': cm}}

             }
    bmps = {}

    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT, size=(260, 350))
        for col, text in enumerate(self.cols):
            self.InsertColumn(col, text)

        self.lastuid = 0
        self.rlist = [{'id': self.lastuid,
                       'type': 'drift',
                       'name': 'D',
                       'cval': 'l',
                       'cmax': 5.0 * m,
                       'cmin': 0.0,
                       'pset': {'r': 4.0 * cm},
                       'citation': 0}]
        self.lastuid += 1

        for x in self.rform.keys():
            self.bmps[x] = wx.Bitmap(ICONPATH + x + ".png", wx.BITMAP_TYPE_PNG)

        self.GenerateList()

        # bind some interesting events
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)

    def rid_by_name(self, name):
        for i, ri in enumerate(self.rlist):
            if name == ri['name']:
                return i
        return -1

    def ru(self, ruid):
        for ri in self.rlist:
            if ruid == ri['id']:
                return ri
        return -1

    def GenerateList(self):
        self.il = wx.ImageList(16, 16, True)
        for idic in self.rlist:
            self.il.Add(self.bmps[idic['type']])
        self.AssignImageList(self.il, wx.IMAGE_LIST_SMALL)

        for i, item in enumerate(self.rlist):
            index = self.InsertItem(MAXINT, str(i))
            self.SetItem(index, 1, item['name'])
            self.SetItem(index, 2, item['type'])
            self.SetItem(index, 3, str(item['citation']))
            self.SetItemImage(index, i, i)

        # set the width of the columns in various ways
        self.SetColumnWidth(0, 40)
        self.SetColumnWidth(1, 60)
        self.SetColumnWidth(2, 100)
        self.SetColumnWidth(3, wx.LIST_AUTOSIZE_USEHEADER)

    def ReGenerateList(self):
        self.DeleteAllItems()
        self.GenerateList()

    def OnAddItem(self, e):
        d1 = wx.SingleChoiceDialog(None, "Pick a reference item", "Choices", self.rform.keys())
        d1_is_deleted = False
        if d1.ShowModal() == wx.ID_OK:
            itemtype = d1.GetStringSelection()
            # add dummy item
            dpset = {}
            for rfp in self.rform[itemtype]['pset']:
                dpset[rfp['name']] = 0.0

            newitem = {'id': self.lastuid,
                       'type': itemtype,
                       'name': 'NEW_%s' % itemtype,
                       'cval': 0.0,
                       'cmax': 0.0,
                       'cmin': 0.0,
                       'pset': dpset,
                       'citation': 0}
            self.lastuid += 1
            self.rlist.append(newitem)
            i = len(self.rlist) - 1

            self.InsertItem(MAXINT, str(i))
            self.SetItem(i, 1, 'NEW_%s' % itemtype)
            self.SetItem(i, 2, itemtype)
            self.SetItem(i, 3, '0')
            self.SetItemImage(i, i, i)
            self.il.Add(self.bmps[itemtype])

            d1.Destroy()
            d1_is_deleted = True
            d2 = GReferenceItemModifier(self, i)
            d2.ShowModal()
            d2.Destroy()

        if not d1_is_deleted:
            d1.Destroy()

    def OnRemoveItem(self, e):
        selid = self.GetFirstSelected()
        # if not selected
        if selid == -1:
            d = wx.MessageDialog(self, "Select reference first",
                                 'Removing Ref. Component Error',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return
        # if it is still referenced
        if self.rlist[selid]['citation'] != 0:
            d = wx.MessageDialog(self, "It is still referenced",
                                 'Removing Ref. component error',
                                 wx.OK | wx.ICON_ERROR)
            d.ShowModal()
            d.Destroy()
            return
        # remove it
        self.rlist.remove(self.rlist[selid])
        self.DeleteItem(selid)
        self.Refresh()

    def OnModifyItem(self, e):
        rid = self.GetFirstSelected()
        if rid != -1:
            tuner = GReferenceItemModifier(self, rid)
            tuner.ShowModal()
            tuner.Destroy()

    def OnItemActivated(self, event):
        selid = self.GetFirstSelected()
        self.comlc.AppendItem(rid=selid)
        self.Select(selid, True)
        self.Parent.Plot()


class GComponentList(wx.ListCtrl):
    cols = ('No.', 'Name', 'Control', 'Value', 'Unit',
            'X', 'Y', 'Z',
            'S.start', 'S.center')
    clist = []
    scale = 0.1
    lastuid = 0

    def __init__(self, parent, reflc):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT | wx.LC_AUTOARRANGE, size=(380, 350))
        for col, text in enumerate(self.cols):
            self.InsertColumn(col, text)
        self.reflc = reflc
        self.Parent.comlc = self
        self.GenerateList()
        self.AppendItem(rid=0)

        self.rmenu = [['Edit', wx.NewId(), self.OnItemActivated],
                      ['Mirror', wx.NewId(), self.OnItemMirrored],
                      ['Delete', wx.NewId(), self.OnDeleteItem],
                      ['Move Up', wx.NewId(), self.OnPushupItem],
                      ['Move Down', wx.NewId(), self.OnPushdownItem],
                      ['Tune Variable', wx.NewId(), self.OnItemSetTune],
                      ['Information', wx.NewId(), self.OnItemInfo]]

        # selected color
        self.MIRRORING_COLOR = (200, 200, 200, 255)
        self.NORMAL_COLOR = (0, 0, 0, 255)
        self.SELECTED_BGCOLOR = (255, 255, 0, 255)
        self.UNSELECTED_BGCOLOR = (255, 255, 255, 255)
        self.MATCHINGPOINT_BGCOLOR = (255, 100, 100, 255)

        # bind some interesting events
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnItemActivated)
        self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressed)
        self.Bind(wx.EVT_LIST_ITEM_RIGHT_CLICK, self.OnItemRightClick)

    def cid_by_uid(self, uid):
        """
    returns index of list(cid) using uid.
    if such uid is not exist, returns -1

    cid : list index
    uid : id of clist['id']
    """
        for i, ci in enumerate(self.clist):
            if ci['id'] == uid:
                return i
        return -1

    def set_ci_cval(self, cid, value, unit):
        """
    modify ci
    """
        self.clist[cid]['cval'] = value
        self.SetItem(cid, 3, f2s(value / unit))
        for uid in self.clist[cid]['mirroredby']:
            mid = self.cid_by_uid(uid)
            self.clist[mid]['cval'] = value
            self.SetItem(mid, 3, f2s(value / unit))

    def GenerateList(self):

        cp = self.Parent.GetCenterPositions()
        sstart = self.Parent.start_svalues
        scenter = self.Parent.center_svalues

        self.il = wx.ImageList(16, 16, True)
        for ci in self.clist:
            self.il.Add(self.reflc.bmps[self.reflc.ru(ci['ruid'])['type']])
        self.AssignImageList(self.il, wx.IMAGE_LIST_SMALL)

        for i, ci in enumerate(self.clist):
            ri = self.reflc.ru(ci['ruid'])
            rfo = self.reflc.rform[ri['type']]
            index = self.InsertItem(MAXINT, str(i))
            self.SetItem(index, 1, ri['name'])
            self.SetItem(index, 2, rfo['cval']['info'])
            self.SetItem(index, 3, f2s(ci['cval'] / rfo['cunit']['value']))
            self.SetItem(index, 4, rfo['cunit']['name'])
            self.SetItem(index, 5, str(cp[i].x / m))
            self.SetItem(index, 6, str(cp[i].y / m))
            self.SetItem(index, 7, str(cp[i].z / m))
            self.SetItem(index, 8, str(sstart[i] / m))
            self.SetItem(index, 9, str(scenter[i] / m))
            self.SetItemImage(index, i)
        # color setting
        for i, ci in enumerate(self.clist):
            if ci['mirror'] != -1:
                ri = self.reflc.ru(ci['ruid'])
                self.SetItemTextColour(i, self.MIRRORING_COLOR)
                self.SetItem(i, 1, ri['name'] + str(self.cid_by_uid(ci['mirror'])))
            # elif ci['tune']:
            #  self.SetItemBackgroundColour(i, self.SELECTED_BGCOLOR)
            #  self.SetItem(i, 1, ci['refname'] + str(self.cid_by_uid(ci['mirror'])))
            # elif ci['tune']:
            #  self.SetItemBackgroundColour(i, self.selclr)

        # set the width of the columns in various ways
        self.SetColumnWidth(0, 40)
        self.SetColumnWidth(1, 60)
        self.SetColumnWidth(2, 120)
        self.SetColumnWidth(3, 100)
        self.SetColumnWidth(4, 50)
        self.SetColumnWidth(5, 50)
        self.SetColumnWidth(6, 50)
        self.SetColumnWidth(7, 50)
        self.SetColumnWidth(8, 55)
        self.SetColumnWidth(9, 70)

    def AppendItem(self, **kwd):
        """
    cid = clist id : mirroring
        or
    rid = rlist id : normal
    """
        next_cid = len(self.clist)
        if kwd.has_key('cid'):  # cid mirroring
            ci = self.clist[kwd['cid']]
            ri = self.reflc.ru(ci['ruid'])
            ri['citation'] += 1
            rid = self.reflc.rid_by_name(ri['name'])
            self.reflc.SetItem(rid, 3, str(ri['citation']))
            self.clist.append({'id': self.lastuid,
                               'ruid': ci['ruid'],
                               'mirroredby': [],
                               'cval': ci['cval'],
                               'mirror': ci['id'],
                               'tune': False,
                               'mp': False})
            mirroring_uid = self.lastuid
            self.lastuid += 1

            name = "%s_M%d" % (ri['name'], kwd['cid'])

        elif kwd.has_key('rid'):  # rid
            ri = self.reflc.rlist[kwd['rid']]
            ri['citation'] += 1
            self.reflc.SetItem(kwd['rid'], 3, str(ri['citation']))
            self.clist.append({'id': self.lastuid,
                               'ruid': ri['id'],
                               'mirroredby': [],
                               'cval': 0.5 * (ri['cmax'] + ri['cmin']),
                               'mirror': -1,
                               'tune': False,
                               'mp': False})
            self.lastuid += 1
            name = ri['name']
            ci = self.clist[-1]

        # just component list update

        cp = self.Parent.GetCenterPositions()
        sstart = self.Parent.start_svalues
        scenter = self.Parent.center_svalues

        # print len(cp), len(self.clist), len(self.reflc.rlist)
        self.il.Add(self.reflc.bmps[self.reflc.ru(ci['ruid'])['type']])
        rfo = self.reflc.rform[ri['type']]
        self.InsertItem(next_cid, str(next_cid))
        self.SetItem(next_cid, 1, name)
        self.SetItem(next_cid, 2, rfo['cval']['info'])
        self.SetItem(next_cid, 3, f2s(self.clist[-1]['cval'] / rfo['cunit']['value']))
        self.SetItem(next_cid, 4, rfo['cunit']['name'])
        self.SetItem(next_cid, 5, str(cp[-1].x / m))
        self.SetItem(next_cid, 6, str(cp[-1].y / m))
        self.SetItem(next_cid, 7, str(cp[-1].z / m))
        self.SetItem(next_cid, 8, str(sstart[-1] / m))
        self.SetItem(next_cid, 9, str(scenter[-1] / m))

        for i in xrange(next_cid + 1):
            self.SetItemImage(i, i)
        if kwd.has_key('cid'):  # cid mirroring
            self.SetItemTextColour(next_cid, self.MIRRORING_COLOR)
            ci['mirroredby'].append(mirroring_uid)

        self.Refresh()

    def ReGenerateList(self):
        self.DeleteAllItems()

        self.GenerateList()
        for i, ci in enumerate(self.clist):
            if ci['tune']:
                self.SetItemBackgroundColour(i, self.SELECTED_BGCOLOR)
            if ci['mp']:
                self.SetItemBackgroundColour(i, self.MATCHINGPOINT_BGCOLOR)

    def OnItemSetTune(self, e):
        # self.Get
        sel = self.GetFirstSelected()
        # print sel
        while sel != -1:
            if self.clist[sel]['tune']:
                self.clist[sel]['tune'] = False
                self.SetItemBackgroundColour(sel, self.UNSELECTED_BGCOLOR)
            else:
                self.clist[sel]['tune'] = True
                self.SetItemBackgroundColour(sel, self.SELECTED_BGCOLOR)
            sel = self.GetNextSelected(sel)

    def OnItemSetMatchingPoint(self, e):
        # self.Geta
        sel = self.GetFirstSelected()
        # print sel
        if sel != -1:
            if self.clist[sel]['mp']:
                self.clist[sel]['mp'] = False
                self.SetItemBackgroundColour(sel, self.UNSELECTED_BGCOLOR)
            else:
                for i, ci in enumerate(self.clist):
                    if ci['mp']:
                        ci['mp'] = False
                        self.SetItemBackgroundColour(i, self.UNSELECTED_BGCOLOR)
                self.clist[sel]['mp'] = True
                self.SetItemBackgroundColour(sel, self.MATCHINGPOINT_BGCOLOR)

    def OnItemMirrored(self, e):
        sel = self.GetFirstSelected()
        while sel != -1:
            if self.clist[sel]['mirror'] != -1:
                # if it is mirrored already, find its mother's cid
                mirrored_uid = self.clist[sel]['mirror']
                mirrored_cid = self.cid_by_uid(mirrored_uid)
            else:
                mirrored_cid = sel

            self.AppendItem(cid=mirrored_cid)
            sel = self.GetNextSelected(sel)
        self.Parent.Plot()

    def OnItemRightClick(self, e):
        sel = self.GetFocusedItem()
        menu = wx.Menu()
        for mitem in self.rmenu:
            itm = menu.Append(mitem[1], mitem[0])
            self.Bind(wx.EVT_MENU, mitem[2], itm)
        if self.clist[sel]['mirror'] != -1:
            menu.Enable(self.rmenu[0][1], False)
            menu.Enable(self.rmenu[1][1], False)
            menu.Enable(self.rmenu[5][1], False)
        self.PopupMenu(menu, e.GetPoint())
        menu.Destroy()

    def OnItemInfo(self, e):
        cid = self.GetFocusedItem()
        ci = self.clist[cid]
        ri = self.reflc.ru(ci['ruid'])
        typ = ri['type']
        text = " Type : %s \n" % typ

        if typ == 'dipole':
            text += " Name : %s\n" % ri['name']
            b0 = self.Parent.beam[0].br / ci['cval']
            text += " Brho = %f (T.m)\n" % self.Parent.beam[0].br
            text += " B0 = %.5f (T)\n" % (b0 / tesla)
            text += " rho = %.1f (mm)\n" % ci['cval']
            text += " K0 = %.6f (1/m)\n" % (b0 / self.Parent.beam[0].br * m)
            leng = ci['cval'] * ri['pset']['a']
            text += " L_arc = %.1f (mm)\n" % leng
            leng_t = 2.0 * ci['cval'] * sin(0.5 * ri['pset']['a'])
            text += " L_rect = %.1f (mm)\n" % leng_t
            text += " full gap = %.1f (mm)\n" % ri['pset']['g']
            bend = ri['pset']['a'] / deg
            e1 = ri['pset']['ef'] / deg
            e2 = ri['pset']['eb'] / deg
            text += " bending angle = %.3f (deg), %f (rad)\n" % (bend, ri['pset']['a'])
            text += " edge angle (in) = %.3f (deg), %f (rad)\n" % (e1, ri['pset']['ef'])
            # HARD CODED K1=0.45, K2=2.8
            K1 = ri['pset']['k1']
            K2 = ri['pset']['k2']
            COMTERM = (1. + sin(ri['pset']['ef']) ** 2) / cos(ri['pset']['ef'])
            I2_in = K1 * (1.0 - K1 * K2 * tan(ri['pset']['ef']) / COMTERM)
            text += "         I2 (in) = %.6f\n" % I2_in
            text += " edge angle (out)= %.3f (deg), %f (rad)\n" % (e2, ri['pset']['eb'])
            COMTERM = (1. + sin(ri['pset']['eb']) ** 2) / cos(ri['pset']['eb'])
            I2_out = K1 * (1.0 - K1 * K2 * tan(ri['pset']['eb']) / COMTERM)
            text += "        I2 (out) = %.6f\n" % I2_out

            matrix = 1
            otmp = Optic(self.Parent.beam[0])
            psettmp = dict(**ri['pset'])
            psettmp['r'] = ci['cval']
            otmp.add_by_devname(ri['type'], psettmp)

            Mi = gen_identity_matrix()
            for matr in otmp.devs:
                Mi = dot(matr.R(matr.arg['l']), Mi)

            text += "Tr. Matrix : \n"
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[0][0], Mi[0][1] / m, Mi[0][2], Mi[0][3] / m, Mi[0][4], Mi[0][5] / m)
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[1][0] * m, Mi[1][1], Mi[1][2] * m, Mi[1][3], Mi[1][4] * m, Mi[1][5])
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[2][0], Mi[2][1] / m, Mi[2][2], Mi[2][3] / m, Mi[2][4], Mi[2][5] / m)
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[3][0] * m, Mi[3][1], Mi[3][2] * m, Mi[3][3], Mi[3][4] * m, Mi[3][5])
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[4][0], Mi[4][1] / m, Mi[4][2], Mi[4][3] / m, Mi[4][4], Mi[4][5] / m)
            text += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
            Mi[5][0] * m, Mi[5][1], Mi[5][2] * m, Mi[5][3], Mi[5][4] * m, Mi[5][5])

            tmptxt = ''
            for matr in otmp.devs:
                Mi = matr.R(matr.arg['l'])
                tmptxt += "Tr. Matrix : \n"
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
                Mi[0][0], Mi[0][1] / m, Mi[0][2], Mi[0][3] / m, Mi[0][4], Mi[0][5] / m)
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
                Mi[1][0] * m, Mi[1][1], Mi[1][2] * m, Mi[1][3], Mi[1][4] * m, Mi[1][5])
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
                Mi[2][0], Mi[2][1] / m, Mi[2][2], Mi[2][3] / m, Mi[2][4], Mi[2][5] / m)
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
                Mi[3][0] * m, Mi[3][1], Mi[3][2] * m, Mi[3][3], Mi[3][4] * m, Mi[3][5])
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n" % (
                Mi[4][0], Mi[4][1] / m, Mi[4][2], Mi[4][3] / m, Mi[4][4], Mi[4][5] / m)
                tmptxt += "%3.6f  %3.6f  %3.6f  %3.6f  %3.6f  %3.6f\n\n" % (
                Mi[5][0] * m, Mi[5][1], Mi[5][2] * m, Mi[5][3], Mi[5][4] * m, Mi[5][5])
            print tmptxt



        elif typ == 'quadrupole':
            text += " Name : %s\n" % ri['name']
            fgrad = ci['cval'] / (tesla / m)
            text += " G = %.5f (T/m)\n" % fgrad
            leng = ri['pset']['l'] / cm
            text += " Eff. length = %.3f (cm)\n" % leng
            bradius = ri['pset']['r'] / cm
            pfld = ci['cval'] * ri['pset']['r'] / tesla
            text += " bore radius = %.3f (cm)\n" % bradius
            text += " field at poletip = %.3f (T)\n" % pfld
            kvalue = ci['cval'] / self.Parent.beam[0].br * m * m  # / (1/(meter*meter))
            text += " K1 = %.6f (1/m2)\n" % kvalue

        elif typ == 'solenoid':
            text += " Name : %s\n" % ri['name']
            b0 = ci['cval'] / tesla
            text += " B0 = %.5f (T)\n" % b0
            leng = ri['pset']['l'] / cm
            text += " Eff. length = %.3f (cm)\n" % leng
            bradius = ri['pset']['r'] / cm
            text += " bore radius = %.3f (cm)\n" % bradius

        wx.MessageBox(text, 'Information', wx.OK | wx.ICON_INFORMATION)

    def OnItemActivated(self, e):
        cid = self.GetFocusedItem()
        if self.clist[cid]['mirror'] != -1:
            cid = self.cid_by_uid(self.clist[cid]['mirror'])
        tuner = GSlidingTuner(self, cid)
        tuner.ShowModal()
        tuner.Destroy()

    def OnKeyPressed(self, e):
        """
    for key tuning by hand
    """
        selid = self.GetFirstSelected()
        key = e.GetKeyCode()

        if e.ControlDown():
            ci = self.clist[selid]
            if ci['mirror'] != -1:
                ci = self.cid_by_uid(ci['mirror'])
            ri = self.reflc.ru(ci['ruid'])
            ci_unit = self.reflc.rform[ri['type']]['cunit']['value']

            if key == wx.WXK_DOWN:
                cval = ci['cval'] - self.scale * ci_unit
                if cval < ri['cmin']:
                    cval = ri['cmin']
                self.set_ci_cval(selid, cval, ci_unit)
                self.Parent.Plot()

            elif key == wx.WXK_UP:
                cval = ci['cval'] + self.scale * ci_unit
                if cval > ri['cmax']:
                    cval = ri['cmax']
                self.set_ci_cval(selid, cval, ci_unit)
                self.Parent.Plot()

            elif key == wx.WXK_RIGHT:
                if self.scale > 0.0001:
                    self.scale = self.scale * 0.1
                    self.Parent.UpdateStatusbar('scale : %.4f' % self.scale)

            elif key == wx.WXK_LEFT:
                if self.scale < 1.0:
                    self.scale = self.scale * 10.0
                    self.Parent.UpdateStatusbar('scale : %.4f' % self.scale)


        elif key == wx.WXK_DELETE:
            self.OnDeleteItem('dummy')

        elif key == wx.WXK_UP:
            if selid > 0:
                self.Select(selid, False)
                self.Select(selid - 1, True)

        elif key == wx.WXK_DOWN:
            if selid != -1 and selid != len(self.clist) - 1:
                self.Select(selid, False)
                self.Select(selid + 1, True)

        elif key == wx.WXK_SPACE:
            if selid != -1 and self.clist[selid]['mirror'] == -1:
                self.OnItemSetTune('dummy')

        elif key == 77:  # 77 => m
            if selid != -1:
                self.OnItemSetMatchingPoint('dummy')

    def OnAddComponent(self, e):
        rselid = self.reflc.GetFirstSelected()
        if rselid == -1:
            return
        self.AppendItem(rid=rselid)
        self.Parent.Plot()

    def OnDeleteItem(self, e):
        delcid_list = []
        sel1 = self.GetFirstSelected()
        selid = sel1
        while selid != -1:
            if len(self.clist[selid]['mirroredby']) != 0:
                d = wx.MessageDialog(self, "Selected component still mirrored by other component!",
                                     'Delete Error',
                                     wx.OK | wx.ICON_ERROR)
                d.ShowModal()
                d.Destroy()
                return

            delcid_list.append(selid)
            selid = self.GetNextSelected(selid)

        for cid in delcid_list:
            delc = self.clist[cid]
            ri = self.reflc.ru(delc['ruid'])
            if delc['mirror'] != -1:
                # remove this uid from mirroredby list of mirroring component
                mirroring_cid = self.cid_by_uid(delc['mirror'])
                # print self.clist[mirroring_cid]['mirroredby'], delc['id']
                self.clist[mirroring_cid]['mirroredby'].remove(delc['id'])
            ri['citation'] -= 1
            self.clist.remove(delc)
            rid = self.reflc.rid_by_name(ri['name'])
            self.reflc.SetItem(rid, 3, str(ri['citation']))
        self.ReGenerateList()
        self.Parent.Plot()

        if len(self.clist) != 0:
            if sel1 == 0:
                self.Select(0, True)
            else:
                self.Select(sel1 - 1, True)

    def OnPushupItem(self, e):
        selid = self.GetFirstSelected()
        if selid == 0 or selid == -1:
            return
        self.SwapItem(selid)
        self.Select(selid, False)
        self.Select(selid - 1, True)
        # self.Parent.Plot()

    def OnPushdownItem(self, e):
        selid = self.GetFirstSelected()
        if selid == len(self.clist) - 1 or selid == -1:
            return
        self.SwapItem(selid + 1)
        self.Select(selid, False)
        self.Select(selid + 1, True)
        # self.Parent.Plot()

    def SwapItem(self, sid):
        """
    swap sid into sid-1, it means moving it up!!
    """
        self.clist[sid], self.clist[sid - 1] = self.clist[sid - 1], self.clist[sid]

        self.ReGenerateList()
        t = """  
    for cid in (sid-1, sid):
      ci = self.clist[cid]
      ri = self.reflc.ru(ci['ruid'])
      rfo = self.reflc.rform[ri['type']]
      bmp = self.reflc.bmps[ri['type']]
      self.il.Replace(cid, bmp)
      name = ri['name']
      if ci['mirror'] != -1:
        name = '%s_M%d' % (name, self.cid_by_uid(ci['mirror']))
        self.SetItemTextColour(cid, self.MIRRORING_COLOR)
      else:
        self.SetItemTextColour(cid, self.NORMAL_COLOR)
        for uid in ci['mirroredby']:
          mcid = self.cid_by_uid(uid)
          #mci = self.clist[mcid]
          self.SetItem(mcid, 1, '%s_M%d' % (name, cid))

      self.SetItem(cid, 1, name)
      self.SetItem(cid, 2, rfo['cval']['info'])
      self.SetItem(cid, 3, f2s(ci['cval']/rfo['cunit']['value']))
      self.SetItem(cid, 4, rfo['cunit']['name'])
      self.SetItemImage(cid, cid, cid)
      """


class GSlidingTuner(wx.Dialog):
    def __init__(self, comlc, cid):
        self.cid = cid
        ci = comlc.clist[cid]
        ri = comlc.reflc.ru(ci['ruid'])
        self.rfo = comlc.reflc.rform[ri['type']]
        self.unit = self.rfo['cunit']['value']
        infotxt = '{name:s} -> {param:s}'.format(param=self.rfo['cval']['name'],
                                                 name=ri['name'])

        wx.Dialog.__init__(self, comlc, -1, infotxt)

        # psld = GCSlider(self, comlc, cid)
        self.maxv = ri['cmax']
        self.minv = ri['cmin']
        if self.minv == self.maxv:
            value = nTICKS / 2
        else:
            value = int(((ci['cval'] - self.minv) /
                         (self.maxv - self.minv)) * nfTICKS)
        self.sld = wx.Slider(self,
                             id=-1,
                             value=nTICKS / 2,
                             minValue=0,
                             maxValue=nTICKS,
                             pos=(10, 10), size=(400, -1),
                             style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS)
        if OSNAME == 'Windows':
            self.sld.SetTickFreq(nTICKS / 10)
        self.sld.Bind(wx.EVT_SCROLL_ENDSCROLL, self.OnScrolled)
        self.sld.SetValue(value)

        # self.txt = GCTxtCtrl(self, comlc, cid)
        self.txt = wx.TextCtrl(self, -1, f2s(ci['cval'] / self.unit))
        self.txt.Bind(wx.EVT_KILL_FOCUS, self.OnKillFocus)

        # Layout with sizers
        sizer = wx.BoxSizer(wx.VERTICAL)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(hbox1, 0, wx.ALL, 5)
        hbox1.Add(self.txt, 1, wx.ALL)
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(hbox2, 0, wx.ALL, 5)
        unit = self.rfo['cunit']['value']
        hbox2.Add(wx.StaticText(self, -1, f2s(ri['cmin'] / unit)))
        hbox2.Add(self.sld)
        hbox2.Add(wx.StaticText(self, -1, f2s(ri['cmax'] / unit)))
        self.SetSizer(sizer)
        sizer.Fit(self)
        self.sld.SetFocus()

    def OnScrolled(self, event):
        try:
            val = float(self.sld.GetValue())
        except:
            val = 0.0
        value = self.minv + val * fTICKS * (self.maxv - self.minv)

        self.txt.SetValue(f2s(value / self.unit))
        self.txt.Refresh()
        self.Parent.set_ci_cval(self.cid, value, self.unit)
        wx.GetApp().TopWindow.Plot()

    def OnKillFocus(self, event):
        try:
            val = float(self.txt.GetValue())
        except:
            val = 0.0

        nvalue = val * self.unit
        if nvalue < self.minv:
            nvalue = self.minv
            self.txt.SetValue(f2s(nvalue / self.unit))
        elif nvalue > self.maxv:
            nvalue = self.maxv
            self.txt.SetValue(f2s(nvalue / self.unit))

        if self.minv == self.maxv:
            ivue = nTICKS / 2
        else:
            ivue = int(((nvalue - self.minv) /
                        (self.maxv - self.minv)) * float(nTICKS))

        self.sld.SetValue(ivue)
        self.sld.Refresh()
        self.Parent.set_ci_cval(self.cid, nvalue, self.unit)
        wx.GetApp().TopWindow.Plot()


class GReferenceItemModifier(wx.Dialog):
    def __init__(self, parent, rid):
        wx.Dialog.__init__(self, parent, -1, "Reference Item Modifier")

        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, 1111)
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.rid = rid
        ri = self.Parent.rlist[rid]
        rf = self.Parent.rform[ri['type']]
        rfp = rf['pset']
        pnames = list()
        pctrls = list()
        txtSizer = wx.GridBagSizer(hgap=5, vgap=5)

        self.nametxt = wx.StaticText(panel, -1, 'Ref. Name:')
        self.namectl = wx.TextCtrl(panel, -1, ri['name'])
        txtSizer.Add(self.nametxt,
                     pos=(0, 0),
                     flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        txtSizer.Add(self.namectl,
                     pos=(0, 1),
                     flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
        txtSizer.Add(wx.StaticLine(panel),
                     pos=(1, 0), span=(1, 2),
                     flag=wx.EXPAND | wx.TOP | wx.BOTTOM)
        i = 2
        for r in rfp:
            pnames.append(wx.StaticText(panel, -1,
                                        "%s, %s(%s)" % (r['name'], r['ninfo'], r['uinfo'])))

            pctrls.append(wx.TextCtrl(panel, -1, str(ri['pset'][r['name']] / r['unit'])))
            txtSizer.Add(pnames[-1],
                         pos=(i, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
            txtSizer.Add(pctrls[-1],
                         pos=(i, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)
            i = i + 1

        self.pnames = pnames
        self.pctrls = pctrls
        self.ri = ri
        self.rfp = rfp
        txtSizer.Add(wx.StaticLine(panel),
                     pos=(i, 0), span=(1, 2),
                     flag=wx.EXPAND | wx.TOP | wx.BOTTOM)
        i = i + 1
        self.cunit = rf['cunit']['value']
        self.cmaxtxt = wx.StaticText(panel, -1,
                                     "MAX %s, %s(%s)" % (rf['cval']['name'], rf['cval']['info'], rf['cunit']['name']))
        self.cmaxctl = wx.TextCtrl(panel, -1, f2s(ri['cmax'] / self.cunit))
        txtSizer.Add(self.cmaxtxt,
                     pos=(i, 0),
                     flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        txtSizer.Add(self.cmaxctl,
                     pos=(i, 1),
                     flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)

        btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        setBtn = wx.Button(panel, -1, "Set")
        setBtn.Bind(wx.EVT_BUTTON, self.OnButtonPressed)
        btnSizer.Add(setBtn)

        # layout the widgets
        mainSizer.Add(txtSizer)
        mainSizer.Add(btnSizer, 1, wx.EXPAND)
        panel.SetSizer(mainSizer)
        mainSizer.Fit(self)

    def OnButtonPressed(self, e):

        # pset
        tpset = {}
        for i, r in enumerate(self.rfp):
            tpset[r['name']] = float(self.pctrls[i].GetValue()) * r['unit']

        # cmax, cmin define
        cmax = float(self.cmaxctl.GetValue()) * self.cunit
        typ = self.ri['type']
        if typ == 'drift' or \
                typ == 'monitor' or \
                typ == 'rfgap' or \
                typ == 'thinlens':
            if cmax < 0.0: cmax = 0.0
            cmin = 0.0

        if typ == 'dipole':
            cmin = cmax

        elif typ == 'solenoid' or \
                typ.startswith('quadrupole'):
            if cmax < 0.0: cmax = 0.0
            cmin = -cmax

        elif typ.startswith('col'):
            if cmax < 0.0: cmax = -cmax
            cmin = 0.0

        self.ri['cmin'] = cmin
        self.ri['cmax'] = cmax
        self.ri['pset'] = tpset
        self.ri['name'] = self.namectl.GetValue().strip()

        self.Parent.SetItem(self.rid, 1, self.ri['name'])
        self.Refresh()

        mainframe = wx.GetApp().TopWindow
        mainframe.comlc.ReGenerateList()
        mainframe.Plot()

        self.Destroy()


class GFitter(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, wx.GetApp().TopWindow, 1111, 'Fitting')

        self.b = self.Parent.FitResult()
        self.fitv = self.Parent.fitv
        # self.real_last = len(self.Parent.comlc.clist)-1

        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, 1111)

        # beam phase canvas
        self.beamcanv_fig = Figure((5., 5.), dpi=60)
        self.beamcanv = FigCanvas(panel, -1, self.beamcanv_fig)
        self.beamcanv_ax = self.beamcanv_fig.add_subplot(111)

        # beam informations
        txtSizer = wx.GridBagSizer(hgap=5, vgap=5)

        genfparam = lambda ___name, ___unit, ___val, ___inst: \
            {'name': ___name, 'unit': ___unit, 'value': ___val, 'instance': ___inst}
        self.plist = [genfparam("alpx",     1.0, self.fitv.ax,  0),
                      genfparam("betx (m):",  m, self.fitv.bx,  0),
                      genfparam("alpy :",   1.0, self.fitv.ay,  0),
                      genfparam("bety (m):",  m, self.fitv.by,  0),
                      genfparam("Dx  (m):",   m, self.fitv.dx,  0),
                      genfparam("Dx' (m):", 1.0, self.fitv.dpx, 0)]
        self.fitv_org = (self.fitv.ax,
                         self.fitv.bx,
                         self.fitv.ay,
                         self.fitv.by,
                         self.fitv.dx,
                         self.fitv.dpx)

        for i, p in enumerate(self.plist):
            txtSizer.Add(wx.StaticText(panel, -1, p['name']),
                         pos=(i, 0),
                         flag=wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
            p['instance'] = wx.TextCtrl(panel,
                                        value=f2s(p['value'] / p['unit']))


            txtSizer.Add(p['instance'],
                         pos=(i, 1),
                         flag=wx.ALIGN_LEFT | wx.ALIGN_CENTER_VERTICAL)

        btnSizer = wx.BoxSizer(wx.VERTICAL)
        bTryfit = wx.Button(panel, -1, "Try fit")
        bTryfit.Bind(wx.EVT_BUTTON, self.OnTryfit)
        btnSizer.Add(bTryfit)
        bApply = wx.Button(panel, -1, "Apply")
        bApply.Bind(wx.EVT_BUTTON, self.OnApply)
        btnSizer.Add(bApply)
        bDone = wx.Button(panel, -1, "Exit")
        bDone.Bind(wx.EVT_BUTTON, self.OnDone)
        btnSizer.Add(bDone)

        # layout the widgets
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        mainSizer.Add(self.beamcanv)
        # mainSizer.Add(self.lbox)
        mainSizer.Add(txtSizer)
        mainSizer.Add(btnSizer, 1, wx.EXPAND)
        panel.SetSizer(mainSizer)
        mainSizer.Fit(self)

        self.DrawBeamsOnLeft()
        self.Show()

    def DrawBeamsOnLeft(self):
        self.beamcanv_ax.clear()
        xmin, xpmin = 0.0, 0.0

        x, xp = twiss_points(self.b['ax'], self.b['bx'], 1.0, 100)  # /1000.0 == /(mm*mrad)
        y, yp = twiss_points(self.b['ay'], self.b['by'], 1.0, 100)
        fx, fxp = twiss_points(self.fitv.ax, self.fitv.bx, 1.0, 100)  # /1000.0 == /(mm*mrad)
        fy, fyp = twiss_points(self.fitv.ay, self.fitv.by, 1.0, 100)

        xmin = min(min(xmin, x[0]), y[0])
        xmin = min(min(xmin, fx[0]), fy[0])
        xpmin = min(xpmin, min(min(xp), min(yp)))
        xpmin = min(xpmin, min(min(fxp), min(fyp)))

        self.beamcanv_ax.add_line(mlines.Line2D(x, xp, color='#ff0000'))
        self.beamcanv_ax.add_line(mlines.Line2D(y, yp, color='#0000ff'))
        self.beamcanv_ax.add_line(mlines.Line2D(fx, fxp, color='#ff0000', linestyle='--'))
        self.beamcanv_ax.add_line(mlines.Line2D(fy, fyp, color='#0000ff', linestyle='--'))

        self.beamcanv_ax.grid(True)
        self.beamcanv_ax.set_xlabel('size [mm]')
        self.beamcanv_ax.set_ylabel('div [mrad]')
        self.beamcanv_ax.set_xlim(1.1 * xmin, -1.1 * xmin)
        self.beamcanv_ax.set_ylim(1.1 * xpmin, -1.1 * xpmin)

        self.beamcanv.draw()

    def OnDone(self, e):
        self.Destroy()

    def OnTryfit(self, e):
        for ci in self.Parent.comlc.clist:
            if ci['tune']:
                self.Parent.init_arr.append(ci['cval'])
            if ci['mp']: break

        self.fitv.ax = float(self.plist[0]['instance'].GetValue()) * self.plist[0]['unit']
        self.fitv.bx = float(self.plist[1]['instance'].GetValue()) * self.plist[1]['unit']
        self.fitv.ay = float(self.plist[2]['instance'].GetValue()) * self.plist[2]['unit']
        self.fitv.by = float(self.plist[3]['instance'].GetValue()) * self.plist[3]['unit']
        self.fitv.dx = float(self.plist[4]['instance'].GetValue()) * self.plist[4]['unit']
        self.fitv.dpx = float(self.plist[5]['instance'].GetValue()) * self.plist[5]['unit']
        if self.fitv_org[4] != self.fitv.dx or \
           self.fitv_org[5] != self.fitv.dpx:
            self.fitv.dx_status = True

        self.DrawBeamsOnLeft()
        self.Parent.OnFitTry()
        self.b = self.Parent.FitResult()
        self.DrawBeamsOnLeft()

    def OnApply(self, e):
        self.Parent.OnFitApply()


