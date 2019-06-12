#!/usr/bin/env /Users/zschillaci/Software/miniconda3/envs/pyenv/bin/python
import os
import glob
import collections
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

def setPlotStyle():
    plt.style.use('ggplot')
    plt.rcParams['lines.linewidth'] = 2.15
    plt.rcParams['lines.markeredgewidth'] = 0.0
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 18
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.facecolor'] = 'whitesmoke'
    plt.rcParams['grid.color'] = 'black'
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = '0.25'
    plt.rcParams['grid.alpha'] = '0.25'
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['figure.titlesize'] = 'large'
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['patch.edgecolor'] = 'none'
    plt.rcParams.update({'figure.max_open_warning': 0})
setPlotStyle()

colors = ['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e', '#9a0eea', '#af884a', '#8c564b', '#ff474c',
          '#0a481e', '#9467bd', '#be6400', '#25a36f', '#ff5b00', '#dbb40c', '#a2cffe', '#acbf69',
          '#6c3461', '#ffcfdc', '#c20078', '#047495']

### FUNDAMENTAL CONSTANTS ###

# SI Units #

#vacuum permittivity
eps0 = 8.85 * 10**-12

#silicon permittivity factor
eps_si = 11.8

#electric charge
q = 1.60 * 10**-19

# Scaling #
i_scale = 10**6  # uA
c_scale = 10**12 # pF
w_scale = 10**6  # um

#############################

def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_list_of_infiles(path):
    infiles, sensors = [], []

    if os.path.isdir(path):
        allfiles = glob.glob(path + '*.txt')
        for afile in allfiles:
            if (('log' not in afile) and ('note' not in afile) and ('noring' not in afile) and (afile.replace('v2','') not in infiles)):
                index = len(afile) - afile[::-1].find('/')

                infiles.append(afile)
                sensors.append(afile[index : -4])

    return infiles, sensors

def autoscale_y(ax, margin=0.1):
    # Rescale y-axis based on data that is visible given the current xlim of the axis
    # ax -- a matplotlib axes object
    # margin -- the fraction of the total height of the y-data to pad the upper and lower ylims

    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo, hi = ax.get_xlim()
        y_displayed = yd[((xd > lo) & (xd < hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed) - margin*h
        top = np.max(y_displayed) + margin*h
        return bot, top

    lines = ax.get_lines()
    bot, top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if (new_bot < bot): bot = new_bot
        if (new_top > top): top = new_top

    ax.set_ylim(bot, top)

def calculate_width(C, A):
    #C = eps * A / d
    return (0 if (C == 0) else ((eps_si * eps0 * A) / C))

def get_doping_profile(CV, area):
    #N = (-C**3 / q * eps * A**2) / [dC/dV]
    #N = (2 / q * eps * A**2) / [d(C**-2)/dV]

    constant = 1 / (q * eps_si * (10**-2 * eps0) * (10**4 * area)**2) # cm^-3

    width, profile = [], []
    for n in range(len(CV['voltage'])):

        dV, dC = 0, 0
        N = 0
        if ((n > 1) and (n < len(CV['voltage']) - 1)):

            dV = CV['voltage'][n + 1] - CV['voltage'][n - 1] # V

            # Method 1 #
            dC = (1 / c_scale) * (CV['capacitance'][n + 1] - CV['capacitance'][n - 1]) # F
            N = -1 * constant * ((1 / c_scale) * CV['capacitance'][n])**3 * (dV / dC)  # cm^-3

            # Method 2 #
            '''
            dC = (c_scale**2) * (CV['invsquarecapacitance'][n + 1] - CV['invsquarecapacitance'][n - 1]) # F^-2
            N = 2 * constant * (dV / dC)  # cm^-3
            '''

        w = w_scale * calculate_width((1 / c_scale) * CV['capacitance'][n], area) # um

        width.append(w)
        profile.append(N)

    CV['width'] = width
    CV['profile'] = profile

    return CV

class Sensor(object):
    def __init__(self, infile, area=None):
        self.infile = infile
        self.area = area

    def load_IV(self):
        #voltage, pad, ring, current
        data = np.loadtxt(self.infile, delimiter=',', skiprows=4)

        voltage, pad, ring, current = [], [], [], []
        for row in data:
            voltage.append(abs(row[0]))
            pad.append(i_scale * row[1])
            ring.append(i_scale * row[2])
            current.append(i_scale * row[3])

        self.IV = pd.DataFrame({'voltage': voltage, 'pad': pad, 
                                'ring': ring, 'current': current})

        return self.IV

    def load_CV(self):
        #voltage, current, capacitance, dvalue
        data = np.loadtxt(self.infile, delimiter=',', skiprows=1)

        voltage, current, capacitance, dvalue, invsquarecapacitance = [], [], [], [], []
        for n, row in enumerate(data):
            voltage.append(abs(row[0]))
            current.append(row[1])
            capacitance.append(c_scale * row[2])
            dvalue.append(row[3])
            invsquarecapacitance.append(0 if (capacitance[n] == 0) else (1 / capacitance[n]**2))

        self.CV = pd.DataFrame({'voltage': voltage, 'current': current,
                                'capacitance': capacitance, 'dvalue': dvalue,
                                'invsquarecapacitance': invsquarecapacitance})
        self.CV = get_doping_profile(self.CV, self.area)

        return self.CV        

class Wafer(object):
    def __init__(self, path, site, wafer, area=None):
        self.path = path
        self.site = site
        self.wafer = wafer
        self.area = area

        self.inpath_to_IV = self.path + self.wafer + '/IV/'
        self.inpath_to_CV = self.path + self.wafer + '/CV/'

        self.outpath = 'results/' + datetime.today().strftime("%Y-%m-%d")
        mkdir(self.outpath)

        self.outpath_to_IV = self.outpath + '/' + self.site + '/IV/'
        self.outpath_to_CV = self.outpath + '/' + self.site + '/CV/'

        self.load_IVs()
        self.load_CVs()

    def load_IVs(self):
        infiles, sensors = get_list_of_infiles(self.inpath_to_IV)

        self.IVs = collections.OrderedDict()
        for infile, sensor in zip(infiles, sensors):
            self.IVs[sensor] = Sensor(infile).load_IV()

    def load_CVs(self):
        infiles, sensors = get_list_of_infiles(self.inpath_to_CV)

        self.CVs = collections.OrderedDict()
        for infile, sensor in zip(infiles, sensors):
            self.CVs[sensor] = Sensor(infile, self.area).load_CV()

    def plot_IVs(self, ylog=True):
        mkdir(self.outpath_to_IV)

        ### pad ###
        fig = plt.figure('pad current vs. voltage' , figsize=(12,8))
        ax = fig.add_subplot(111)
        for n, (sensor, IV) in enumerate(self.IVs.items()):
            plt.plot(IV['voltage'], IV['pad'], 'o--', label=sensor, color=colors[n])

        if ylog:
            plt.yscale('log')
        autoscale_y(ax, margin=5)
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('I$_{Pad}$ [uA]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_IV + '/' + self.wafer + '-padIV.pdf')
        plt.clf()
        plt.close()

        ### ring ###
        fig = plt.figure('ring current vs. voltage' , figsize=(12,8))
        ax = fig.add_subplot(111)
        for n, (sensor, IV) in enumerate(self.IVs.items()):
            plt.plot(IV['voltage'], IV['ring'], 'o--', label=sensor, color=colors[n])

        if ylog:
            plt.yscale('log')
        autoscale_y(ax, margin=5)
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('I$_{Ring}$ [uA]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_IV + '/' + self.wafer + '-ringIV.pdf')
        plt.clf()
        plt.close()

    def plot_CVs(self):
        mkdir(self.outpath_to_CV)

        ### capacitance ###
        fig = plt.figure('cv analysis', figsize=(12, 8))
        ax = fig.add_subplot(111)

        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['voltage'], CV['capacitance'], 'o--', label=sensor, color=colors[n])

        autoscale_y(ax, margin=0.1)
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('C [pF]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-CV.pdf')
        plt.clf()
        plt.close()

        ### invsquarecapacitance ###
        fig = plt.figure('cv analysis', figsize=(12, 8))
        ax = fig.add_subplot(111)

        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['voltage'], CV['invsquarecapacitance'], 'o--', label=sensor, color=colors[n])

        autoscale_y(ax, margin=0.1)
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('C$^{-2}$ [pF]$^{-2}$')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-invsquareCV.pdf')
        plt.clf()
        plt.close()

        ### profile ###
        fig = plt.figure('cv analysis', figsize=(12, 8))
        ax = fig.add_subplot(111)

        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['width'], CV['profile'], 'o--', label=sensor, color=colors[n])

        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        ax.set_xscale('log')
        ax.set_xlim(0.1, 100.0)
        ax.set_xticks([0.1, 1.0, 10.0, 100.0])
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.tick_params(axis='x', reset=False, which='major', length=6, width=2)
        ax.tick_params(axis='x', reset=False, which='minor', length=3, width=2)

        autoscale_y(ax, margin=0.1)
        plt.xlabel('Width [um]')
        plt.ylabel('Doping Concentration [cm]$^{-3}$')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-profileCV.pdf')
        plt.clf()
        plt.close()

class Site(object):
    def __init__(self, path, site, wafers, area=None):
        self.path = path
        self.site = site
        self.wafers = wafers
        self.area = area

    def show(self, IV=False, CV=False):
        print('Site:', site)
        for wafer in self.wafers:
            print('Wafer:', wafer)
            theWafer = Wafer(self.path, self.site, wafer, self.area)
            if IV:
                theWafer.plot_IVs()
            if CV:
                theWafer.plot_CVs()

'''
# December 2018
directory = '/Users/zschillaci/Google Drive/Sharable Silicon R&D Info/Measurements/JSI_Irradiation_Nov-Dec2018/'

### CNM-AIDA ###

path = directory + 'CNM-aida-sensors/'
site = 'CNM'
wafers = ['W3', 'W5', 'W7', 'W9', 'W11', 'W14']
area = (1.0 * 10**-3)**2  # 1.0mm by 1.0mm

CNM = Site(path, site, wafers, area)
CNM.show(IV=True, CV=True)

### HAMA ###

path = directory + 'hamamatsu-sensors/'
site = 'HAMA'
wafers = ['EXX28995-WNo8', 'EDX28976-WNo2']
area = (1.3 * 10**-3)**2 # 1.3mm by 1.3mm

HAMA = Site(path, site, wafers, area)
HAMA.show(IV=True, CV=True)
'''

# May 2019
directory = '/Users/zschillaci/Google Drive/LGAD Backup/'

path = directory + 'hamamatsu-sensors/'
site = 'HAMA'
wafers = ['EDX30329-WNo12',     'EXX30330-WNo11',       'EXX30330-WNo12',
          'EXX30327-WNo1',      'EXX30327-WNo3',        'EXX30327-WNo4',
          'EXX30327-WNo5',      'EXX30327-WNo6',        'EXX30328-WNo1',
          'EXX30328-WNo2',      'EXX30328-WNo3',        'EXX30328-WNo4',
          'EXX30328-WNo5',      'EDX30329-WNo9',        'EDX30329-WNo9-5x5',
          'EDX30329-WNo10',     'EDX30329-WNo11',       'EDX30329-WNo11-2x2',
          'EDX30329-WNo11-5x5', 'EXX30327-WNo1-5x5',    'EXX30327-WNo4-2x2',
          'EXX30327-WNo4-5x5',  'EXX30327-WNo5-5x5',    'EXX30327-WNo6-5x5',
          'EXX30328-WNo1-5x5',  'EXX30328-WNo2-5x5',    'EXX30328-WNo3-2x2',
          'EXX30328-WNo3-5x5',  'EXX30328-WNo5-5x5-P4', 'EXX30328-WNo5-5x5-P5',
          'EXX30330-WNo11-5x5', 'EXX30330-WNo12-5x5',   'EXX30330-WNo13',
          'EXX30330-WNo14',     'EXX30330-WNo14-5x5']

area = None

HAMA = Site(path, site, wafers, area)
HAMA.show(IV=True)
