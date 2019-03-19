#!/usr/bin/env /Users/zschillaci/Software/miniconda3/envs/pyenv/bin/python
import os
import glob
import collections
import numpy as np
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
            if ('log' not in afile):
                index = len(afile) - afile[::-1].find('/')

                infiles.append(afile)
                sensors.append(afile[index : -4])
                
    return infiles, sensors

def setYlim(ylim,yticks):
    #Set y-limits of plot based upon min. and max. of plots
    tick = abs(yticks[0][1] - yticks[0][0])

    low = ylim[0] - 0.15 * tick
    high = ylim[-1] + 1.25 * tick

    plt.ylim(low, high)

def calculate_width(C, A):
    #C = eps * A / d
    return (0 if (C == 0) else ((eps_si * eps0 * A) / C))

def analyze_CV(CV, area):
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

    return width, profile

class Sensor(object):
    def __init__(self, infile, area):
        self.infile = infile
        self.area = area

    def load_IV(self):
        #voltage, pad, ring, current
        data = np.loadtxt(self.infile, delimiter=',', skiprows=2)

        voltage, pad, ring, current = [], [], [], []
        for row in data:
            voltage.append(abs(row[0]))
            pad.append(i_scale * row[1])
            ring.append(i_scale * row[2])
            current.append(i_scale * row[3])

        IV = {'voltage': voltage, 'pad': pad,
              'ring': ring, 'current': current}

        return IV
        
    def load_CV(self):
        #voltage, current, capacitance, dvalue
        data = np.loadtxt(self.infile, delimiter=',', skiprows=1)

        voltage, current, capacitance, dvalue = [], [], [], []
        invsquarecapacitance = []
        for n, row in enumerate(data):
            voltage.append(abs(row[0]))
            current.append(row[1])
            capacitance.append(c_scale * row[2])
            dvalue.append(row[3])

            invsquarecapacitance.append(0 if (capacitance[n] == 0) else (1 / capacitance[n]**2))

        CV = {'voltage': voltage, 'current': current,
              'capacitance': capacitance, 'dvalue': dvalue,
              'invsquarecapacitance': invsquarecapacitance}
        
        width, profile = analyze_CV(CV, self.area)
        CV['width'] = width
        CV['profile'] = profile

        return CV
    
class Wafer(object):
    def __init__(self, path, site, wafer, area):
        self.path = path
        self.site = site
        self.wafer = wafer
        self.area = area
        
        self.inpath_to_IV = self.path + self.wafer + '/IV/'
        self.inpath_to_CV = self.path + self.wafer + '/CV/'

        # self.outpath_to_IV = './' + self.site + '/IV/'
        # self.outpath_to_CV = './' + self.site + '/CV/'

        self.outpath = self.path + 'Results/'
        mkdir(self.outpath)

        self.outpath_to_IV = self.outpath + '/IV/'
        self.outpath_to_CV = self.outpath + '/CV/'

        self.load_IVs()
        self.load_CVs()

    def load_IVs(self):
        infiles, sensors = get_list_of_infiles(self.inpath_to_IV)

        self.IVs = collections.OrderedDict()
        for infile, sensor in zip(infiles, sensors):
            self.IVs[sensor] = Sensor(infile, self.area).load_IV()

    def load_CVs(self):
        infiles, sensors = get_list_of_infiles(self.inpath_to_CV)

        self.CVs = collections.OrderedDict()
        for infile, sensor in zip(infiles, sensors):
            self.CVs[sensor] = Sensor(infile, self.area).load_CV()

    def plot_IVs(self):
        mkdir(self.outpath_to_IV)

        ### pad ###
        plt.figure('pad current vs. voltage' , figsize=(12,8))
        for n, (sensor, IV) in enumerate(self.IVs.items()):
            plt.plot(IV['voltage'], IV['pad'], 'o--', label=sensor, color=colors[n])

        if (self.site == 'HAMA'):
            plt.yscale('log')

        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('I$_{Pad}$ [uA]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)
     
        plt.savefig(self.outpath_to_IV + '/' + self.wafer + '-padIV.pdf')
        plt.clf()
        plt.close()

        ### ring ###
        plt.figure('ring current vs. voltage' , figsize=(12,8))
        for n, (sensor, IV) in enumerate(self.IVs.items()):
            plt.plot(IV['voltage'], IV['ring'], 'o--', label=sensor, color=colors[n])

        plt.yscale('log')

        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('I$_{Ring}$ [uA]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_IV + '/' + self.wafer + '-ringIV.pdf')
        plt.clf()
        plt.close()

    def plot_CVs(self):
        mkdir(self.outpath_to_CV)

        plt.figure('cv analysis' , figsize=(12,8))

        ### capacitance ###
        plt.subplot(111)
        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['voltage'], CV['capacitance'], 'o--', label=sensor, color=colors[n])

        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('C [pF]')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-CV.pdf')
        plt.clf()
        plt.close()

        plt.figure('cv analysis', figsize=(12, 8))

        ### invsquarecapacitance ###
        plt.subplot(111)
        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['voltage'], CV['invsquarecapacitance'], 'o--', label=sensor, color=colors[n])

        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('C$^{-2}$ [pF]$^{-2}$')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-invsquareCV.pdf')
        plt.clf()
        plt.close()

        ### dvalue ###
        '''
        plt.subplot(111)
        for n, (sensor, CV) in enumerate(self.CVs.items()):
            plt.plot(CV['voltage'], CV['dvalue'], 'o--', label=sensor, color=colors[n])

        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('V$_{Rev}$ [V]')
        plt.ylabel('D')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)
        '''

        fig = plt.figure('cv analysis', figsize=(12, 8))

        ### profile ###
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

        plt.ylim(-0.25 * 10**16, 3.0 * 10**16)
        setYlim(plt.ylim(), plt.yticks())
        plt.xlabel('Width [um]')
        plt.ylabel('Doping Concentration [cm]$^{-3}$')
        plt.legend(loc=9, ncol=3)
        plt.grid(False)

        plt.savefig(self.outpath_to_CV + '/' + self.wafer + '-profileCV.pdf')
        plt.clf()
        plt.close()

class Site(object):
    def __init__(self, path, site, wafers, area):
        self.path = path
        self.site = site
        self.wafers = wafers
        self.area = area

    def show(self):
        for wafer in self.wafers:
            theWafer = Wafer(self.path, self.site, wafer, self.area)
            theWafer.plot_IVs()
            theWafer.plot_CVs()

# directory = '/Users/zschillaci/Google Drive/TANDEM - Diodes for dosimeters and LGADs - Sept 2018/JSI_Irradiation/'
directory = '/Users/zschillaci/Google Drive/Sharable Silicon R&D Info/Measurements/JSI_Irradiation_Nov-Dec2018/'

### CNM-AIDA ###

path = directory + 'CNM-aida-sensors/'
site = 'CNM'
wafers = ['W3', 'W5', 'W7', 'W9', 'W11', 'W14']
area = (1.0 * 10**-3)**2  # 1.0mm by 1.0mm

CNM = Site(path, site, wafers, area)
CNM.show()

### HAMA ###

path = directory + 'hamamatsu-sensors/'
site = 'HAMA'
wafers = ['EXX28995-WNo8', 'EDX28976-WNo2']
area = (1.3 * 10**-3)**2 # 1.3mm by 1.3mm

HAMA = Site(path, site, wafers, area)
HAMA.show()

