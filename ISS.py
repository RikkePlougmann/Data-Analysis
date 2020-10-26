import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.stats import linregress

class Experiment(object):
    """Load an ISS experiment exported as text or VAMAS file.

Author: Jakob Ejler Sorensen
Version: 3.0
Date: 2017 June 23
    """

    def __init__(self, filename, mass=4, theta=146.7, E0=1000):
        """Initialize the class"""
        # Constants
        self.settings = dict()
        self.settings['mass'] = mass
        self.settings['theta'] = theta
        self.settings['E0'] = E0

        # Initialize variables
        self.energy = dict()
        self.cps = dict()
        self.dwell = dict()
        self.mode = dict()
        self.mode_value = dict()
        self.filename = filename

        # Read data from textfile:
        if filename.endswith('.txt'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            self.format = 'Text file'

            start_points = [i for i in range(len(lines)) if lines[i] == 'Energy\tCounts\r\n']
            self.scans = len(start_points)
            if self.scans == 0:
                raise ImportError('File apparently empty!')

            if lines[0].lower().startswith('note'):
                self.note = lines[0].split('=')[1].lstrip(' ')
            # Copy data points
            counter = 0
            for start in start_points:
                line = lines[start-4].split('\t')
                for tab in range(len(line)):
                    if line[tab].lower() == 'dwell':
                       self.dwell[counter] = float(lines[start-3].split('\t')[tab])
                if not start == start_points[-1]:
                    interval = range(start+1, start_points[counter+1]-4)
                else:
                    interval =  range(start+1, len(lines))
                self.energy[counter] = np.zeros(len(interval))
                self.cps[counter] = np.zeros(len(interval))
                counter_inner = 0
                for index in interval:
                    line = lines[index].rstrip('\r\n').split('\t')
                    self.energy[counter][counter_inner] = float(line[0])
                    self.cps[counter][counter_inner] = float(line[1])
                    counter_inner += 1
                self.cps[counter] = self.cps[counter]/self.dwell[counter]
                counter += 1
        # Read data from old VAMAS block file
        elif filename.endswith('.vms'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
	    # Old format:
            if lines[6].lower().startswith('experiment type'):
                self.format = 'Old VAMAS'
                print('Loading file: ' + filename)
                blocks_4 = [i for i in range(len(lines)) if (lines[i] == '-1\r\n') and (lines[i+1].lower() == 'kinetic energy\r\n')]
                blocks_2_ISS = [i for i in range(len(lines)) if (lines[i] == 'ISS\r\n') and (lines[i+1] == '\r\n')]
                print(lines[9].rstrip('\r\n'))
                self.scans = len(blocks_4)
                if len(blocks_4) == int(lines[9].rstrip('\r\n')) and len(blocks_4) == len(blocks_2_ISS):
                    self.scans = len(blocks_4)
                else:
                    raise ImportError('Error: Identified {} "Block 4", {} "Block 2", but "Block 1" says: {}'.format(len(blocks_4), len(blocks_2_ISS), int(lines[9].rstrip('\r\n'))))

                # Copy data points
                self.note = dict()
                for counter in range(len(blocks_4)):
                    i = blocks_4[counter]
                    if not len(lines[blocks_2_ISS[counter] - 1]) == 5:
                        self.note[counter] = lines[blocks_2_ISS[counter] - 1].rstrip('\r\n')
                    else:
                        self.note[counter] = ''
                    self.mode[counter] = lines[i-11].rstrip('\r\n')
                    self.mode_value[counter] = float( lines[i-10].rstrip('\r\n') )
                    self.dwell[counter] = float( lines[i+9].rstrip('\r\n') )
                    data_points = int( lines[i+16] )
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float( lines[i+4].rstrip('\r\n') )
                    E_start = float( lines[i+3].rstrip('\r\n') )
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float( lines[i+19+counter_inner] )/self.dwell[counter]
		    #self.note[counter] = ''
                print(self.energy.viewkeys())
                print('Comments: {}'.format(self.note))
                print('Dwell time: {}'.format(self.dwell))
                print('Modes: {}'.format(self.mode))
                print('Mode values: {}'.format(self.mode_value))
	    # New format:
            if lines[6].lower().startswith('created with'):
                self.format = 'New VAMAS'
                #ENDING = '_1-Detector_Region.vms'
                #filen = filename.rstrip(ENDING)
                #counter = 0
                #while True:
            #        try:
            #            f = open(filen + '--' + str(counter+1) + ENDING)
        #                counter += 1
    #                except:
#                        #print('{} files detected of series:'.format(counter))
#                        #print('* ' + filen + '--' + str(1) + ENDING + ' *')
#                        COUNTER = counter
#                        break
                # Open filename with ISS data
                #self.scans = COUNTER
                COUNTER = 1
                for counter in range(COUNTER):
                    #new_filename = filen + '--' + str(counter+1) + ENDING
                    new_filename = filename
                    f = open(new_filename, 'r')
                    lines = f.readlines()
                    f.close()
                    print('Loading file: ' + new_filename)
                    blocks_4 = [i for i in range(len(lines)) if (lines[i] == '-1\n') and (lines[i+1].lower() == 'kinetic energy\n')]
                    print(lines[9].rstrip('\r\n'))
                    print(blocks_4)
                    if len(blocks_4) > 1:
                        print('*** Interesting! More than 1 scan has been detected in above single file!')
                    # Copy data points
                    i = blocks_4[0]
                    self.mode[counter] = lines[i-11].rstrip('\r\n')
                    self.mode_value[counter] = float( lines[i-10].rstrip('\r\n') )
                    self.dwell[counter] = float( lines[i+9].rstrip('\r\n') )
                    data_points = int( lines[i+16] )
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float( lines[i+4].rstrip('\r\n') )
                    E_start = float( lines[i+3].rstrip('\r\n') )
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float( lines[i+19+counter_inner] )/self.dwell[counter]



    def ConvertEnergy(self, mass):
        """Converts a measured energy to mass of surface atom
corresponding the settings stored in the experiment.
        """
        angle = self.settings['theta'] * np.pi/180
        return self.settings['E0'] * ( (self.settings['mass']*np.cos(angle) + np.sqrt(mass**2 - self.settings['mass']**2*np.sin(angle)**2))/(mass + self.settings['mass']) )**2





    def PlotAllScans(self, exclude=[None], color=None):
        """Plot all elements in file in single figure."""
        selection = [i for i in range(self.scans) if not i in exclude]
        if not color:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i])
        else:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i], color=color)
        plt.xlabel('Kinetic energy (eV)')
        plt.ylabel('Counts per second')





    def Normalize(self, interval, exclude=[None], unit='Mass'):
        """Normalize to highest value in interval=[value1, value2]"""
        if type(interval) == int:
            self.normalization_criteria = interval
        elif type(interval) == str:
            if interval == 'Au':
                self.normalization_criteria = 196.
        if not type(interval) == list:
            interval = [0, 0]
            interval[0] = self.ConvertEnergy(self.normalization_criteria) - 10
            interval[1] = self.ConvertEnergy(self.normalization_criteria) + 10
            selection = [i for i in range(self.scans) if (not i in exclude) and (not interval[0] > max(self.energy[i])) and (not interval[1] < min(self.energy[i]))]
        for __counter in selection:
            range_1 = np.where(self.energy[__counter] < interval[1])[0]
            range_2 = np.where(self.energy[__counter] > interval[0])[0]
            energy_range = np.intersect1d(range_1, range_2)
            value = max( self.cps[__counter][energy_range] )
            self.cps[__counter] = self.cps[__counter]/value




    def AddMassLines(self, masses, offset=0, color='k'):
        """Add vertical lines for mass references."""
        energies = self.ConvertEnergy(np.array(masses))
        ax = plt.gca()
        [x1, x2, y1, y2] = ax.axis()
        for i in range(len(energies)):
            ax.axvline(x=energies[i]-offset, ymin=0, ymax=1, linestyle='dotted', color=color)
            ax.text(float(energies[i])/x2, 0.95, 'm-{}'.format(masses[i]), transform=ax.transAxes)

    def AddRegions(self, ax=plt.gca()):
        """Add regions indicating the whereabouts of 3d, 4d, 5d metals and the lanthanides and actinides."""
        d3 = [45, 65]
        d4 = [89, 112]
        d5 = [178, 201]
        lant = [139, 175]
        act = [227, 260]
        for i in [d3, d4, d5]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]), color='k', alpha=0.2)
        for i in [lant, act]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]), color='y', alpha=0.2)
