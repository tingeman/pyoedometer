
# embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases
#
# Copyright (C) 2005 Florent Rougon
#               2006 Darren Dale
#               2015 Jens H Nielsen
#
# This file is an example program for matplotlib. It may be used and
# modified with no restriction; raw copies as well as modified versions
# may be distributed without limitation.


# SOME GOOD INFORMATION:
# https://helpful.knobs-dials.com/index.php/Qt_and_PyQt_notes
# https://stackoverflow.com/questions/5304570/how-to-save-new-data-in-tree-model-view
# http://rowinggolfer.blogspot.dk/2010/05/qtreeview-and-qabractitemmodel-example.html
# https://stackoverflow.com/questions/31342228/pyqt-tree-widget-adding-check-boxes-for-dynamic-removal
# https://stackoverflow.com/questions/31600462/spanning-multiple-columns-using-qtreeview-and-qabstractitemmodel

# Integrating Matplotlib graphs
# https://sukhbinder.wordpress.com/2013/12/16/simple-pyqt-and-matplotlib-example-with-zoompan/
# https://pythonspot.com/en/pyqt5-matplotlib/

# pyqt threading
# https://nikolak.com/pyqt-threading-tutorial/

from __future__ import unicode_literals
import pdb
import sys
import os
import random
import datetime as dt
import numpy as np
import pandas as pd
import matplotlib

# Make sure that we are using QT5
import yaml

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib import gridspec


from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import treemodel as tm
from sqrtlogscale import SqrtLogScale

progname = os.path.basename(sys.argv[0])
progversion = "0.1"


# ----------------------------------------------------------------------------------------------------------------------
# Data definitions
# ----------------------------------------------------------------------------------------------------------------------

# TEST INFORMATION
# applied load [kPa]
# room set-point temperature [C]
# Times are server time (UTC, no dst)

history = []
hist_fname = 'hist.yml'
data_path = './sample_1_2'

# history = [
#     {'step':  1, 'date': dt.date(2017, 10,  3), 'load':    1.5, 'temp': -10, 'time1': dt.time(16,  0,  0), 'time2': dt.time(16,  0,  0), 'comment': 'Set the sample and LVDTs around 16:00'},
#     {'step':  2, 'date': dt.date(2017, 10,  4), 'load':   10.0, 'temp': -10, 'time1': dt.time(14, 18, 46), 'time2': dt.time(14, 20,  0), 'comment': ''},
#     {'step':  3, 'date': dt.date(2017, 10,  9), 'load':   20.0, 'temp': -10, 'time1': dt.time(12, 53,  1), 'time2': dt.time(12, 53,  6), 'comment': ''},
#     {'step':  4, 'date': dt.date(2017, 10, 10), 'load':   50.0, 'temp': -10, 'time1': dt.time(12, 49, 46), 'time2': dt.time(12, 49, 56), 'comment': ''},
#     {'step':  5, 'date': dt.date(2017, 10, 11), 'load':  100.0, 'temp': -10, 'time1': dt.time(13,  4, 12), 'time2': dt.time(13,  4, 38), 'comment': ''},
#     {'step':  6, 'date': dt.date(2017, 10, 12), 'load':  200.0, 'temp': -10, 'time1': dt.time(13, 20, 11), 'time2': dt.time(13, 20, 24), 'comment': ''},
#     {'step':  7, 'date': dt.date(2017, 10, 13), 'load':  200.0, 'temp':  -9, 'time1': dt.time(13, 58, 59), 'time2': dt.time(13, 58, 59), 'comment': ''},
#     {'step':  8, 'date': dt.date(2017, 10, 14), 'load':  200.0, 'temp':  -8, 'time1': dt.time(14, 20, 39), 'time2': dt.time(14, 20, 39), 'comment': ''},
#     {'step':  9, 'date': dt.date(2017, 10, 15), 'load':  200.0, 'temp':  -7, 'time1': dt.time(14, 25, 51), 'time2': dt.time(14, 25, 51), 'comment': ''},
#     {'step': 10, 'date': dt.date(2017, 10, 16), 'load':  200.0, 'temp':  -6, 'time1': dt.time(13, 56, 13), 'time2': dt.time(13, 56, 13), 'comment': ''},
#     {'step': 11, 'date': dt.date(2017, 10, 17), 'load':  200.0, 'temp':  -5, 'time1': dt.time(14, 54, 16), 'time2': dt.time(14, 54, 16), 'comment': ''},
#     {'step': 12, 'date': dt.date(2017, 10, 18), 'load':  200.0, 'temp':  -4, 'time1': dt.time(15,  7,  9), 'time2': dt.time(15,  7,  9), 'comment': ''},
#     {'step': 13, 'date': dt.date(2017, 10, 19), 'load':  200.0, 'temp':  -3, 'time1': dt.time(15, 38, 37), 'time2': dt.time(15, 38, 37), 'comment': ''},
#     {'step': 14, 'date': dt.date(2017, 10, 23), 'load':  200.0, 'temp':  -2, 'time1': dt.time(14, 41, 49), 'time2': dt.time(14, 41, 49), 'comment': ''},
#     {'step': 15, 'date': dt.date(2017, 10, 26), 'load':  200.0, 'temp':  -1, 'time1': dt.time(15, 50, 38), 'time2': dt.time(15, 50, 38), 'comment': ''},
#     {'step': 16, 'date': dt.date(2017, 10, 27), 'load':  200.0, 'temp':   0, 'time1': dt.time(14, 14, 39), 'time2': dt.time(14, 14, 39), 'comment': ''},
#     {'step': 17, 'date': dt.date(2017, 10, 28), 'load':  200.0, 'temp':   5, 'time1': dt.time(14, 47, 33), 'time2': dt.time(14, 47, 33), 'comment': ''},
#     {'step': 18, 'date': dt.date(2017, 10, 29), 'load':  100.0, 'temp':   5, 'time1': dt.time(16, 48, 44), 'time2': dt.time(16, 49,  0), 'comment': ''},
#     {'step': 19, 'date': dt.date(2017, 10, 30), 'load':   50.0, 'temp':   5, 'time1': dt.time(17, 55, 28), 'time2': dt.time(17, 55, 38), 'comment': ''},
#     {'step': 20, 'date': dt.date(2017, 10, 31), 'load':   10.0, 'temp':   5, 'time1': dt.time(11, 12, 56), 'time2': dt.time(11, 13,  5), 'comment': ''},
#     {'step': 21, 'date': dt.date(2017, 11,  1), 'load':   20.0, 'temp':   5, 'time1': dt.time(11, 38, 23), 'time2': dt.time(11, 38, 26), 'comment': ''},
#     {'step': 22, 'date': dt.date(2017, 11,  2), 'load':   50.0, 'temp':   5, 'time1': dt.time(10, 57, 47), 'time2': dt.time(10, 57, 51), 'comment': ''},
#     {'step': 23, 'date': dt.date(2017, 11,  2), 'load':  120.0, 'temp':   5, 'time1': dt.time(17,  1, 50), 'time2': dt.time(17,  2, 16), 'comment': 'Check times'},
#     {'step': 24, 'date': dt.date(2017, 11,  3), 'load':  300.0, 'temp':   5, 'time1': dt.time(15, 47, 42), 'time2': dt.time(15, 47, 48), 'comment': ''},
#     {'step': 25, 'date': dt.date(2017, 11,  4), 'load':  600.0, 'temp':   5, 'time1': dt.time(19, 11, 35), 'time2': dt.time(19, 11, 41), 'comment': ''},
#     {'step': 26, 'date': dt.date(2017, 11,  6), 'load': 1200.0, 'temp':   5, 'time1': dt.time(14, 36, 47), 'time2': dt.time(14, 36, 53), 'comment': ''},
#     {'step': 27, 'date': dt.date(2017, 11,  7), 'load': 2400.0, 'temp':   5, 'time1': dt.time(12,  0,  0), 'time2': dt.time(12,  0,  0), 'comment': ''},
# ]

# Sample info
sample_1_h0 = 20.50  # mm
sample_2_h0 = 21.12  # mm

# LVDT calibrations
lvdt_1_a = np.mean([-992.649, -1005.437])
lvdt_1_b = np.mean([4339.163, 4404.368])
lvdt_1_h0 = 9.0234  # mm

lvdt_2_a = np.mean([-985.457, -1004.581])
lvdt_2_b = np.mean([4125.023, 4227.782])
lvdt_2_h0 = 9.0234  # mm

lvdt_3_a = np.mean([-985.514, -1000.657])
lvdt_3_b = np.mean([4816.420, 4888.629])
lvdt_3_h0 = 0  # mm

lvdt_labels = ['Sample 1', 'sample 2']

# HOBO CALIBRATIONS:   T = T_meas + delta_T
cal_hobo_1 = -0.02
cal_hobo_2 = -0.02
cal_hobo_3 = 0.00

hobo_labels = ['Isolated chamber', 'Frost chamber', 'Computer room']

# PT100 calibrations:  multiplicative factor
cal_pt100_1 = 1.000320692
cal_pt100_2 = 1.001201753

pt100_labels = ['Sample 1', 'sample 2']


class PT100_converter():
    # Callendar-Van Dusen Coefficients
    # FOR TEMPERATURES < 0
    g = 1.79090E+00
    j = 2.55819E+02
    i = 9.14550E+00
    h = -2.92363E+00

    # FOR TEMPERATURES > = 0
    d = -2.310000E-06
    e = 1.758481E-05
    a = 3.908300E-03
    f = -1.155000E-06

    def __init__(self):
        pass

    @classmethod
    def calc_temperature(cls, X):
        return np.where(X < 1,
                        cls.g * (X - 1) ** 4 + cls.h * (X - 1) ** 3 + cls.i * (X - 1) ** 2 + cls.j * (X - 1),
                        (np.sqrt(cls.d * X + cls.e) - cls.a) / cls.f)




class MyMplAxCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax1 = fig.add_subplot(111)

        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyMplTwinAxCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax1 = fig.add_subplot(111)
        self.ax2 = self.ax1.twinx()

        self.compute_initial_figure()

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyMplTrippleAxCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)

        ax1 = fig.add_subplot(311)
        ax1a = ax1.twinx()
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313, sharex=ax1)

        self.axes = [ax1, ax1a, ax2, ax3]

        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plot(self, step=0, start_dt=None, end_dt=None):

        def get_closest(array, values):
            # make sure array is a numpy array
            array = np.array(array)

            # get insert positions
            idxs = np.searchsorted(array, values, side="left")

            # find indexes where previous index is closer
            prev_idx_is_less = ((idxs == len(array)) | (np.fabs(values - array[np.maximum(idxs - 1, 0)]) < np.fabs(
                values - array[np.minimum(idxs, len(array) - 1)])))

            idxs[prev_idx_is_less] -= 1

            return idxs

        def reduce(time_array):
            log_t0 = np.log10(time_array[1])
            log_tmax = np.log10(time_array[-1])
            nomvals = 10 ** np.linspace(log_t0, log_tmax, np.floor((log_tmax - log_t0) * 20))
            idx = get_closest(time_array, nomvals)
            return idx


        for ax in self.axes:
            ax.clear()

        mw = self.window()
        lvdts = mw.data['lvdts']
        pt100T = mw.data['pt100T']
        hoboT = mw.data['hoboT']

        id = step
        hist = history[step]
        if hist['date'] is None: return  # Only plot if the step has a date-value assigned

        # Get the start date and time
        if start_dt is None:
            step_starttime = dt.datetime.combine(hist['date'], hist['time1'])
        else:
            step_starttime = start_dt

        if end_dt is None:
            if id + 1 >= len(history) or history[id + 1]['time1'] is None:
                # If this is the last step, use the last time in the dataset
                step_endtime = max(max(lvdts.index), max(lvdts.index))
            else:
                # If this is not the last step, use the start time for the next step minus 1 sec
                step_endtime = dt.datetime.combine(history[id + 1]['date'], history[id + 1]['time1'])
                step_endtime -= dt.timedelta(seconds=1)
        else:
            step_endtime = end_dt

        # select data for plotting: only data within the present load/temperature step start and end times
        lvdt_step = lvdts[(lvdts.index >= step_starttime) & (lvdts.index <= step_endtime)]
        hobo_step = hoboT[(hoboT.index >= step_starttime) & (hoboT.index <= step_endtime)]
        pt100_step = pt100T[(pt100T.index >= step_starttime) & (pt100T.index <= step_endtime)]

        print('lvdt_step length: {}'.format(len(lvdt_step)))

        if not lvdt_step.empty:
            lvdt_step['tmin'] = (lvdt_step.index - lvdt_step.index[0]).total_seconds() / 60.
        if not pt100_step.empty:
            pt100_step['tmin'] = (pt100_step.index - pt100_step.index[0]).total_seconds() / 60.
        if not hobo_step.empty:
            hobo_step['tmin'] = (hobo_step.index - hobo_step.index[0]).total_seconds() / 60.

        # reduce the dataset
        if mw.data['plot_options']['reduce']:
            if not lvdt_step.empty:
                idx = np.unique(reduce(lvdt_step['tmin']))
                lvdt_step = lvdt_step.ix[idx]
#                if mw.data['plot_options']['scale'] == 'log':
#                    lvdt_step.index = tmin[idx]

            if not pt100_step.empty:
                idx = np.unique(reduce(pt100_step['tmin']))
                pt100_step = pt100_step.ix[idx]
#                if mw.data['plot_options']['scale'] == 'log':
#                    pt100_step.index = tmin[idx]

            if not hobo_step.empty:
                idx = np.unique(reduce(hobo_step['tmin']))
                hobo_step = hobo_step.ix[idx]
#                if mw.data['plot_options']['scale'] == 'log':
#                    hobo_step.index = tmin[idx]

        print('lvdt_step length: {}'.format(len(lvdt_step)))
        print('tmin:  {0:.3f} - {1:.3f}'.format(lvdt_step['tmin'][0],lvdt_step['tmin'][-1]))
        # Plot the first lvdt
        
        if not pt100_step.empty:
            if mw.data['plot_options']['scale'] == 'lin':
                pt100_step['T(1)'].plot(style='.r', ax=self.axes[2], label=pt100_labels[0], zorder=10)
                pt100_step['T(2)'].plot(style='.b', ax=self.axes[2], label=pt100_labels[1], zorder=10)
            else:
                pt100_step['T(1)'].plot(x=pt100_step['tmin'], style='.r', ax=self.axes[2], label=pt100_labels[0], zorder=10)
                pt100_step['T(2)'].plot(x=pt100_step['tmin'], style='.b', ax=self.axes[2], label=pt100_labels[1], zorder=10)

        if not hobo_step.empty:
            if mw.data['plot_options']['scale'] == 'lin':
                hobo_step['T(1)'].plot(style='.g', ax=self.axes[2],  label=hobo_labels[0], zorder=2)
                hobo_step['T(2)'].plot(style='.y', ax=self.axes[2],  label=hobo_labels[1], zorder=1)
                hobo_step['T(3)'].plot(style='.k', ax=self.axes[3],  label=hobo_labels[1], zorder=1)
            else:
                hobo_step['T(1)'].plot(x=hobo_step['tmin'], style='.g', ax=self.axes[2], label=hobo_labels[0], zorder=2)
                hobo_step['T(2)'].plot(x=hobo_step['tmin'], style='.y', ax=self.axes[2], label=hobo_labels[1], zorder=1)
                hobo_step['T(3)'].plot(x=hobo_step['tmin'], style='.k', ax=self.axes[3], label=hobo_labels[1], zorder=1)

        if not lvdt_step.empty:
            if mw.data['plot_options']['scale'] == 'lin':
                lvdt_step['Volt(1)'].plot(style='.r', ax=self.axes[0], label=lvdt_labels[0], zorder=10)
                lvdt_step['Volt(2)'].plot(style='.b', ax=self.axes[1], label=lvdt_labels[1], zorder=10)
            else:
                lvdt_step['Volt(1)'].plot(x=lvdt_step['tmin'], style='.r', ax=self.axes[0], label=lvdt_labels[0], zorder=10)
                lvdt_step['Volt(2)'].plot(x=lvdt_step['tmin'], style='.b', ax=self.axes[1], label=lvdt_labels[1], zorder=10)

        #pdb.set_trace()

        print(mw.data['plot_options']['scale'])
        if mw.data['plot_options']['scale'] == 'log':
            for ax in self.axes:
                ax.set_xscale('log')
            print('set logx scale')
        elif mw.data['plot_options']['scale'] == 'sqrt':
            for ax in self.axes:
                ax.set_xscale('sqrtlog', intersect=1.)
            print('set sqrtlogx scale')

        self.setScale()

        # invert the direction of both y-axes
        self.axes[0].invert_yaxis()
        self.axes[1].invert_yaxis()

        self.axes[0].set_ylabel('Strain, Sample 1 [mm]')
        self.axes[1].set_ylabel('Strain, Sample 2 [mm]')
        self.axes[2].set_ylabel('Temperature [C]')
        self.axes[3].set_ylabel('Temperature [C]')

        self.draw()

    def setScale(self):
        mw = self.window()
        print('In setScale method')
        print(mw.data['plot_options']['scale'])
        if mw.data['plot_options']['scale'] == 'log':
            for ax in self.axes:
                ax.set_xscale('log')
                ax.set_xlim(auto=True)
            print('set logx scale')
        elif mw.data['plot_options']['scale'] == 'sqrt':
            for ax in self.axes:
                ax.set_xscale('sqrtlog', intersect=1.)
                ax.set_xlim(left=0.)
            print('set sqrtlogx scale')
        elif mw.data['plot_options']['scale'] == 'lin':
            for ax in self.axes:
                #ax.set_xscale('linear')
                ax.set_xlim(auto=True)
            print('set lin scale')

        self.draw()

    def home(self):
        self.window().zp_toolbar.home()

    def zoom(self):
        self.window().zp_toolbar.zoom()

    def pan(self):
        self.window().zp_toolbar.pan()


class LVDTCanvas(MyMplTwinAxCanvas):
    """Simple canvas with a sine plot."""

    def plot(self, lvdts, step=0):
        #pdb.set_trace()

        self.ax1.clear()
        self.ax2.clear()

        id = step
        hist = history[step]
        if hist['date'] is None: return  # Only plot if the step has a date-value assigned

        # if not (hist['step'] == 15): continue

        # Get the start date and time
        step_starttime = dt.datetime.combine(hist['date'], hist['time1'])

        # Get the end date and time
        if id + 1 > len(history) or history[id + 1]['time1'] is None:
            # If this is the last step, use the last time in the dataset
            step_endtime = max(max(lvdts.index), max(lvdts.index))
        else:
            # If this is not the last step, use the start time for the next step minus 1 sec
            step_endtime = dt.datetime.combine(history[id + 1]['date'], history[id + 1]['time1'])
            step_endtime -= dt.timedelta(seconds=1)

        print('step: {}, start: {}, end: {}'.format(step, step_starttime, step_endtime))

        # select data for plotting: only data within the present load/temperature step start and end times
        lvdt_step = lvdts[(lvdts.index >= step_starttime) & (lvdts.index <= step_endtime)]

        # Plot the first lvdt
        lvdt_step['Volt(1)'].plot(style='-r', ax=self.ax1, label=lvdt_labels[0], zorder=10)
        lvdt_step['Volt(2)'].plot(style='-b', ax=self.ax2, label=lvdt_labels[1], zorder=10)

        # invert the direction of both y-axes
        self.ax1.invert_yaxis()
        self.ax2.invert_yaxis()

        self.ax1.set_ylabel('Strain, Sample 1 [mm]')
        self.ax2.set_ylabel('Strain, Sample 2 [mm]')
        self.draw()


class TempCanvas1(MyMplAxCanvas):
    """Simple canvas with a sine plot."""
    def plot(self, pt100T, hoboT, step=0):

        self.ax1.clear()

        id = step
        hist = history[step]
        if hist['date'] is None: return  # Only plot if the step has a date-value assigned

        # if not (hist['step'] == 15): continue

        # Get the start date and time
        step_starttime = dt.datetime.combine(hist['date'], hist['time1'])

        # Get the end date and time
        if id + 1 > len(history) or history[id + 1]['time1'] is None:
            # If this is the last step, use the last time in the dataset
            step_endtime = max(max(pt100T.index), max(pt100T.index))
        else:
            # If this is not the last step, use the start time for the next step minus 1 sec
            step_endtime = dt.datetime.combine(history[id + 1]['date'], history[id + 1]['time1'])
            step_endtime -= dt.timedelta(seconds=1)

        # select data for plotting: only data within the present load/temperature step start and end times
        hobo_step = hoboT[(hoboT.index >= step_starttime) & (hoboT.index <= step_endtime)]
        pt100_step = pt100T[(pt100T.index >= step_starttime) & (pt100T.index <= step_endtime)]

        if not pt100_step.empty:
            pt100_step['T(1)'].plot(style='-r', ax=self.ax1, label=pt100_labels[0], zorder=10)
            pt100_step['T(2)'].plot(style='-b', ax=self.ax1, label=pt100_labels[1], zorder=10)
        if not hobo_step.empty:
            hobo_step['T(1)'].plot(style='-g', ax=self.ax1,  label=hobo_labels[0], zorder=2)
            hobo_step['T(2)'].plot(style='-y', ax=self.ax1,  label=hobo_labels[1], zorder=1)

        self.ax1.set_ylabel('Temperature [C]')
        self.draw()

class TempCanvas2(MyMplAxCanvas):
    """Simple canvas with a sine plot."""
    def plot(self, hoboT, step=0):

        self.ax1.clear()

        id = step
        hist = history[step]
        if hist['date'] is None: return  # Only plot if the step has a date-value assigned

        # if not (hist['step'] == 15): continue

        # Get the start date and time
        step_starttime = dt.datetime.combine(hist['date'], hist['time1'])

        # Get the end date and time
        if id + 1 > len(history) or history[id + 1]['time1'] is None:
            # If this is the last step, use the last time in the dataset
            step_endtime = max(max(hoboT.index), max(hoboT.index))
        else:
            # If this is not the last step, use the start time for the next step minus 1 sec
            step_endtime = dt.datetime.combine(history[id + 1]['date'], history[id + 1]['time1'])
            step_endtime -= dt.timedelta(seconds=1)

        # select data for plotting: only data within the present load/temperature step start and end times
        hobo_step = hoboT[(hoboT.index >= step_starttime) & (hoboT.index <= step_endtime)]

        if not hobo_step.empty:
            hobo_step['T(1)'].plot(style='-k', ax=self.ax1,  label=hobo_labels[2], zorder=10)

        self.ax1.set_ylabel('Temperature [C]')
        self.draw()


class MyDynamicMplCanvas(MyMplAxCanvas):
    """A canvas that updates itself every second with a new plot."""

    def __init__(self, *args, **kwargs):
        MyDynamicMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.ax1.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')

    def update_figure(self):
        # Build a list of 4 random integers between 0 and 10 (both inclusive)
        l = [random.randint(0, 10) for i in range(4)]
        self.ax1.cla()
        self.ax1.plot([0, 1, 2, 3], l, 'r')
        self.draw()


class GeoLocationWidget(QtWidgets.QWidget):

  __pyqtSignals__ = ("latitudeChanged(double)",
                     "longitudeChanged(double)")

  def __init__(self, parent = None):

     QtWidgets.QWidget.__init__(self, parent)

     latitudeLabel = QtWidgets.QLabel(self.tr("Latitude:"))
     self.latitudeSpinBox = QtWidgets.QDoubleSpinBox()
     self.latitudeSpinBox.setRange(-90.0, 90.0)
     self.latitudeSpinBox.setDecimals(5)

     longitudeLabel = QtWidgets.QLabel(self.tr("Longitude:"))
     self.longitudeSpinBox = QtWidgets.QDoubleSpinBox()
     self.longitudeSpinBox.setRange(-180.0, 180.0)
     self.longitudeSpinBox.setDecimals(5)

     layout = QtWidgets.QGridLayout(self)
     layout.addWidget(latitudeLabel, 0, 0)
     layout.addWidget(self.latitudeSpinBox, 0, 1)
     layout.addWidget(longitudeLabel, 1, 0)
     layout.addWidget(self.longitudeSpinBox, 1, 1)


class DataFilesTreeWidget(QtWidgets.QTreeView):
    def __init__(self, parent=None):
        if parent is None:
            QtWidgets.QTreeView.__init__(self)
        else:
            QtWidgets.QTreeView.__init__(self, parent)

        #model = tm.CheckableTreeModel(headers=['Name','Plot'])
        model = tm.TreeModel(headers=['Name','Plot'])
        it = model.addItem({'Name':'1', 'Plot':'ABC'})
        it2 = model.addItem({'Name': '1', 'Plot': 'DEF'}, it)

        self.setModel(model)
        self.setAlternatingRowColors(True)
        self.expandAll()

class DataSelectorWidget(QtWidgets.QWidget):
    def __init__(self, parent=None):
        if parent is None:
            QtWidgets.QWidget.__init__(self)
        else:
            QtWidgets.QWidget.__init__(self, parent)

        layout = QtWidgets.QVBoxLayout(self)

        self.cb = QtWidgets.QComboBox()
        self.cb.currentIndexChanged.connect(self.comboChange)
        layout.addWidget(self.cb)

        self.stdte = QtWidgets.QDateTimeEdit()
        self.stdte.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.stdte.setCurrentSectionIndex(5)
        self.stdte.dateTimeChanged.connect(self.dateTimeChange)
        layout.addWidget(self.stdte)

        self.etdte = QtWidgets.QDateTimeEdit()
        self.etdte.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.etdte.setCurrentSectionIndex(5)
        self.etdte.dateTimeChanged.connect(self.dateTimeChange)
        layout.addWidget(self.etdte)

        self.linrb = QtWidgets.QRadioButton('Linear(t)')
        self.linrb.setChecked(True)
        self.logrb = QtWidgets.QRadioButton('Log(t)')
        self.sqrtrb = QtWidgets.QRadioButton('Sqrt(t)')

        self.linrb.toggled.connect(lambda: self.radioBtnState(self.linrb))
        self.logrb.toggled.connect(lambda: self.radioBtnState(self.logrb))
        self.sqrtrb.toggled.connect(lambda: self.radioBtnState(self.sqrtrb))

        layout.addWidget(self.linrb)
        layout.addWidget(self.logrb)
        layout.addWidget(self.sqrtrb)

        self.redcb = QtWidgets.QCheckBox('Reduce data')
        self.redcb.setChecked(True)
        self.redcb.stateChanged.connect(lambda: self.cbToggle(self.redcb))
        layout.addWidget(self.redcb)

        self.setLayout(layout)
        #self.setWindowTitle("combo box demo")

    def radioBtnState(self, b):
        mw = self.window()
        if not b.isChecked():
            print('{} not checked!'.format(b.text()))
            return
        elif b.text() == "Linear(t)":
            print('Linear(t)')
            mw.data['plot_options']['scale'] = 'lin'
        elif b.text() == "Log(t)":
            print('Log(t)')
            mw.data['plot_options']['scale'] = 'log'
        elif b.text() == "Sqrt(t)":
            print('Sqrt(t)')
            mw.data['plot_options']['scale'] = 'sqrt'

        mw.widgets['pw'].canvas.setScale()

    def cbToggle(self, cb):
        if not cb.isChecked():
            print('{} not checked!'.format(cb.text()))
            self.window().data['plot_options']['reduce'] = False
        elif cb.text() == "Reduce data":
            print('Reduce data')
            self.window().data['plot_options']['reduce'] = True

    def setupContent(self):

        self.stdte.blockSignals(True)
        self.etdte.blockSignals(True)
        self.cb.blockSignals(True)

        for id, hist in enumerate(history):
            label = 'Step: {step}, Load: {load} kPa, T: {temp} C'.format(**hist)
            self.cb.addItem(label)
        self.cb.setCurrentIndex(0)

        if len(history) >= 1:
            self.stdte.setDateTime(dt.datetime.combine(history[0]['date'],history[0]['time1']))
        if len(history) >=2:
            self.etdte.setDateTime(dt.datetime.combine(history[1]['date'], history[1]['time1']))
        else:
            self.etdte.setDateTime(dt.datetime.combine(history[0]['date'], history[0]['time1'])+
                                   dt.timedelta(hours=24))

        self.stdte.blockSignals(False)
        self.etdte.blockSignals(False)
        self.cb.blockSignals(False)

    def comboChange(self, i):

        self.stdte.blockSignals(True)
        self.etdte.blockSignals(True)
        self.cb.blockSignals(True)

        if len(history) > i:
            self.stdte.setDateTime(dt.datetime.combine(history[i]['date'],history[i]['time1']))
        if len(history) > i+1:
            self.etdte.setDateTime(dt.datetime.combine(history[i+1]['date'], history[i+1]['time1']))
        else:
            self.etdte.setDateTime(dt.datetime.combine(history[i]['date'], history[i]['time1'])+
                                   dt.timedelta(hours=24))

        self.stdte.blockSignals(False)
        self.etdte.blockSignals(False)
        self.cb.blockSignals(False)

        mw = self.window()
        mw.widgets['pw'].canvas.plot(step=i)

    def dateTimeChange(self):
        mw = self.window()
        mw.widgets['pw'].canvas.plot(step=self.cb.currentIndex(),
                                     start_dt=self.stdte.dateTime().toPyDateTime(),
                                     end_dt=self.etdte.dateTime().toPyDateTime())


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Oedometer data processing")

        # File menu
        self.file_menu = QtWidgets.QMenu('&File', self)

        self.file_menu.addAction('&Open', self.fileOpen,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_O)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        # Help menu
        self.help_menu = QtWidgets.QMenu('&Help', self)
        self.help_menu.addAction('&About', self.about)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.data = dict()
        self.widgets = dict()

        # Create plot_widget
        pw = QtWidgets.QWidget(self)
        l = QtWidgets.QVBoxLayout(pw)
        pw.canvas = MyMplTrippleAxCanvas(pw, width=5, height=4, dpi=100)
        l.addWidget(pw.canvas)
        self.widgets['pw'] = pw

        self.data['plot_options'] = dict(reduce=True,
                                         scale='lin')

        # Create and add other widgets
        dftw = DataFilesTreeWidget()
        dsw = DataSelectorWidget()
        self.widgets['dftw'] = dftw
        self.widgets['dsw'] = dsw

        # Create and add zoom pan buttons
        self.zp_toolbar = NavigationToolbar(pw.canvas, self)
        self.zp_toolbar.hide()

        # Just some button
        self.zoom_button = QtWidgets.QPushButton('Zoom')
        self.zoom_button.clicked.connect(pw.canvas.zoom)

        self.pan_button = QtWidgets.QPushButton('Pan')
        self.pan_button.clicked.connect(pw.canvas.pan)

        self.home_button = QtWidgets.QPushButton('Home')
        self.home_button.clicked.connect(pw.canvas.home)

        # set the layout
        pzw = QtWidgets.QWidget(self)
        l2 = QtWidgets.QHBoxLayout(pzw)

        l2.addWidget(self.zoom_button)
        l2.addWidget(self.pan_button)
        l2.addWidget(self.home_button)

        self.widgets['pzw'] = pzw

        # Create splitter
        splitter1 = QtWidgets.QSplitter(Qt.Horizontal)
        splitter2 = QtWidgets.QSplitter(Qt.Vertical)
        splitter3 = QtWidgets.QSplitter(Qt.Vertical)

        splitter3.addWidget(pw)
        splitter3.addWidget(pzw)

        splitter2.addWidget(dftw)
        splitter2.addWidget(dsw)
        splitter2.setSizes([400, 400])

        splitter1.addWidget(splitter2)
        splitter1.addWidget(splitter3)
        splitter1.setSizes([200, 600])

        # set the focus
        pw.setFocus()

        # add widgets to window
        self.setCentralWidget(splitter1)

        self.statusBar().showMessage("All hail matplotlib!", 2000)

    def fileOpen(self):
        self.data['lvdts'] = read_lvdt_data(statusbar = self.statusBar())
        self.data['pt100T'] = read_pt100_data(statusbar = self.statusBar())
        self.data['hoboT'] = read_hobo_data(statusbar = self.statusBar())

        history[:] = []
        hist = read_history(hist_fname)
        if hist is not None: history.extend(hist)

        self.widgets['dsw'].setupContent()
        self.statusBar().showMessage("Done...", 5000)


    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtWidgets.QMessageBox.about(self, "About",
                                    """embedding_in_qt5.py example
Copyright 2005 Florent Rougon, 2006 Darren Dale, 2015 Jens H Nielsen

This program is a simple example of a Qt5 application embedding matplotlib
canvases.

It may be used and modified with no restriction; raw copies as well as
modified versions may be distributed without limitation.

This is modified from the embedding in qt4 example to show the difference
between qt4 and qt5"""
                                )


def read_lvdt_data(statusbar=None):
    fname = 'CR6_SN4983_LVDT_table.dat'
    fname = os.path.join(data_path, fname)

    if statusbar is not None:
        statusbar.showMessage(fname, 2000)  # print to screen so we can see the progress

    # Read LVDT data
    lvdts = pd.read_csv(fname, sep=None, skiprows=4,
                        names=['datetime', 'recnum', 'Volt(1)', 'Volt(2)', 'Volt(3)'],
                        index_col=0, parse_dates=[0], dayfirst=False,
                        na_values=["NAN"])

    # Apply LVDT calibrations

    lvdts['Volt(1)'] = 10 - ((lvdts['Volt(1)'] - lvdt_1_b) / lvdt_1_a)  # - lvdt_1_h0
    lvdts['Volt(2)'] = 10 - ((lvdts['Volt(2)'] - lvdt_2_b) / lvdt_2_a)  # - lvdt_2_h0
    lvdts['Volt(3)'] = 10 - ((lvdts['Volt(3)'] - lvdt_3_b) / lvdt_3_a)  # - lvdt_2_h0

    return lvdts


def read_pt100_data(statusbar=None):
    fname = 'CR6_SN4983_PT100_table.dat'
    fname = os.path.join(data_path, fname)

    if statusbar is not None:
        statusbar.showMessage(fname, 2000)  # print to screen so we can see the progress

    # Read pt100 data
    pt100T = pd.read_csv(fname, sep=None, skiprows=4,
                         names=['datetime', 'recnum', 'T(1)', 'T(2)', 'X(1)', 'X(2)'], index_col=0, parse_dates=[0],
                         dayfirst=False)

    # Apply pt100 calibrations
    pt100T['T(1)'] = PT100_converter.calc_temperature(pt100T['X(1)'] * cal_pt100_1)
    pt100T['T(2)'] = PT100_converter.calc_temperature(pt100T['X(2)'] * cal_pt100_2)

    return pt100T


def read_hobo_data(statusbar=None):
    # Read hobo data
    # This is slightly complicated because data are in several files.
    # So we read each file and append the rows of each file to a list
    # When all data are read, we generate the dataframe...

    rows = []  # holds the combined rows of all files
    index = []  # holds the combined list of datetimes (one for each row)
    # utc_list = []  # holds the list of utc time off sets
    directory = os.path.join(data_path,'hobo') # This is where we read the datafiles from
    for filename in os.listdir(directory):  # Iterate over all files in the directory
        if filename.endswith(".csv"):  # If file is a csv file we should process it
            fname = os.path.join(directory, filename)
            if statusbar is not None:
                statusbar.showMessage(fname, 2000)  # print to screen so we can see the progress

            print(fname)
            # read the first two lines, so we can find the time zone information
            with open(fname, "r", encoding='latin-1') as f:
                line1 = f.readline()
                line2 = f.readline()
                line3 = f.readline()
            print('{0:30s}'.format(line3))

            # Extract the time zone infromation
            time_zone = line2[line2.index('GMT'): line2.index('","Temp')]
            utc_offset = float(time_zone[4:6]) * 3600 + float(time_zone[7:]) * 60  # in seconds
            if time_zone[3] == '-':
                utc_offset = -utc_offset
            utc_offset_dt = dt.timedelta(days=0, seconds=utc_offset)

            # Read the current file
            df = pd.read_csv(fname, sep=',', skiprows=2, index_col=0, usecols=[1, 2, 3, 4],
                             names=['datetime', 'T(1)', 'T(2)', 'T(3)'], parse_dates=[0], dayfirst=False)

            # add rows and datetimes to the respective lists
            rows.extend(df.values.tolist())
            index.extend((df.index + utc_offset_dt).tolist())
            # utc_list.extend([utc_offset_dt for i in xrange(len(df))])

    # apply the timezone offset
    # index = np.array(index) + np.array(utc_list)

    # Now create the final dataframe
    hoboT = pd.DataFrame(rows, index=index, columns=['T(1)', 'T(2)', 'T(3)'])

    # The rows are mixed up, because the files were not necessarily read
    # in the right order (according to time)
    # so now we sort all the rows timewise
    hoboT.sort_index(inplace=True)

    # apply the hobo calibrations
    hoboT['T(1)'] += cal_hobo_1
    hoboT['T(2)'] += cal_hobo_2
    hoboT['T(3)'] += cal_hobo_3

    return hoboT


def read_history(fname):
    fname = os.path.join(data_path, fname)
    with open(fname, 'r') as stream:
        try:
            dat = yaml.load(stream)
            hist = dat['history']
        except yaml.YAMLError as exc:
            hist = []
            print(exc)

    for ix,val in enumerate(hist):
        try:
            val['time1'] = dt.datetime.strptime(hist[ix]['time1'], '%H:%M:%S').time()
            val['time2'] = dt.datetime.strptime(hist[ix]['time2'], '%H:%M:%S').time()
            hist[ix] = val
        except:
            pass

    return hist


def plot_step(lvdts, pt100T, hoboT, step=0):
    id = step
    hist = history[step]
    if hist['date'] is None: return  # Only plot if the step has a date-value assigned

    # if not (hist['step'] == 15): continue

    # Get the start date and time
    step_starttime = dt.datetime.combine(hist['date'], hist['time1'])

    # Get the end date and time
    if id + 1 > len(history) or history[id + 1]['time1'] is None:
        # If this is the last step, use the last time in the dataset
        step_endtime = max(max(hoboT.index), max(pt100T.index))
    else:
        # If this is not the last step, use the start time for the next step minus 1 sec
        step_endtime = dt.datetime.combine(history[id + 1]['date'], history[id + 1]['time1'])
        step_endtime -= dt.timedelta(seconds=1)

    # select data for plotting: only data within the present load/temperature step start and end times
    lvdt_step = lvdts[(lvdts.index >= step_starttime) & (lvdts.index <= step_endtime)]
    hobo_step = hoboT[(hoboT.index >= step_starttime) & (hoboT.index <= step_endtime)]
    pt100_step = pt100T[(pt100T.index >= step_starttime) & (pt100T.index <= step_endtime)]

    f, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False,
                                           gridspec_kw={'height_ratios': [4, 4, 4, 1]},
                                           figsize=(10, 12))

    # make sure the top 3 axes share the same x-axis
    ax1.get_shared_x_axes().join(ax1, ax2, ax3)

    # Plot the first lvdt
    l1, = ax1.plot_date(lvdt_step.index, lvdt_step['Volt(1)'], '-r', label=lvdt_labels[0], zorder=10)

    # Create a second y-axis and plot the second lvdt
    ax1a = ax1.twinx()
    l2, = ax1a.plot_date(lvdt_step.index, lvdt_step['Volt(2)'], '-b', label=lvdt_labels[1], zorder=10)

    # invert the direction of both y-axes
    ax1.invert_yaxis()
    ax1a.invert_yaxis()

    l3, = ax2.plot_date(pt100_step.index, pt100_step['T(1)'], '-r', label=pt100_labels[0], zorder=10)
    l4, = ax2.plot_date(pt100_step.index, pt100_step['T(2)'], '-b', label=pt100_labels[1], zorder=10)
    l5, = ax2.plot_date(hobo_step.index, hobo_step['T(1)'], '-g', label=hobo_labels[0], zorder=2)
    l6, = ax2.plot_date(hobo_step.index, hobo_step['T(2)'], '-y', label=hobo_labels[1], zorder=1)

    l7, = ax3.plot_date(hobo_step.index, hobo_step['T(3)'], '-k', label=hobo_labels[2])

    # Plot figure title
    plt.suptitle('Step: {0},  Date: {1},  Load: {2}, kPa  Temp: {3} C'.format(hist['step'],
                                                                              hist['date'],
                                                                              hist['load'],
                                                                              hist['temp']),
                 fontsize=14)

    ax1.set_ylabel('Strain, Sample 1 [mm]')
    ax1a.set_ylabel('Strain, Sample 2 [mm]')
    ax2.set_ylabel('Temperature [C]')
    ax3.set_ylabel('Temperature [C]')

    # plot common legend for all subplots
    handles = [l3, l4, l5, l6, l7]
    labels = [h.get_label() for h in handles]

    ax4.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fontsize=12)
    ax4.axis('off')

    # f.tight_layout()
    # ax1.legend(zorder=0)
    # ax2.legend(zorder=0)




qApp = QtWidgets.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
