
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
from scipy.optimize import fsolve

# Make sure that we are using QT5
import yaml

import sqrtlogscale

matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from sqrtlogscale import SqrtLogScale

progname = os.path.basename(sys.argv[0])
progversion = "0.1"


# ----------------------------------------------------------------------------------------------------------------------
# Data definitions
# ----------------------------------------------------------------------------------------------------------------------

# TEST INFORMATION

history = []
hist_fname = 'hist.yml'
data_path = './sample_1_2'

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

# HOBO CALIBRATIONS:   T = T_meas + delta_T
cal_hobo_1 = -0.02
cal_hobo_2 = -0.02
cal_hobo_3 = 0.00



# PT100 calibrations:  multiplicative factor
cal_pt100_1 = 1.000320692
cal_pt100_2 = 1.001201753



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


                        
def get_closest(a, vals):
    """Find the element of the first array that is closest to each element of
    the second array.
    
    Parameters
    ----------
    a: array
        First input array
    vals: array
        Second input array, nominal values to be used as basis for comparison
    
    Returns
    -------
    idx: array of ints
        Resulting array of indices into 'a'.
    """
    
    # make sure array is a numpy array
    array = np.array(a)

    # get insert positions
    idx = np.searchsorted(a, vals, side="left")

    # find indexes where previous index is closer
    prev_idx_is_less = ((idx == len(a)) | (np.fabs(vals - array[np.maximum(idx - 1, 0)]) < np.fabs(
        vals - array[np.minimum(idx, len(a) - 1)])))

    idx[prev_idx_is_less] -= 1

    return idx

    
def log_sampling(a, num=20):
    """Samples the passed array of values (typically times measured in seconds
    or minutes) by a certain number of samples per log cycle (default: 
    20 samples per decade).
    
    The function produces an array of log-spaced values from the first to the
    last value in the input array (from second value if the first is zero)
    using 'num' samples per decade. It then finds the values of the input array
    that are closest to the values in the log-spaced array, and returns the
    indices to these values.
    
    Typical use-case is to resample a regularly (linearly) sampled dataset to 
    logarithmic spacing for visualization purposes.
    
    Parameters
    ----------
    a: array-like
        Input array of values/times to sample by certain number of samples
        pr log-cycle
    num: int
        Number of samples to use per log-cycle
    
    Returns
    -------
    idx: array of ints
        Array of indices into the input array
    """
    
    log_t0 = np.log10(a[1])
    log_tmax = np.log10(a[-1])
    nomvals = 10 ** np.linspace(log_t0, log_tmax, np.floor((log_tmax - log_t0) * num))
    idx = get_closest(a, nomvals)
    return np.unique(idx)

def sqrt_sampling(a, num=20):
    """Samples the passed array of values (typically times measured in seconds
    or minutes) equidistantly on a sqrt scale.
    
    The function produces an array of log-spaced values from the first to the
    last value in the input array (from second value if the first is zero)
    using 'num' samples per decade. It then finds the values of the input array
    that are closest to the values in the log-spaced array, and returns the
    indices to these values.
    
    Typical use-case is to reduce a regularly (linearly) sampled dataset to 
    logarithmic spacing for visualization purposes.
    
    Parameters
    ----------
    a: array-like
        Input array of values/times to reduce to certain numper of samples
        pr log-cycle
    num: int
        Number of samples to use per log-cycle
    
    Returns
    -------
    idx: array of ints
        Array of indices into the input array
    """
    
    sqrt_t0 = np.sqrt(a[0])
    sqrt_tmax = np.sqrt(a[-1])
    nomvals = np.linspace(sqrt_t0, sqrt_tmax, np.floor((sqrt_tmax - sqrt_t0) * num)) ** 2
    idx = get_closest(a, nomvals)
    return np.unique(idx)
    

def read_lvdt_data(fname, channel, calibration, lvdt_h0, sample_h0):
    fname = os.path.join(data_path, fname)

    # Read LVDT data
    lvdt_dat = pd.read_csv(fname, sep=None, skiprows=4,
                           names=['datetime', 'recnum', 'Volt'],
                           usecols=[0,1,2+channel],
                           index_col=0, parse_dates=[0], dayfirst=False,
                           na_values=["NAN"],
                           engine='python')

    # Apply LVDT calibration
    lvdt_dat['strain'] = ((lvdt_dat['Volt'] - calibration['coeff'][1]) / calibration['coeff'][0])
    lvdt_dat['eps'] = ((lvdt_dat['strain'] - lvdt_h0) / sample_h0) * 100        
    if calibration['negate']: lvdt_dat['eps'] = lvdt_dat['eps']*-1 
    
    lvdt_dat.drop(labels=['recnum'], axis=1, inplace=True)
    
    return lvdt_dat


def read_pt100_data(fname, channel, calibration):
    fname = os.path.join(data_path, fname)

    # Read pt100 data
    pt100T = pd.read_csv(fname, sep=None, skiprows=[0,2,3], index_col=0, parse_dates=[0],
                         dayfirst=False, engine='python')

    nchannels = int((len(pt100T.columns)-1)/2)
    
    # Apply pt100 calibrations
    out = pd.DataFrame(PT100_converter.calc_temperature(calibration['coeff'][0] * 
                                                        pt100T.iloc[:,1+nchannels+channel] + 
                                                        calibration['coeff'][1]),
                       index = pt100T.index, columns=['T'])

    return out


def read_hobo_data(hobo_conf):
    # Read hobo data
    directory = os.path.join(data_path,'hobo') # This is where we read the datafiles from
    hdf_fname = os.path.join(directory,'hobo_dat.hdf')
    
    if os.path.exists(hdf_fname):
        hoboT = pd.read_hdf(hdf_fname, 'data')
    else:
        # This is slightly complicated because data are in several files.
        # So we read each file and append the rows of each file to a list
        # When all data are read, we generate the dataframe...

        rows = []  # holds the combined rows of all files
        index = []  # holds the combined list of datetimes (one for each row)
        # utc_list = []  # holds the list of utc time off sets

        # global offset in seconds to be added to the time stamp
        global_offset = dt.timedelta(days=0, seconds=hobo_conf.get('deltatime', 0))
        
        for filename in os.listdir(directory):  # Iterate over all files in the directory
            if filename.endswith(".csv"):  # If file is a csv file we should process it
                fname = os.path.join(directory, filename)

                print('\n\n\n\nFile: {0}'.format(fname))
                # read the first two lines, so we can find the time zone information
                with open(fname, "r", encoding='latin-1') as f:
                    line1 = f.readline()
                    line2 = f.readline()
                    line3 = f.readline()
                print('{0:30s}'.format(line3))

                # Extract the time zone infromation
                time_zone = line2[line2.index('GMT'): line2.index('","Temp')]
                utc_offset = (float(time_zone[4:6]) * 3600 + float(time_zone[7:]) * 60)  # in seconds
                if time_zone[3] == '-':
                    utc_offset = -utc_offset
                utc_offset_dt = dt.timedelta(days=0, seconds=utc_offset)
                print('tz: {0}, tz offset: {1}, tz offset dt: {2}'.format(time_zone, utc_offset, utc_offset_dt))
                
                # Read the current file
                df = pd.read_csv(fname, sep=',', skiprows=2, index_col=0, usecols=[1, 2, 3, 4],
                                 names=['datetime', 'T(1)', 'T(2)', 'T(3)'], parse_dates=[0], dayfirst=False)

                id = len(index)
                # add rows and datetimes to the respective lists
                rows.extend(df.values.tolist())
                index.extend((df.index - utc_offset_dt + global_offset).tolist())
                print(index[id])
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
        hoboT['T(1)'] = hobo_conf['config'][0]['calibration']['coeff'][0]*hoboT['T(1)']+hobo_conf['config'][0]['calibration']['coeff'][1]
        hoboT['T(2)'] = hobo_conf['config'][1]['calibration']['coeff'][0]*hoboT['T(2)']+hobo_conf['config'][1]['calibration']['coeff'][1]
        hoboT['T(3)'] = hobo_conf['config'][2]['calibration']['coeff'][0]*hoboT['T(3)']+hobo_conf['config'][2]['calibration']['coeff'][1]

        hoboT.to_hdf(hdf_fname, 'data')
    
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

    for ix,val in enumerate(dat['history']):
        try:
            val['time1'] = dt.datetime.strptime(dat['history'][ix]['time1'], '%H:%M:%S').time()
            val['time2'] = dt.datetime.strptime(dat['history'][ix]['time2'], '%H:%M:%S').time()
            dat['history'][ix] = val
        except:
            pass

    return dat


def select_data(dfs, hist, buffer=0, end_buffer=0, max_len=0):
    """Returns the subset of the DataFrame(s) data that fall within
    the start and end date and time of the specified load step.
    
    If an end date or time is not specified, the function will use 
    the start date and time of the following load step less one second,
    or if this is the last load step, the last timestamp of the dataset.
    
    if a list of DataFrames is passed, all DataFrames will be processed
    and a list of DataFrames returned.
    
    Parameters
    ----------
    dfs: DataFrame or list of DataFrames
        The data to subset from. Must have timestamp index.
    hist: dict
        Dictionary conataining the keys 'date', 'time1', 'date2', 'time2'. 
        Dates must be datetime.date objects and times datetime.time objects.
    buffer: int
        Time in seconds to include before start of load step
    end_buffer: int
        Time in seconds to include after end of load step
    max_len: int
        Maximum length of selected time series in seconds. If 0 (default)
        no maximum length is imposed.
    
    Returns:
    --------
    dfs: DataFrame or list of DataFrames
        The subset of data from the input DataFrame(s)
    """
    
    if hasattr(dfs, 'extend'):
        dfs_islist = True
    else:
        dfs_islist = False
        
    
    if hist['date'] is None or hist['time1'] is None:
        raise ValueError('The load step must have at least the start date and time defined')       # Only plot if the step has a date-value assigned

    # Get the start date and time
    step_starttime = dt.datetime.combine(hist['date'], hist['time1']) - dt.timedelta(seconds=buffer)
    
    # Get the end date and time
    if ('date2' in hist and hist['date2'] is not None) and ('time2' in hist and hist['time2'] is not None):
        step_endtime = dt.datetime.combine(hist['date2'], hist['time2']) + dt.timedelta(seconds=end_buffer)
    else:
        if dfs_islist:
            step_endtime = max([max(dfi) for dfi in [df.index for df in dfs]])
        else:
            step_endtime = max(dfs.index)
    
    if max_len > 0:
        series_length = (step_endtime-step_starttime).days*3600*24+(step_endtime-step_starttime).seconds
        if series_length > max_len:
            step_endtime = step_starttime + dt.timedelta(seconds=max_len)
    
    # select data for plotting: only data within the present load/temperature step start and end times
    if dfs_islist:
        dfs = [df[(df.index >= step_starttime) & (df.index <= step_endtime)] for df in dfs]
    else:
        dfs = dfs[(dfs.index >= step_starttime) & (dfs.index <= step_endtime)]
        
    return dfs    

    
def plot_step_overview_hobo(lvdt_step, pt100_step, hobo_step, step_info):

    hobo_labels = ['Isolated chamber', 'Frost chamber', 'Computer room'] 
    
    
    f, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False,
                                           gridspec_kw={'height_ratios': [4, 4, 4, 1]},
                                           figsize=(10, 12))

    # make sure the top 3 axes share the same x-axis
    ax1.get_shared_x_axes().join(ax1, ax2, ax3)

    # Plot the first lvdt
    l1, = ax1.plot_date(lvdt_step.index, lvdt_step['strain'], '-k', label='LVDT strain', zorder=10)

    # # Create a second y-axis and plot the second lvdt
    # ax1a = ax1.twinx()
    # l2, = ax1a.plot_date(lvdt_step.index, lvdt_step['eps'], '-k', label='Strain', zorder=10)

    # invert the direction of both y-axes
    ax1.invert_yaxis()

    l3, = ax2.plot_date(pt100_step.index, pt100_step['T'], '-b', label='Sample temperature', zorder=10)
    handles = [l1, l3]
    
    if len(hobo_step) > 0:
        #l5, = ax2.plot_date(hobo_step.index, hobo_step['T(1)'], ls='-', c='g', marker=None, label=hobo_labels[0], zorder=2)
        #l6, = ax2.plot_date(hobo_step.index, hobo_step['T(2)'], ls='-', c='orange', marker=None, label=hobo_labels[1], zorder=1)
        l7, = ax3.plot_date(hobo_step.index, hobo_step['T(3)'], ls='-', c='r', marker=None, label=hobo_labels[2])
        #handles.extend([l5, l6, l7])
        handles.extend([l7])
        
    # Plot figure title
    plt.suptitle('Step: {0},  Date: {1},  Load: {2}, kPa  Temp: {3} C'.format(step_info['step'],
                                                                              step_info['date'],
                                                                              step_info['load'],
                                                                              step_info['temp']),
                 fontsize=14)

    ax1.set_ylabel('Strain [mm]')
    #ax1a.set_ylabel('Strain [%]')
    ax2.set_ylabel('Temperature [C]')
    ax3.set_ylabel('Temperature [C]')

    # plot common legend for all subplots
    labels = [h.get_label() for h in handles]

    ax4.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fontsize=12)
    ax4.axis('off')

    axes = [ax1, ax2, ax3]
    for ax in axes:
        ax.fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.grid(True)
    # f.tight_layout()
    # ax1.legend(zorder=0)
    # ax2.legend(zorder=0)
    
    return f

def plot_step_overview_hobo2(lvdt_dat, pt100_dat, hobo_dat, step_info):

    # select the data
    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100T, hoboT], hist, buffer=60, max_len=120)
    step_starttime = dt.datetime.combine(step_info['date'], step_info['time1'])
    step_endtime = dt.datetime.combine(step_info['date2'], step_info['time2'])
    lvdt_start = add_minutes(lvdt_step, t0=step_starttime)
    pt100_start = add_minutes(pt100_step, t0=step_starttime)
    hobo_start = add_minutes(hobo_step, t0=step_starttime)

    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100T, hoboT], hist, buffer=0, end_buffer=0, max_len=0)
    lvdt_step = add_minutes(lvdt_step)
    pt100_step = add_minutes(pt100_step)
    hobo_step = add_minutes(hobo_step)
    
    hobo_labels = ['Isolated chamber', 'Frost chamber', 'Computer room'] 
    
    #f, axs = plt.subplots(nrows=4, ncols=2, sharex=False, sharey=False,
    #                      gridspec_kw={'height_ratios': [4, 4, 4, 1],
    #                                   'width_ratios':  [1, 3]},
    #                      figsize=(10, 12))

    f = plt.figure(figsize=(10,9))
    gs = matplotlib.gridspec.GridSpec(12, 4)
    plt.subplot(gs.new_subplotspec((1, 0), colspan=1, rowspan=4))
    plt.subplot(gs.new_subplotspec((1, 1), colspan=3, rowspan=4))
    plt.subplot(gs.new_subplotspec((5, 1), colspan=3, rowspan=3))
    plt.subplot(gs.new_subplotspec((8, 1), colspan=3, rowspan=3))
    plt.subplot(gs.new_subplotspec((11, 1), colspan=3, rowspan=1))
    
    axs = f.axes
    
    # make sure the top 3 axes share the same x-axis
    axs[1].get_shared_x_axes().join(axs[1],axs[2], axs[3])
    
    # Plot the first lvdt
    l, = axs[0].plot(lvdt_start['minutes'], lvdt_start['eps'], '-', c='0.3', lw=0.5, marker='.', ms=4, mec='k', mfc='k', zorder=10)
    l1, = axs[1].plot_date(lvdt_step.index, lvdt_step['eps'], '-k', label='LVDT strain', zorder=10)
    
    axs[0].axvline(x=0, ls='--', color='0.65', zorder=-100)
    axs[1].axvline(x=step_starttime, ls='--', color='0.65', zorder=-100)
    axs[2].axvline(x=step_starttime, ls='--', color='0.65', zorder=-100)
    axs[3].axvline(x=step_starttime, ls='--', color='0.65', zorder=-100)
    axs[1].axvline(x=step_endtime, ls='--', color='0.65', zorder=-100)
    axs[2].axvline(x=step_endtime, ls='--', color='0.65', zorder=-100)
    axs[3].axvline(x=step_endtime, ls='--', color='0.65', zorder=-100)
    
    y0 = lvdt_step[lvdt_step['minutes']==0]['eps'].values[0]
    axs[0].axhline(y=y0, ls='--', color='0.65', zorder=-100)
    axs[1].axhline(y=y0, ls='--', color='0.65', zorder=-100)
    
    # invert the direction of both y-axes
    axs[1].invert_yaxis()
    axs[0].invert_yaxis()
    
    l3, = axs[2].plot_date(pt100_step.index, pt100_step['T'], '-b', label='Sample temperature', zorder=10)
    handles = [l1, l3]
    
    if len(hobo_step) > 0:
        l5, = axs[2].plot_date(hobo_step.index, hobo_step['T(1)'], ls='-', c='orange', marker=None, label=hobo_labels[0], zorder=2)
        #l6, = ax2.plot_date(hobo_step.index, hobo_step['T(2)'], ls='-', c='orange', marker=None, label=hobo_labels[1], zorder=1)
        l7, = axs[3].plot_date(hobo_step.index, hobo_step['T(3)'], ls='-', c='r', marker=None, label=hobo_labels[2])
        #handles.extend([l5, l6, l7])
        handles.extend([l5, l7])
    else:
        axs[3].plot_date(lvdt_step.index[0], lvdt_step['strain'][0], '-', color='none', zorder=10)
        bbox_props = dict(boxstyle="square,pad=0.3", fc="none", ec="none")
        t = axs[3].text(0.5, 0.5, "No data available", 
                          ha="center", va="center", rotation=0,
                          size=15, color='0.75', bbox=bbox_props,
                          transform=axs[3].transAxes)
        
        
    # Plot figure title
    plt.suptitle('Step: {0},  Date: {1},  Load: {2}, kPa  Temp: {3} C'.format(step_info['step'],
                                                                              step_info['date'],
                                                                              step_info['load'],
                                                                              step_info['temp']),
                 fontsize=14)

    axs[0].set_xlabel('Time [min]')
    axs[0].set_ylabel('Strain [%]')
    axs[1].set_ylabel('Strain [%]')
    #ax1a.set_ylabel('Strain [%]')
    axs[2].set_ylabel('Temperature [C]')
    axs[3].set_ylabel('Temperature [C]')

    # plot common legend for all subplots
    labels = [h.get_label() for h in handles]

    axs[4].axis('off')
    axs[4].legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fontsize=10)
    
    axes = [axs[1], axs[2]]
    if len(hobo_step) > 0:
        axes.append(axs[3])
        
    for ax in axes:
        ax.fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.grid(True)
    # f.tight_layout()
    # ax1.legend(zorder=0)
    # ax2.legend(zorder=0)
    
    def ignore_spikes(ax, buffer, percentiles=[1,5], threshold=1, min_lim=None):
        ylim = ax.get_ylim()
        if ylim[0]>ylim[1]:
            inverted = True
        else:
            inverted = False
            
        yln = min(ylim)
        ylx = max(ylim)
        
        dat = []
        for l in ax.lines:
            dat.extend(list(l.get_data()[1]))
        
        if np.abs(np.nanpercentile(dat,100-percentiles[0])-np.nanpercentile(dat,100-percentiles[1])) < threshold:
            ylx = np.nanpercentile(dat,100-percentiles[0])+buffer
            
        if np.abs(np.nanpercentile(dat,percentiles[1])-np.nanpercentile(dat,percentiles[0])) < threshold:
            yln = np.nanpercentile(dat,percentiles[0])-buffer
        
        if min_lim is not None:
            if ylx < min_lim[1]:
                ylx = min_lim[1]
            if yln > min_lim[0]:
                yln = min_lim[0]
        
        if inverted:
            ax.set_ylim([ylx, yln])
        else:
            ax.set_ylim([yln, ylx])

    def min_ylim(ax, min_lim):
        ylim = ax.get_ylim()
        if ylim[0]>ylim[1]:
            inverted = True
        else:
            inverted = False

        yln = min(ylim)
        ylx = max(ylim)
    
        if ylx-yln < min_lim:
            center = (yln+ylx)/2.
            yln = center-min_lim/2.
            ylx = center+min_lim/2.
        
        if inverted:
            ax.set_ylim([ylx, yln])
        else:
            ax.set_ylim([yln, ylx])
            
    min_lims = [0.1, 1, 0.4, 2]
    
    min_ylim(axs[0], min_lims[0])
    
    #if step_info['step'] == 27:
    #    pdb.set_trace()
    
    min_lim = np.sort([lvdt_step['eps'][0], lvdt_step['eps'][-1]])
    min_lim[0] -= 0.2
    min_lim[1] += 0.2
    
    ignore_spikes(axs[1], min_lims[1]/2., min_lim=min_lim)
    min_ylim(axs[1], min_lims[1])
    
    ignore_spikes(axs[2], min_lims[2]/2.)
    min_ylim(axs[2], min_lims[2])            
    
    min_ylim(axs[3], min_lims[3])            
    
    f.tight_layout()
    
    return f    
    

    
    
def plot_step_overview(lvdt_step, pt100_step, hobo_step, step_info):

    f, (ax1, ax2, ax4) = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False,
                                           gridspec_kw={'height_ratios': [4, 4, 1]},
                                           figsize=(10, 10))

    # make sure the top 2 axes share the same x-axis
    ax1.get_shared_x_axes().join(ax1, ax2)

    # Plot the first lvdt
    l1, = ax1.plot_date(lvdt_step.index, lvdt_step['strain'], '-k', label='LVDT strain', zorder=10)

    # invert the direction of the y-axes
    ax1.invert_yaxis()

    l3, = ax2.plot_date(pt100_step.index, pt100_step['T'], '-b', label='Sample temperature', zorder=10)
    handles = [l1, l3]
    
    # Plot figure title
    plt.suptitle('Step: {0},  Date: {1},  Load: {2}, kPa  Temp: {3} C'.format(step_info['step'],
                                                                              step_info['date'],
                                                                              step_info['load'],
                                                                              step_info['temp']),
                 fontsize=14)

    ax1.set_ylabel('Strain [mm]')
    ax2.set_ylabel('Temperature [C]')

    # plot common legend for all subplots
    labels = [h.get_label() for h in handles]

    ax4.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fontsize=12)
    ax4.axis('off')

    axes = [ax1, ax2]
    for ax in axes:
        ax.fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.grid(True)
    # f.tight_layout()
    
    return f
    
    
    
def plot_time_curve(strain, step_info, temp=None, intersect=1, markers=False):
    
    strain = add_minutes(strain)
    
    if temp is not None:
        f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False,
                                               gridspec_kw={'height_ratios': [3, 1]},
                                               figsize=(10, 5))
        temp = add_minutes(temp)
    else:
        f = plt.figure(figsize=(15,5))
        ax1 = plt.axis()
        
    if markers:
        ls = '.-k'
    else:
        ls = '-k'
    
    ax1.plot(strain['minutes'][1:], strain['eps'][1:], ls)
       
    if intersect is None:
        ax1.set_xscale('log')
    else:
        ax1.set_xscale('sqrtlog', intersectx=intersect)

        #ax.set_ylim([-2, np.ceil(np.sqrt(intersect))+1])
    #plt.gca().set_xlim([0, 2300])
    ax1.set_ylabel('Strain [%]')
    ax1.grid(True)
    ax1.grid(b=True, which='minor', ls='--')

    xlim = ax1.get_xlim()
    new_xlim = [xl for xl in xlim]
    base = np.power(10, np.floor(np.log10(strain['minutes'][-1])))
    new_xlim[1] = np.ceil(strain['minutes'][-1]/base)*base
    ax1.set_xlim(new_xlim)
    
    ax1.invert_yaxis()
    sqrtlogscale.annotate_xaxis(ax=ax1, intersect=intersect)
    
    if temp is not None:       
        if intersect is None:
            ax2.set_xscale('log')
            temp['minutes'].iloc[0] = ax1.get_xlim()[0]
        else:
            ax2.set_xscale('sqrtlog', intersectx=intersect)

        ax2.plot(temp['minutes'], temp['T'], '-k')

        ax2.set_ylabel('Temperature [C]')
        ax2.grid(True)
        ax2.grid(b=True, which='minor', ls='--')
    
        yl2 = ax2.get_ylim()
        if np.diff(yl2) < 0.4:
            midpoint = np.round(np.mean(yl2))
            ax2.set_ylim([midpoint-0.2, midpoint+0.2])
        
        sqrtlogscale.annotate_xaxis(ax=ax2, intersect=intersect, arrows=None)
        
        ax2.set_xlim(ax1.get_xlim())
        
        # make sure the axes share the same x-axis
        ax1.get_shared_x_axes().join(ax2)
        ax2.set_xlabel('Time (min)')
    else:
        ax1.set_xlabel('Time (min)')
    
    
    
    plt.show(block=False)
    return ax1
    
def add_minutes(ts, t0=None):
    """Adds timestamp in minutes with the first sample as zero minutes"""
#    if not 'minutes' in ts.columns:
    if len(ts) > 0:
        if t0 is None:
            ts['minutes'] = (ts.index-ts.index[0]).days*24*60 + (ts.index-ts.index[0]).seconds/60.
        else:
            ts['minutes'] = (ts.index-t0).days*24*60 + (ts.index-t0).seconds/60.
    
    return ts
    

def plot_fit_line(a, b, ax, num=10, ls='--r', type='linear', intersect=None):
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    if type == 'sqrt':
        func = np.sqrt
    elif type == 'log10':
        func = np.log10
    elif type == 'linear':
        func = lambda x: x    # Do nothing

        
    t0 = func(xlim[0])
    tmax = func(xlim[1])

    if type == 'sqrt':
        if intersect is not None:
            tmax = func(intersect*1.1)
        xvals = np.linspace(t0, tmax, int(np.floor((tmax - t0) * num))) ** 2
    elif type == 'log10':
        if np.isinf(t0): t0 = -6
        if intersect is not None:
            t0 = func(intersect*0.9)
        xvals = np.logspace(t0, tmax, int(np.floor((tmax - t0) * num)))
    elif type == 'linear':
        xvals = np.linspace(t0, tmax, int(np.floor((tmax - t0) * num)))
    
    y = a*func(xvals)+b
    ax.plot(xvals, y, ls)
    
    new_ylim = [yl for yl in ylim]
    
    y0 = np.floor(b*10)/10
    if y0 < new_ylim[1]:
        new_ylim[1] = y0
    
    ax.set_ylim(new_ylim)
    return xvals, y
    

def fit_line(times, values, t1, t2, use_points='endpoints', type='linear'):
    
    if type == 'sqrt':
        func = np.sqrt
    elif type == 'log10':
        func = np.log10
    elif type == 'linear':
        func = lambda x: x    # Do nothing
        
    if use_points == 'endpoints':
        id1 = np.argmin(np.abs(times-t1))
        id2 = np.argmin(np.abs(times-t2))
        
        if id1==id2:
            raise ValueError('Requested time range is not within the range of timestamps passed.')
        
        p = np.polyfit(func([times[id1], times[id2]]), [values[id1], values[id2]], 1)
        
    elif use_points == 'subsample':
        raise NotImplementedError()
    else:
        raise NotImplementedError()

    return p[0], p[1]
    

def find_t100(sqrt_a, sqrt_b, log_a, log_b):

    def func(t):
        return sqrt_a * np.sqrt(t) + sqrt_b - (log_a * np.log10(t) + log_b)
        
    t100 = fsolve(func, 0.4)
    
    return t100[0]
    

def sign_changes(a):
    """
    https://stackoverflow.com/a/2652425/1760389
    """
    asign = np.sign(a)
    sz = asign == 0
    while sz.any():
        asign[sz] = np.roll(asign, 1)[sz]
        sz = asign == 0
    signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
    signchange[0] = 0
    
    return signchange

def intersect_lines(p1, p2, p3, p4):
    x = ((p2[0]*p1[1]-p1[0]*p2[1])*(p4[0]-p3[0]) - (p4[0]*p3[1]-p3[0]*p4[1])*(p2[0]-p1[0])) / ((p2[0]-p1[0])*(p4[1]-p3[1]) - (p4[0]-p3[0])*(p2[1]-p1[1]))
    y = ((p2[0]*p1[1]-p1[0]*p2[1])*(p4[1]-p3[1]) - (p4[0]*p3[1]-p3[0]*p4[1])*(p2[1]-p1[1])) / ((p2[0]-p1[0])*(p4[1]-p3[1]) - (p4[0]-p3[0])*(p2[1]-p1[1]))
    return (x, y)    
    
    
def interpret_BrinchHansen(times, values, h0, t1, t2, t3, t4):    
    params = {}
    sqrt_a, sqrt_b = fit_line(times, values, t1, t2, use_points='endpoints', type='sqrt')
    log_a, log_b = fit_line(times, values, t3, t4, use_points='endpoints', type='log10')
    
    params['sqrt_a'] = sqrt_a
    params['sqrt_b'] = sqrt_b
    params['log_a'] = log_a
    params['log_b'] = log_b
    params['eps0'] = sqrt_b
    params['epss'] = log_a
    
    
    # Find (t100, eps100)
    params['t100'] = find_t100(sqrt_a, sqrt_b, log_a, log_b)
    idx = sign_changes(values-params['t100']).nonzero()[0]
    params['eps100'] = (values[idx[0]]-values[idx[0]-1])/(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1])) * (params['t100']-times[idx[0]-1]) + values[idx[0]-1]
    
    # Find (t50, eps50)
    params['eps50'] =  (params['eps100']+params['eps0'])/2.
    idx = sign_changes(values-params['eps50']).nonzero()[0]
    params['t50'] = ((params['eps50']-values[idx[0]-1])*(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1]))/(values[idx[0]]-values[idx[0]-1]) + np.sqrt(times[idx[0]-1]))**2.    
    
    params['Cv'] = np.NaN
    
    return params
    
    
def interpret_iso17892_5(times, values, h0, t1, t2, t3, t4):
    params = {}
    sqrt_a, sqrt_b = fit_line(times, values, t1, t2, use_points='endpoints', type='sqrt')
    
    sqrt115_vals = sqrt_a/1.15 * np.sqrt(times) + sqrt_b
    idx = sign_changes(sqrt115_vals-values).nonzero()[0]
    (x, y) = intersect_lines((np.sqrt(times[idx[-1]-1]), values[idx[-1]-1]),
                                   (np.sqrt(times[idx[-1]]), values[idx[-1]]),
                                   (np.sqrt(times[idx[-1]-1]), sqrt115_vals[idx[-1]-1]),
                                   (np.sqrt(times[idx[-1]]), sqrt115_vals[idx[-1]]))
    
    params['sqrt_a'] = sqrt_a
    params['sqrt_b'] = sqrt_b
    params['eps0'] = sqrt_b
    params['t90'] = x**2
    params['eps90'] = y
    params['eps100'] = (params['eps90']- params['eps0'])/0.9+params['eps0']
    params['eps50'] =  (params['eps100']+params['eps0'])/2.
    

#    x = (y-y1)*(x2-x1)/(y2-y1)+x1    
    idx = sign_changes(values-params['eps50']).nonzero()[0]
    params['t50'] = ((params['eps50']-values[idx[0]-1])*(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1]))/(values[idx[0]]-values[idx[0]-1]) + np.sqrt(times[idx[0]-1]))**2.
    
    idx = sign_changes(values-params['eps100']).nonzero()[0]
    params['t100'] = ((params['eps100']-values[idx[0]-1])*(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1]))/(values[idx[0]]-values[idx[0]-1]) + np.sqrt(times[idx[0]-1]))**2.

    
    Cv = 0.848 * h0**2 / params['t90']
    params['Cv'] = Cv
    
    log_a, log_b = fit_line(times, values, t3, t4, use_points='endpoints', type='log10')
    
    params['log_a'] = log_a
    params['log_b'] = log_b
    params['epss'] = log_a
   
    return params
    
    
def annotate_yaxis(ax, params, linecolor='black', linewidth=1, fontsize=10):
    
    xloc = 1.02
    
    vals = []
    keys = ['eps0','eps50','eps90','eps100']
    nom_labels = [r'$\varepsilon_0$',r'$\varepsilon_{50}$',r'$\varepsilon_{90}$',r'$\varepsilon_{100}$']
    labels = []
    for id, key in enumerate(keys):
        if key in params:
            vals.append(params[key])
            labels.append(nom_labels[id])
            
    #ax.annotate('', xy=(1.03, vals[0]), xytext=(1.03, vals[1]), xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
    #            arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth}, zorder=-100)
    for id in range(1,len(vals)):
        print('{0},   {1}'.format(vals[id-1], vals[id]))
        ax.annotate('', xy=(xloc, vals[id-1]), xytext=(xloc, vals[id]), xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
                arrowprops={'arrowstyle': '-', 'color':linecolor, 'linewidth':linewidth, 'shrinkA':0, 'shrinkB':0}, zorder=-100)
    
    for id, val in enumerate(vals):
        ax.annotate('', xy=(xloc-0.01, val), xytext=(xloc+0.01, val), xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
                arrowprops={'arrowstyle': '-', 'color':linecolor, 'linewidth':linewidth}, zorder=-100)
        
        ax.annotate(labels[id], xy=(xloc+0.01, val), xytext=(xloc+0.01, val), xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
                arrowprops={'arrowstyle': '-', 'color':'none', 'linewidth':linewidth}, zorder=-100, va='center')


def plot_legend(ax, info, step):
    
    txt = 'Sample: {0}\n'.format(info['sample']['name'])
    txt += 'Depth: {0}\n'.format(info['sample']['depth'])
    txt += 'Load: {0} kPa\n'.format(info['history'][step]['load'])
    txt += 'Temperature: {0} C'.format(info['history'][step]['temp'])
    
    at = AnchoredText(txt,
                  prop=dict(size=10), frameon=True,
                  loc='upper right')
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    
def get_key_index(lst, key, val):
    key_list = [d[key] for d in lst]
    try:
        return key_list.index(val)
    except: 
        return None

        
if __name__ == '__main__':
    info = read_history(hist_fname)
    history = info['history']
    interpret = info['interpret']
    
    
    lvdt_dat = read_lvdt_data(info['lvdt']['filename'],
                              info['lvdt']['channel'],
                              info['lvdt']['calibration'],
                              info['lvdt']['h0'],
                              info['sample']['h0'])
    
    pt100T = read_pt100_data(info['pt100']['filename'],
                             info['pt100']['channel'],
                             info['pt100']['calibration'])
    
    hoboT = read_hobo_data(info['hobo'])
    
    
    steps_to_interpret = [d['step'] for d in interpret]
    
    steps_to_plot = [] # [d['step'] for d in history]
    steps_to_overview = [d['step'] for d in history]
    
    for step_id, hist in enumerate(history):
        
        if hist['step'] in steps_to_overview or hist['step'] in steps_to_plot:
            # select the data
            [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100T, hoboT], hist)
            lvdt_step = add_minutes(lvdt_step)
            pt100_step = add_minutes(pt100_step)
            hobo_step = add_minutes(hobo_step)
        else:
            continue
        
        print('Plotting {0:.0f} kPa, {1:.0f} C'.format(hist['load'], hist['temp']))
        
        if hist['step'] in steps_to_overview:
            #f = plot_step_overview_hobo(lvdt_step, pt100_step, hobo_step, hist)
            f = plot_step_overview_hobo2(lvdt_dat, pt100T, hoboT, hist)
            ovname = '{0}_raw_{1:02.0f}_{2}kPa_{3}C.png'.format(info['sample']['name'].replace(' ','-'),
                                                  hist['step'], hist['load'], hist['temp'])
            f.savefig(os.path.join(data_path, ovname), dpi=200)
        
        if hist['step'] in steps_to_plot:
            params = {}
            if hist['step'] in steps_to_interpret:
                
                int_id = get_key_index(interpret, 'step', hist['step'])
                if int_id is not None:
                    # do interpretation
                    h0 = 20   # mm     What thickness should be used here?
                    if interpret[int_id]['type'] == 'iso_sqrt':
                        params = interpret_iso17892_5(lvdt_step['minutes'].values, lvdt_step['eps'].values, h0, 
                                                      interpret[int_id]['t1'], interpret[int_id]['t2'], 
                                                      interpret[int_id]['t3'], interpret[int_id]['t4'])
                    else:
                        raise ValueError('Unknown interpretation method!')
                
            if params:
                ax = plot_time_curve(lvdt_step, hist, temp=pt100_step, intersect=params['t100'], markers=False)
                plot_fit_line(params['sqrt_a'], params['sqrt_b'], ax, type='sqrt', num=100, intersect=params['t100'])
                plot_fit_line(params['sqrt_a']/1.15, params['sqrt_b'], ax, type='sqrt', num=100, ls='--r', intersect=params['t100'])
                plot_fit_line(params['log_a'], params['log_b'], ax, type='log10', num=100, ls='--r', intersect=params['t100'])
                annotate_yaxis(ax, params)
            else:
                ax = plot_time_curve(lvdt_step, hist, temp=pt100_step, intersect=None, markers=False)
                
            plot_legend(ax, info, step_id)
            
            plt.draw()
            
            tcname = '{0}_timecurve_{1:02.0f}_{2}kPa_{3}C.png'.format(info['sample']['name'].replace(' ','-'),
                                                  hist['step'], hist['load'], hist['temp'])
            f.savefig(os.path.join(data_path, tcname), dpi=200)
    
