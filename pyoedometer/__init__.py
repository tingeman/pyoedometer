import os
import numpy as np
import datetime as dt
import pandas as pd
import yaml
import matplotlib
from matplotlib import pyplot as plt
from . import sqrtlogscale
import pdb    
#
## Sample info
#sample_1_h0 = 20.50  # mm
#sample_2_h0 = 21.12  # mm
#
## LVDT calibrations
#lvdt_1_a = np.mean([-992.649, -1005.437])
#lvdt_1_b = np.mean([4339.163, 4404.368])
#lvdt_1_h0 = 9.0234  # mm
#
#lvdt_2_a = np.mean([-985.457, -1004.581])
#lvdt_2_b = np.mean([4125.023, 4227.782])
#lvdt_2_h0 = 9.0234  # mm
#
#lvdt_3_a = np.mean([-985.514, -1000.657])
#lvdt_3_b = np.mean([4816.420, 4888.629])
#lvdt_3_h0 = 0  # mm
#
## HOBO CALIBRATIONS:   T = T_meas + delta_T
#cal_hobo_1 = -0.02
#cal_hobo_2 = -0.02
#cal_hobo_3 = 0.00
#
#
#
## PT100 calibrations:  multiplicative factor
#cal_pt100_1 = 1.000320692
#cal_pt100_2 = 1.001201753



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
    

def read_lvdt_data(lvdt_info, sample_info):

    fname =       lvdt_info['filename']
    channel =     lvdt_info['channel']
    calibration = lvdt_info['calibration']
    lvdt_h0 =     lvdt_info['h0']
    sample_h0 =   sample_info['h0']

    # Read LVDT data
    lvdt_dat = pd.read_csv(fname, sep=None, skiprows=4,
                           names=['datetime', 'recnum', 'Volt'],
                           usecols=[0,1,2+channel],
                           index_col=0, parse_dates=[0], dayfirst=False,
                           na_values=["NAN"],
                           engine='python')

    # Apply LVDT calibration
    if calibration['type'].lower() == "linear":
        lvdt_dat['raw_strain'] = (calibration['coeff'][0]*lvdt_dat['Volt'] + calibration['coeff'][1])
    elif calibration['type'].lower() == "invlinear":
        lvdt_dat['raw_strain'] = ((lvdt_dat['Volt'] - calibration['coeff'][1]) / calibration['coeff'][0])
    else:
        raise ValueError("Calibration type specifier not recognized (should be 'linear' or 'invlinear')")
        
    lvdt_dat['strain'] = lvdt_dat['raw_strain'] - lvdt_h0
    if calibration['negate']: lvdt_dat['strain'] = lvdt_dat['strain']*-1 
    
    if 'corrections' in lvdt_info and lvdt_info['corrections'] is not None:
        for corr in lvdt_info['corrections']:
            
            t0 = dt.datetime.strptime(corr['t0'], '%Y-%m-%d %H:%M:%S')
            t1 = dt.datetime.strptime(corr['t1'], '%Y-%m-%d %H:%M:%S')
            idx = (lvdt_dat.index >= t0) & (lvdt_dat.index >= t0)
            lvdt_dat.loc[idx,'strain'] += corr['dh']
    
    lvdt_dat['eps'] = ((lvdt_dat['strain']) / sample_h0) * 100        
    
    
    lvdt_dat.drop(labels=['recnum'], axis=1, inplace=True)
    
    return lvdt_dat


def read_pt100_data(pt100_info, sample_info=None):
    fname = pt100_info['filename']
    channel = pt100_info['channel']
    calibration = pt100_info['calibration']

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
    directory = hobo_conf['datapath'] # This is where we read the datafiles from
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
        sep = hobo_conf.get('separator',',')
        
        for filename in os.listdir(directory):  # Iterate over all files in the directory
            if filename.endswith("."+hobo_conf.get('file_extension', 'csv')):  # If file is a csv file we should process it
                fname = os.path.join(directory, filename)

                print('\n\n\n\nFile: {0}'.format(fname))
                # read the first two lines, so we can find the time zone information
                with open(fname, "r", encoding='latin-1') as f:
                    line1 = f.readline()
                    line2 = f.readline()
                    line3 = f.readline()
                print('{0:30s}'.format(line3))

                # Extract the time zone infromation
                time_zone = line2[line2.index('GMT'): line2.index('"{0}"Temp'.format(sep))]
                utc_offset = (float(time_zone[4:6]) * 3600 + float(time_zone[7:]) * 60)  # in seconds
                if time_zone[3] == '-':
                    utc_offset = -utc_offset
                utc_offset_dt = dt.timedelta(days=0, seconds=utc_offset)
                print('tz: {0}, tz offset: {1}, tz offset dt: {2}'.format(time_zone, utc_offset, utc_offset_dt))
                
                # Read the current file
                df = pd.read_csv(fname, sep=sep, skiprows=2, index_col=0, usecols=[1, 2, 3, 4],
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


def read_config(fname):
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
            print('Problem with times in Step {0}'.format(val['step']))

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
        no maximum length is imposed. If negative, the maximum length will
        be counted from the end of the time series.
    
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
    if max_len < 0:
        series_length = (step_endtime-step_starttime).days*3600*24+(step_endtime-step_starttime).seconds
        if series_length > np.abs(max_len):
            step_starttime = step_endtime + dt.timedelta(days=0, seconds=max_len)
    
    # select data for plotting: only data within the present load/temperature step start and end times
    if dfs_islist:
        dfs = [df.loc[(df.index >= step_starttime).copy() & (df.index <= step_endtime),:] for df in dfs]
    else:
        dfs = dfs.loc[(dfs.index >= step_starttime).copy() & (dfs.index <= step_endtime),:]
        
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
        l5, = ax2.plot_date(hobo_step.index, hobo_step['T(1)'], ls='-', c='g', marker=None, label=hobo_labels[0], zorder=2)
        l6, = ax2.plot_date(hobo_step.index, hobo_step['T(2)'], ls='-', c='orange', marker=None, label=hobo_labels[1], zorder=1)
        l7, = ax3.plot_date(hobo_step.index, hobo_step['T(3)'], ls='-', c='r', marker=None, label=hobo_labels[2])
        handles.extend([l5, l6, l7])
        #handles.extend([l7])
        
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
    
    
def plot_step_overview_hobo2(lvdt_dat, pt100_dat, hobo_dat, step_info, plot_mV=False):

    plot_end_of_step = True
    
    if plot_mV:
        lvdt_col = 'Volt'
        lvdt_unit = 'mV'
    else:
        lvdt_col = 'eps'
        lvdt_unit = '%'
    
    step_starttime = dt.datetime.combine(step_info['date'], step_info['time1'])
    step_endtime = dt.datetime.combine(step_info['date2'], step_info['time2'])
    step_duration_h = (step_endtime-step_starttime).days*24. + (step_endtime-step_starttime).seconds/3600.    
    
    # select the data
    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100_dat, hobo_dat], step_info, buffer=10*60, max_len=20*60)
    lvdt_start = add_minutes(lvdt_step, t0=step_starttime)
    pt100_start = add_minutes(pt100_step, t0=step_starttime)
    hobo_start = add_minutes(hobo_step, t0=step_starttime)

    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100_dat, hobo_dat], step_info, end_buffer=10*60, max_len=-20*60)
    lvdt_end = add_minutes(lvdt_step, t0=step_starttime)
    pt100_end = add_minutes(pt100_step, t0=step_starttime)
    hobo_end = add_minutes(hobo_step, t0=step_starttime)    
    
    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100_dat, hobo_dat], step_info, buffer=0, end_buffer=0, max_len=0)
    lvdt_step = add_minutes(lvdt_step)
    pt100_step = add_minutes(pt100_step)
    hobo_step = add_minutes(hobo_step)

    if len(lvdt_step)>0:
        end_minutes = lvdt_step['minutes'][-1]
    elif len(pt100_step)>0:
        end_minutes = pt100_step['minutes'][-1]
    else:
        end_minutes = 0
        
    hobo_labels = ['Isolated chamber', 'Frost chamber', 'Computer room'] 
    
    #f, axs = plt.subplots(nrows=4, ncols=2, sharex=False, sharey=False,
    #                      gridspec_kw={'height_ratios': [4, 4, 4, 1],
    #                                   'width_ratios':  [1, 3]},
    #                      figsize=(10, 12))

    f = plt.figure(figsize=(15,9))
    gs = matplotlib.gridspec.GridSpec(12, 4)
    plt.subplot(gs.new_subplotspec((1, 0), colspan=1, rowspan=4))      # axes for plot of beginning of step
    plt.subplot(gs.new_subplotspec((1, 1), colspan=3, rowspan=4))      # axes for plot of overview of strain
    plt.subplot(gs.new_subplotspec((5, 1), colspan=3, rowspan=3))      # axes for plot of overview frost room temperatures
    plt.subplot(gs.new_subplotspec((8, 1), colspan=3, rowspan=3))      # axes for plot of overview lab temperatures
    plt.subplot(gs.new_subplotspec((11, 1), colspan=3, rowspan=1))     # axes for legend
    plt.subplot(gs.new_subplotspec((6, 0), colspan=1, rowspan=4))      # axes for plot of end of step
    
    axs = f.axes
    
    # make sure the top 3 axes share the same x-axis
    axs[1].get_shared_x_axes().join(axs[1],axs[2], axs[3])
    
    handles = []
    
    if len(lvdt_step) > 0:
        # Plot the first lvdt
        l1, = axs[1].plot_date(lvdt_step.index, lvdt_step[lvdt_col], '-k', label='LVDT strain', zorder=10)
        handles.extend([l1])
        
        # plot start of time series
        #l, = axs[0].plot(lvdt_start['minutes'], lvdt_start[lvdt_col], '-', c='0.3', lw=0.5, marker='.', ms=4, mec='k', mfc='k', zorder=10)
        l, = axs[0].plot_date(lvdt_start.index, lvdt_start[lvdt_col], '-', c='0.3', lw=0.5, marker='.', ms=4, mec='k', mfc='k', zorder=10)
        
        # plot end of time series
        l, = axs[5].plot_date(lvdt_end.index, lvdt_end[lvdt_col], '-', c='0.3', lw=0.5, marker='.', ms=4, mec='k', mfc='k', zorder=10)
        
        #axs[0].axvline(x=0, ls='--', color='0.65', zorder=-100)
        axs[0].axvline(x=step_starttime, ls='--', color='0.65', zorder=-100)
        axs[5].axvline(x=step_endtime, ls='--', color='0.65', zorder=-100)
        
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
        axs[5].invert_yaxis()

    if len(pt100_step) > 0:
        l3, = axs[2].plot_date(pt100_step.index, pt100_step['T'], '-b', label='Sample temperature', zorder=10)
        handles.extend([l3])
    
    if len(hobo_step) > 0:
        l5, = axs[2].plot_date(hobo_step.index, hobo_step['T(1)'], ls='-', c='orange', marker=None, label=hobo_labels[0], zorder=2)
        #l6, = ax2.plot_date(hobo_step.index, hobo_step['T(2)'], ls='-', c='orange', marker=None, label=hobo_labels[1], zorder=1)
        l7, = axs[3].plot_date(hobo_step.index, hobo_step['T(3)'], ls='-', c='r', marker=None, label=hobo_labels[2])
        #handles.extend([l5, l6, l7])
        handles.extend([l5, l7])
    else:
        if len(lvdt_step)>0:
            axs[3].plot_date(lvdt_step.index[0], lvdt_step[lvdt_col][0], '-', color='none', zorder=10)
        elif len(pt100_step)>0:
            axs[3].plot_date(pt100_step.index[0], pt100_step['T'][0], '-', color='none', zorder=10)
        else:
            pass
        bbox_props = dict(boxstyle="square,pad=0.3", fc="none", ec="none")
        t = axs[3].text(0.5, 0.5, "No data available", 
                          ha="center", va="center", rotation=0,
                          size=15, color='0.75', bbox=bbox_props,
                          transform=axs[3].transAxes)
        
        
    # Plot figure title
    plt.suptitle('Name: {5}, Step: {0},  Date: {1},  Load: {2} kPa,  Temp: {3} C, Duration: {4:.1f} h'.format(step_info['step'],
                                                                              step_info['date'],
                                                                              step_info['load'],
                                                                              step_info['temp'],
                                                                              step_duration_h,
                                                                              step_info['name']),
                 fontsize=14)
                 
    #axs[0].set_xlabel('Time [min]')
    axs[0].set_ylabel('Strain [{0}]'.format(lvdt_unit))
    axs[1].set_ylabel('Strain [{0}]'.format(lvdt_unit))
    axs[5].set_ylabel('Strain [{0}]'.format(lvdt_unit))
    axs[2].set_ylabel('Temperature [C]')
    axs[3].set_ylabel('Temperature [C]')

    # plot common legend for all subplots
    labels = [h.get_label() for h in handles]

    axs[4].axis('off')
    axs[4].legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=4, fontsize=10)
    
    axes = [axs[1], axs[2]]
    if len(hobo_step) > 0:
        axes.append(axs[3])
        
    for ax in axes:
        ax.fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
        ax.grid(True)

    axs[0].fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
    plt.setp(axs[0].xaxis.get_majorticklabels(), rotation=45)   
    
    axs[5].fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
    plt.setp(axs[5].xaxis.get_majorticklabels(), rotation=45)        

    bbox_props = dict(boxstyle="square,pad=0", fc="none", ec="none")
    t = axs[0].text(0.03, 0.97, 'Beginning of step'.format(step_duration_h), 
                          ha="left", va="top", rotation=0,
                          size=9, color='k', bbox=bbox_props,
                          transform=axs[0].transAxes)
    t = axs[5].text(0.03, 0.97, 'End of step', 
                          ha="left", va="top", rotation=0,
                          size=9, color='k', bbox=bbox_props,
                          transform=axs[5].transAxes)
            
    min_lims = [0.1, 1, 0.4, 2, 0.1]
    
    min_ylim(axs[0], min_lims[0])
    
    #if step_info['step'] == 27:
    #    pdb.set_trace()
    
    if len(lvdt_step)>0:
        min_lim = np.sort([lvdt_step[lvdt_col][0], lvdt_step[lvdt_col][-1]])
    elif len(pt100_step)>0:
        min_lim = np.sort([pt100_step['T'][0], pt100_step['T'][-1]])
    min_lim[0] -= 0.2
    min_lim[1] += 0.2
    
    ignore_spikes(axs[1], min_lims[1]/2., min_lim=min_lim)
    min_ylim(axs[1], min_lims[1])
    
    ignore_spikes(axs[2], min_lims[2]/2.)
    min_ylim(axs[2], min_lims[2])            
    
    min_ylim(axs[3], min_lims[3])            
    
    min_ylim(axs[5], min_lims[4])  
    
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
    
    
    
def plot_time_curve(strain, step_info, temp=None, hobo=None, intersect=1, markers=False):
    step_starttime = dt.datetime.combine(step_info['date'], step_info['time1'])
    step_endtime = dt.datetime.combine(step_info['date2'], step_info['time2'])
    step_duration_h = (step_endtime-step_starttime).days*24. + (step_endtime-step_starttime).seconds/3600.
    
    strain = add_minutes(strain, t0=step_starttime)
    
    if hobo is not None:
        if len(hobo) == 0:
            hobo = None
    
    if temp is not None:
        if hobo is not None:
            f, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=False, sharey=False,
                                              gridspec_kw={'height_ratios': [8, 3, 3, 1]},
                                              figsize=(10, 12))
            hobo = add_minutes(hobo, t0=step_starttime)
        else:
            f, (ax1, ax2, ax4) = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False,
                                         gridspec_kw={'height_ratios': [8, 3, 1]},
                                         figsize=(10, 12./14.*11))
            
        temp = add_minutes(temp, t0=step_starttime)
    else:
        f, (ax1, ax4) = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False,
                                         gridspec_kw={'height_ratios': [8, 1]},
                                         figsize=(10, 12./14.*8))
        
    handles = []

    if len(strain)>0:
        l, = ax1.plot(strain['minutes'][1:], strain['eps'][1:], ls='-', color='k', label='LVDT strain')
        handles.append(l)   

        if markers is not None:
            if isinstance(markers, bool) & markers:
                ax1.plot(strain['minutes'][1:], strain['eps'][1:], ls='None', color='k', marker='.', ms=5)
            else:
                ax1.plot(strain['minutes'][1:markers], strain['eps'][1:markers], ls='None', color='k', marker='.', ms=5)
       
    if intersect is None:
        ax1.set_xscale('log')
    else:
        ax1.set_xscale('sqrtlog', intersectx=intersect)

    if len(strain)>0:
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
        
        min_lim = np.sort([strain['eps'][0], strain['eps'][-1]])
        min_lim[0] -= 0.1
        min_lim[1] += 0.1
        ignore_spikes(ax1, 0.2, min_lim=min_lim)
        min_ylim(ax1, 0.1)
        
    if temp is not None:       
        if intersect is None:
            ax2.set_xscale('log')
            if temp.loc[temp.index[0],'minutes'] <= 0:
                temp.loc[temp.index[0],'minutes'] = ax1.get_xlim()[0]
        else:
            ax2.set_xlim(ax1.get_xlim())
            ax2.set_xscale('sqrtlog', intersectx=intersect)

        l, = ax2.plot(temp['minutes'], temp['T'], '-b', label='Sample temperature')
        handles.append(l)
        
        ax2.set_ylabel('Temperature [C]')
        ax2.grid(True)
        ax2.grid(b=True, which='minor', ls='--')
        
        min_ylim(ax2, 0.4)
        
        sqrtlogscale.annotate_xaxis(ax=ax2, intersect=intersect, arrows=None)
        
        ax2.set_xlim(ax1.get_xlim())
        
        # make sure the axes share the same x-axis
        ax1.get_shared_x_axes().join(ax1, ax2)
    
    if hobo is not None:
        if intersect is None:
            ax3.set_xscale('log')
            if hobo['minutes'].iloc[0] <= 0:
                hobo['minutes'].iloc[0] = ax1.get_xlim()[0]
        else:
            ax3.set_xlim(ax1.get_xlim())
            ax3.set_xscale('sqrtlog', intersectx=intersect)

        l, = ax3.plot(hobo['minutes'], hobo['T(3)'], '-g', label='Computer room')
        handles.append(l)
        
        ax3.set_ylabel('Temperature [C]')
        ax3.grid(True)
        ax3.grid(b=True, which='minor', ls='--')
        
        min_ylim(ax3, 2)
        
        sqrtlogscale.annotate_xaxis(ax=ax3, intersect=intersect, arrows=None)
        
        ax3.set_xlim(ax1.get_xlim())
        
        # make sure the axes share the same x-axis
        ax1.get_shared_x_axes().join(ax1, ax2, ax3)
    
    if hobo is not None:
        ax3.set_xlabel('Time (min)')
    elif temp is not None:
        ax2.set_xlabel('Time (min)')
    else:
        ax1.set_xlabel('Time (min)')
    
    # Plot figure title
    plt.suptitle('Name: {0}, Step: {1},  Date: {2},  Load: {3} kPa,  Temp: {4} C, Duration: {5:.1f} h'.format(
                                                                              step_info['name'],
                                                                              step_info['step'],
                                                                              step_info['date'],
                                                                              step_info['load'],
                                                                              step_info['temp'],
                                                                              step_duration_h),
                 fontsize=12)
    
    # plot common legend for all subplots
    labels = [h.get_label() for h in handles]
    ax4.axis('off')
    ax4.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1), ncol=3, fontsize=10)
    
    f.tight_layout(rect=[0,0,0.95,0.9])
    
    sqrtlogscale.annotate_xaxis(ax=ax1, intersect=intersect, arrow_margin=1.05)
    
    return f, ax1
    
def add_minutes(ts, t0=None):
    """Adds timestamp in minutes with the first sample as zero minutes"""
#    if not 'minutes' in ts.columns:
    if len(ts) > 0:
        if t0 is None:
            new_col = (ts.index-ts.index[0]).days*24*60 + (ts.index-ts.index[0]).seconds/60.
        else:
            new_col = (ts.index-t0).days*24*60 + (ts.index-t0).seconds/60.
    else:
        new_col = None
    
    return ts.assign(minutes=new_col)
    
    

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
    
def basic_interpretation(strain, interpret, history, params=None):
    if params is None:
        params = {}
    
    step = interpret['step']
    
    step_id = get_key_index(history, 'step', step)
    
    # Get load of previous step
    if step_id == 0:
        params['sig_n-1'] = 0
    else:
        params['sig_n-1'] = history[step_id-1]['load']
        
    # Get load of current step
    params['sig_n'] = history[step_id]['load']
    
    if len(strain)>0:
        # strain in mm and % at beginning of step
        params['delta_n-1'] = strain['strain'][0]
        params['eps_n-1'] = strain['strain'][0]
        
        # strain in mm and % at end of step
        if 'epsf' in interpret:
            params['delta_n'] = get_strainf(strain, interpret['epsf'])
            params['epsf'] = get_epsf(strain, interpret['epsf'])
        else:
            params['delta_n'] = get_strainf(strain)
            params['epsf'] = get_epsf(strain)

    return params
    
def interpret_k0(params):
    gamma_w = 10 #  kN/m3
    if 'Cv' in params:
        params['K'] = (params['sig_n']-params['sig_n-1'])/((params['eps100']-params['eps0'])/100)
        params['k0'] = params['Cv']*gamma_w/params['K']
        #h_n_1 = (params['h0']-params['delta_n-1'])/1000
        #params['k0'] = (h_n_1*params['delta_n-1']/1000 - (params['delta_n-1']/1000)**2) * gamma_w / (2*params['t50']*60*(params['sig_n']-params['sig_n-1']))
    return params
    
def interpret_iso17892_5(times, values, interpret, history, params=None):
    
    if params is None:
        params = {}
    
    if 'timec' in interpret:
        timec_info = interpret['timec']
    else:
        timec_info = {}
    
    t1 = timec_info.get('t1', None)
    t2 = timec_info.get('t2', None)
    t3 = timec_info.get('t3', None)
    t4 = timec_info.get('t4', None)
    
    if (t1 is not None) & (t2 is not None):
        sqrt_a, sqrt_b = fit_line(times, values, t1, t2, use_points='endpoints', type='sqrt')
        sqrt115_vals = sqrt_a * np.sqrt(times)/1.15 + sqrt_b
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
        params['delta100'] = params['eps100']*params['h0']/100
        params['eps50'] =  (params['eps100']+params['eps0'])/2.
        params['delta50'] = params['eps50']*params['h0']/100
    
#    x = (y-y1)*(x2-x1)/(y2-y1)+x1    
        idx = sign_changes(values-params['eps50']).nonzero()[0]
        params['t50'] = ((params['eps50']-values[idx[0]-1])*(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1]))/(values[idx[0]]-values[idx[0]-1]) + np.sqrt(times[idx[0]-1]))**2.
        
        idx = sign_changes(values-params['eps100']).nonzero()[0]
        params['t100'] = ((params['eps100']-values[idx[0]-1])*(np.sqrt(times[idx[0]])-np.sqrt(times[idx[0]-1]))/(values[idx[0]]-values[idx[0]-1]) + np.sqrt(times[idx[0]-1]))**2.
        
        params['L'] = params['h0']/2*(1-params['eps50']/100)    # Sample height at 50% consolidtation
        Cv = 0.848 * (params['L']/1000)**2 / (params['t90']*60)
        
        params['Cv'] = Cv

    if (t3 is not None) & (t4 is not None):        
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


def plot_legend(ax, info, step, loc='upper right'):
    
    txt = 'Sample: {0}\n'.format(info['sample']['name'])
    txt += 'Depth: {0}\n'.format(info['sample']['depth'])
    txt += 'Load: {0} kPa\n'.format(info['history'][step]['load'])
    txt += 'Temperature: {0} C'.format(info['history'][step]['temp'])
    
    at = AnchoredText(txt,
                  prop=dict(size=10), frameon=True,
                  loc=loc)
    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at)

    
def get_epsf(strain, conf=None, deltat=10):
    if conf is not None:
        deltat = conf['dt']
            
    deltat = dt.timedelta(seconds=deltat*60)
    start_average = strain.index[-1]-deltat
    
    return strain[strain.index >= start_average]['eps'].mean()

    
def get_strainf(strain, conf=None, deltat=10):
    if conf is not None:
        deltat = conf['dt']
            
    deltat = dt.timedelta(seconds=deltat*60)
    start_average = strain.index[-1]-deltat
    
    return strain[strain.index >= start_average]['strain'].mean()
    
    
def get_step_temp(temp, conf):
    if conf['t1'] == -1:
        t1 = temp['minutes'].max()
    else:
        t1 = conf['t1']
        
    return temp[(temp['minutes'] >= conf['t0']) & (temp['minutes'] <= t1)]['T'].mean()

    
def get_raw_strain(strain, t0=0, dt=0.25):
    if t0 == -1:
        t0 = strain['minutes'].max()-dt
    
    avg = strain[(strain['minutes'] >= t0) & (strain['minutes'] <= t0+dt)]['raw_strain'].mean()
    std = strain[(strain['minutes'] >= t0) & (strain['minutes'] <= t0+dt)]['raw_strain'].std()    
        
    return (avg, std)
    
    
def get_key_index(lst, key, val):
    key_list = [d[key] for d in lst]
    try:
        return key_list.index(val)
    except: 
        return None

def plot_consolidation(df):
    f = plt.figure()
    
    # find valid rows of the interpretation file:
    thisdf = df.dropna(subset=['load'], how='any')                # remove any rows with no load specified
    thisdf = thisdf.dropna(subset=['eps100','epsf'], how='all')   # remove rows where neither eps100 nor epsf is specified.

    if 'linestyle' in df.columns:
        # if linestyles are specified
        # plot each segment with the specified line style
        for id, row in thisdf.iterrows():
            if id < len(df)-1:
                l1, = plt.semilogx(thisdf.iloc[id:id+2]['load'], thisdf.iloc[id:id+2]['epsf'], color='k', ls=row['linestyle'], lw=1)
                l2, = plt.semilogx(thisdf.iloc[id:id+2]['load'], thisdf.iloc[id:id+2]['eps100'], color='b', ls=row['linestyle'], lw=1)
            
            if ('marker' in row):
                if (isinstance(row['marker'], str)):
                    marker = row['marker']
                else:
                    marker = None
            else:    
                marker = '.'

            try:
                if np.isnan(marker):
                    pdb.set_trace()
            except:
                pass
                
            plt.semilogx(thisdf.iloc[id]['load'], thisdf.iloc[id]['epsf'], color='k', ls='none', lw=1, marker=marker)
            plt.semilogx(thisdf.iloc[id]['load'], thisdf.iloc[id]['eps100'], color='b', ls='none', lw=1, marker=marker)
        
    else:
        # if linestyles are not specified, just plot default plot.
        l1, = plt.semilogx(thisdf['load'], thisdf['epsf'], '.-k', label=r'$\varepsilon_f$', lw=1)
        l2, = plt.semilogx(thisdf['load'], thisdf['eps100'], '.-b', label=r'$\varepsilon_{100}$', lw=1)
    
    ax = plt.gca()
    ax.invert_yaxis()
    ax.grid(True)
    ax.set_xlabel(r'Stress, $\sigma$ [kPa]')
    ax.set_ylabel(r'Strain, $\varepsilon$ [%]')
    
    for id, row in thisdf.iterrows():
        if row['annotate']:
            ax.annotate(row['txt'], xy=(row['load'], row['epsf']),
                        xytext = (row['load']+row['offset_x'], row['epsf']+row['offset_y']),
                        ha='left', va='center', size=7,
                        #bbox=dict(boxstyle='square', fc='w', ec='w'),
                        arrowprops={'arrowstyle': '-', 'color':'k', 'linewidth':0.5, 'shrinkA':0.05, 'shrinkB':0.05})
    
    ax.legend(handles=[l1,l2], labels=[r'$\varepsilon_f$', r'$\varepsilon_{100}$'], loc='lower left')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    return f
    
    
    
def plot_full_overview(lvdt_dat, pt100_dat, hobo_dat, 
                       history, sample_info, 
                       plot_mV=False, plot_temp=False):

    t0_list = [dt.datetime.combine(h['date'], h['time1']) for h in history]
    tend_list = [dt.datetime.combine(h['date2'], h['time2']) for h in history]
    
    total_duration_h = (tend_list[-1]-t0_list[0]).days + (tend_list[-1]-t0_list[0]).seconds/3600/24.    
    
    steps = [d['step'] for d in history]
    
    # select the data
    #[lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100_dat, hobo_dat], step_info, buffer=60, max_len=120)
    #lvdt_start = add_minutes(lvdt_step, t0=step_starttime)
    #pt100_start = add_minutes(pt100_step, t0=step_starttime)
    #hobo_start = add_minutes(hobo_step, t0=step_starttime)
    
    hobo_labels = ['Isolated chamber', 'Frost chamber', 'Computer room'] 
    
    if plot_mV:
        lvdt_col = 'Volt'
        lvdt_unit = 'mV'
    else:
        lvdt_col = 'eps'
        lvdt_unit = '%'
    
    if plot_temp:
        f = plt.figure(figsize=(12.8,7))
        gs = matplotlib.gridspec.GridSpec(3, 1)
        plt.subplot(gs.new_subplotspec((0, 0), colspan=1, rowspan=2))
        plt.subplot(gs.new_subplotspec((2, 0), colspan=1, rowspan=1))
        axs = f.axes    
        # make sure the top 3 axes share the same x-axis
        axs[0].get_shared_x_axes().join(axs[0],axs[1])
    else:
        f = plt.figure(figsize=(6.4*2,4.8))
        gs = matplotlib.gridspec.GridSpec(1, 1)
        plt.subplot(gs.new_subplotspec((0, 0), colspan=1, rowspan=2))
        axs = f.axes
    
    
    # Plot the first lvdt
    l1, = axs[0].plot_date(lvdt_dat.index, lvdt_dat[lvdt_col], '-k', label='LVDT strain', zorder=10)
    
    if plot_temp:
        lt1, = axs[1].plot_date(pt100_dat.index, pt100_dat['T'], '-k', label='PT100 temp', zorder=10)
        lt1a, = axs[1].plot_date(hobo_dat.index, hobo_dat['T(1)'], ls='-', c='orange', marker=None, label=hobo_labels[0], zorder=2)
    
    #ax = plt.gca()
    
    for id, t in enumerate(t0_list):
        l2 = axs[0].axvline(x=t, ls='--', color='0.65', zorder=+10, label='Step designation')
        if plot_temp:
            lt2 = axs[1].axvline(x=t, ls='--', color='0.65', zorder=+10, label='Step designation')
        
        axs[0].annotate('{0:.0f}'.format(steps[id]), xy=((tend_list[id]-t)/2.+t,1.03),
                    xytext = ((tend_list[id]-t)/2.+t,1.03),
                    xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
                    arrowprops={'arrowstyle': '-', 'color':'none', 'linewidth': 1}, zorder=-100, 
                    ha='center', va='center', size=9
                    )
    #for t in tend_list:
    #    l3 = ax.axvline(x=t, ls='--', color='0.65', zorder=-100, label='End of step')
    axs[0].axvline(x=tend_list[-1], ls='--', color='0.65', zorder=-100)
    if plot_temp:
        axs[1].axvline(x=tend_list[-1], ls='--', color='0.65', zorder=-100)
    
    # invert the direction of both y-axes
    axs[0].invert_yaxis()

    plt.suptitle('Name: {0},  Duration: {1:.1f} days'.format(sample_info['name'],
                                                          total_duration_h),
                 fontsize=14)

    axs[0].set_ylabel('Strain [{0}]'.format(lvdt_unit))
    
    axs[0].set_xlim([t0_list[0] - dt.timedelta(seconds=24*3600),
                 tend_list[-1] + dt.timedelta(seconds=24*3600)])
        
    axs[0].fmt_xdata = matplotlib.dates.DateFormatter('%Y-%m-%d %H:%M:%S')
    axs[0].grid(True)

    min_lim = np.sort([lvdt_dat.loc[lvdt_dat.index>=t0_list[0]].iloc[0][lvdt_col], 
                       lvdt_dat.loc[lvdt_dat.index<=tend_list[-1]].iloc[-1][lvdt_col]])
    min_lim[0] -= 1
    min_lim[1] += 1
    
    axs[0].set_ylim([min_lim[1],min_lim[0]])

    if plot_temp:
        axs[1].set_ylabel('Temperature [C]')
        axs[1].grid(True)
        
    
    #ignore_spikes(ax, 1., min_lim=min_lim)
#    min_ylim(axs[1], min_lims[1])

    # plot common legend for all subplots
    handles = [l1,l2]
    labels = [h.get_label() for h in handles]

    l = axs[0].legend(handles, labels, loc='upper right', fontsize=10)
    l.set_zorder(20)
    #f.tight_layout()
    
    return f    
