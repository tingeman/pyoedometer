from importlib import reload
import os
import matplotlib
matplotlib.use('Qt5Agg') # Make sure that we are using QT5
import matplotlib.pyplot as plt

import pandas as pd

try:
    import ipdb as pdb
except ImportError:
    import pdb

import pyoedometer
import pyoedometer.pyoedio as pyoedio
reload(pyoedometer)
reload(pyoedio)

from pyoedometer import *



# Read information from files
info = read_config('sample_data/config_sample.yml')
#config = info['config']          # configuration of the measurement setup
sample_info = info['sample']     # information about the sample 
history = info['history']        # history of load steps and nominal (set point) temperatures
interpret = info['interpret']    # interpretation parameters for each step
latex_info = info['latex']       # information for LaTeX file generation

lvdt_dat = read_lvdt_data(info['lvdt'], info['sample'])    # Read LVDT data
pt100T = read_pt100_data(info['pt100'])                    # Read PT100 temperature data
hoboT = read_hobo_data(info['hobo'])                       

#------------------------------------------------------------
# Set up processing parameters
#------------------------------------------------------------

if True:
    # Get processing parameters from the configuration file
    steps_to_interpret = info['processing']['steps_to_interpret']  
    steps_to_plot =  info['processing']['steps_to_plot_timecurve'] 
    steps_to_overview = info['processing']['steps_to_plot_overview'] 
    steps_to_save = info['processing']['steps_to_save']            

    plot_full_timeseries = info['processing']['plot_full_timeseries']
    plot_full_timeseries_temp = info['processing']['plot_temperature_with_full_timeseries'] 
    close_figs_after_save = info['processing']['close_figs_after_saving']       

else:
    # or specify processing parameters manually
    #
    # These lists can be empty, in which case no steps will be processed.
    # If you want to process all steps, you can use the commented lines below.
    # The steps are specified by their step number, which is defined in the history section of the configuration file.
    #
    # Often you will want to work on single step, and iteratively change the parameters in the 
    # interpretation section, until you finde the best fit for the data.

    steps_to_interpret = []
    steps_to_plot =  [] 
    steps_to_overview = []
    steps_to_save = []

    # Specify which steps to process or plot
    #steps_to_interpret = []
    steps_to_interpret = [22, 25] 
    #steps_to_interpret = [d['step'] for d in interpret]   # interpret all steps defined in `interpret`

    #steps_to_plot = []
    steps_to_plot =  [22, 25]
    #steps_to_plot = [7,9,10,11,12,17]
    #steps_to_plot = [d['step'] for d in history]          # plot all steps defined in `history`

    #steps_to_overview = [] # range(26,37,1) # []
    steps_to_overview = [22, 25]
    #steps_to_overview = [d['step'] for d in history]      # plot overview of all steps defined in `history`

    #steps_to_save = []
    steps_to_save = [22, 25]
    #steps_to_save = [d['step'] for d in history]          # save all steps defined in `history`

    plot_full_timeseries = False         # If True, plot the full time series of the LVDT data
    plot_full_timeseries_temp = False    # If True, add temperature data to the full time series plot
    close_figs_after_save = False        # If True, close figures after saving them to file


def convert_input_to_list_of_steps(steps, step_list):
    """Convert argument to a list of integers."""
    if isinstance(steps, str):
        if steps.lower() == 'none':
            return []
        elif steps.lower() == 'all':
            return [d['step'] for d in step_list]
        else:
            raise ValueError(f"Unknown value for steps: {steps}. "
                                "Use 'none', 'all', or a list of step numbers.")
    elif isinstance(steps, list):
        return steps
    else:
        raise TypeError("Steps must be a string or a list of integers.")


steps_to_interpret = convert_input_to_list_of_steps(steps_to_interpret, interpret)
steps_to_plot =  convert_input_to_list_of_steps(steps_to_plot, history)
steps_to_overview = convert_input_to_list_of_steps(steps_to_overview, history)
steps_to_save = convert_input_to_list_of_steps(steps_to_save, history)

#------------------------------------------------------------

# create all output dirs specified in the config section of the configuration file.
pyoedio.create_dirs([v for k,v in info['paths'].items() if 'path' in k])

# ---------------------------------------------------------------
# DO:
#     Data selection
#     Save step
#     Interpret time curve
#     Plot step overview
#     Plot time curve
# for each specified step.
# ---------------------------------------------------------------

# Iterate over all steps
step_params = {}
step_timeseries = {}


# ---------------------------------------------------------------
# Do basic processing of all time steps
# ---------------------------------------------------------------

for step_id, hist in enumerate(history):

    print('Basic processing {0:.0f} kPa, {1:.0f} C'.format(hist['load'], hist['temp']))

    # Estract data for the current step
    hist['name'] = info['sample']['name']
    
    # Select data for the current step
    [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100T, hoboT], hist)

    # add a new column to each dataframe, containing time in minutes since the start of the load step
    lvdt_step = add_minutes(lvdt_step)
    pt100_step = add_minutes(pt100_step)
    hobo_step = add_minutes(hobo_step)

    params = {'step': hist['step'], 
              'load': hist['load'], 
              'temp': hist['temp'], 
              'h0': sample_info['h0']}

    h0 = sample_info['h0']   # mm     What thickness should be used here?

    int_id = get_key_index(interpret, 'step', hist['step']) # get the index of the 
    if int_id is not None:
        # Extract basic parameters from the strain data
        params.update(basic_interpretation(lvdt_step, interpret[int_id], history, params))

        # Calculate the actual temperature for the current step
        if 'temp' in interpret[int_id]:
            params['temp'] = get_step_temp(pt100_step, interpret[int_id]['temp'])

    # store the parameters for later use and processing
    step_params[hist['step']] = params
    # store the extracted timeseries for later use and processing
    step_timeseries[hist['step']] = {'lvdt': lvdt_step,
                                     'pt100': pt100_step,
                                     'hobo': hobo_step}

# ----------------------------------------------------------------------
# Do detailed processing, plotting, and saving of specified time steps
# ----------------------------------------------------------------------

for step_id, hist in enumerate(history):
    
    # Estract data for the current step
    hist['name'] = info['sample']['name']
    
    # Check if the step is in the lists of steps to process
    # if it is, extract relevant data for that step
    if (hist['step'] in steps_to_overview or 
        hist['step'] in steps_to_plot or 
        hist['step'] in steps_to_interpret or
        hist['step'] in steps_to_save):

        # Obtain the data for the current step
        lvdt_step = step_timeseries[hist['step']]['lvdt']
        pt100_step = step_timeseries[hist['step']]['pt100']
        hobo_step = step_timeseries[hist['step']]['hobo']

        # Obtain the basic parameters for the current step
        params = step_params[hist['step']]
    else:
        continue
    
    print('Detailed processing {0:.0f} kPa, {1:.0f} C'.format(hist['load'], hist['temp']))
   
    if hist['step'] in steps_to_save:
        pyoedio.save_step(info['paths'], sample_info, hist, lvdt_step, pt100_step, hobo_step)
    
    if hist['step'] in steps_to_interpret:
        # do interpretation and store interpreted 
        # parameters in dictionary params
        
        h0 = sample_info['h0']   # mm     What thickness should be used here?
    
        int_id = get_key_index(interpret, 'step', hist['step']) # get the index of the current step in the interpretation dictionary
        if int_id is not None:
            if 'timec' in interpret[int_id]: 

                # Get the intial strain for the current step
                # If this is the first step of the experiment, eps0 should be 0.0
                # If this is not the first step, use epsf from the previous step
                if hist['step'] == history[0]['step']:
                    params['eps0'] = 0
                else:
                    prev_step = history[step_id - 1]['step']  # Get the previous step number
                    if prev_step in step_params:
                        # Use epsf from the previous step as eps0 for the current step
                        params['eps0'] = step_params[prev_step]['epsf']
                    else:
                        raise ValueError(f"Previous step {prev_step} not found in step_params. Could not determine eps0.")
                

                # ---------------------------------------------------------------
                # Do the actual interpretation of the time curve 
                # ---------------------------------------------------------------
                
                # Get the interpretation parameters for the current step
                timec_info = interpret[int_id]['timec']
            
                if timec_info['type'] == 'iso_sqrt':
                    # Interpret the time curve using ISO 17892-5 method
                    params = interpretation_iso17892_5(lvdt_step['minutes'].values, 
                                                       lvdt_step['eps'].values,
                                                       interpret[int_id], 
                                                       history, 
                                                       params)

                    # Adjust the minutes in the LVDT step data to account for small (sub-second)
                    # initial time offset, as estimated by the sqrt(t) interpretation.
                    lvdt_step['minutes'] = lvdt_step['minutes'] - params.get('dt_eps0', 0)
                    pt100_step['minutes'] = pt100_step['minutes'] - params.get('dt_eps0', 0)
                    hobo_step['minutes'] = hobo_step['minutes'] - params.get('dt_eps0', 0)
                    
                elif timec_info['type'] == 'bh':
                    # Interpret the time curve using Brinch Hansen method
                    t1 = interpret[int_id]['timec'].get('t1')
                    t2 = interpret[int_id]['timec'].get('t2')
                    t3 = interpret[int_id]['timec'].get('t3')
                    t4 = interpret[int_id]['timec'].get('t4')
                    h0 = params.get('h0', sample_info['h0'])

                    params = interpretation_BrinchHansen(lvdt_step['minutes'].values, 
                                                         lvdt_step['eps'].values, 
                                                         h0,t1, t2, t3, t4)

                else:
                    raise ValueError('Unknown interpretation method!')
                
                params = calculate_k0(params)


    if hist['step'] in steps_to_overview:
        f = plot_step_overview_hobo2(lvdt_dat, pt100T, hoboT, hist)
        ovname = '{0:02.0f}_{1}_raw_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                info['sample']['name'].replace(' ','-'),
                                                                hist['load'], hist['temp'])
        f.savefig(os.path.join(info['paths']['savefigpath'], ovname), dpi=200)
        if close_figs_after_save:
            plt.close(f)

            
    if hist['step'] in steps_to_plot:
        intersect = None
        if interpret is not None:
            int_id = get_key_index(interpret, 'step', hist['step'])
            if (int_id is not None):
                if ('timec' in interpret[int_id]):
                    intersect = interpret[int_id]['timec']['intersect']
            
        if 't100' in params :
            f, ax = plot_time_curve(lvdt_step, hist, temp=pt100_step, intersect=params['t100'], markers=30)
        else:
            f, ax = plot_time_curve(lvdt_step, hist, temp=pt100_step, hobo=hobo_step, intersect=intersect, markers=30)
        
        if ('sqrt_a' in params) & ('sqrt_b' in params):
            plot_fit_line(params['sqrt_a'], params['sqrt_b'], ax, type='sqrt', num=100, intersect=params['t100'])
            plot_fit_line(params['sqrt_a']/1.15, params['sqrt_b'], ax, type='sqrt', num=100, ls='--r', intersect=params['t100'])
            annotate_yaxis(ax, params)
            
            if ('log_a' in params) & ('log_b' in params):
                plot_fit_line(params['log_a'], params['log_b'], ax, type='log10', num=100, ls='--r', intersect=params['t100'])
        elif ('log_a' in params) & ('log_b' in params):
            plot_fit_line(params['log_a'], params['log_b'], ax, type='log10', num=100, ls='--r', intersect=intersect)
            
        tcname = '{0:02.0f}_{1}_timec_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                  info['sample']['name'].replace(' ','-'),
                                                                  hist['load'], hist['temp'])
        f.savefig(os.path.join(info['paths']['savefigpath'], tcname), dpi=200)
        if close_figs_after_save:
            plt.close(f)
    
    step_params[hist['step']].update(params)  # Update the step parameters with the interpreted values


# ---------------------------------------------------------------
# SAVE time curve interpretations to excel file
# ---------------------------------------------------------------    
if len(steps_to_interpret) > 0:
    #params = pd.DataFrame.from_dict(params_list)
    #params = params.set_index('step')

    params = pd.DataFrame.from_dict(step_params, orient='index')
    params = params.reset_index(drop=False)

    columns = params.columns
    order = ['load', 'temp', 'eps0', 'eps50', 'eps90', 'eps100', 'epsf', 'epss', 
             'Cv', 'k0', 't50', 't90', 't100', 'sqrt_a', 'sqrt_b', 'log_a', 'log_b', 
             't_start', 't_end']

    order.extend([col for col in columns if col not in order])
    columns = [col for col in order if col in columns]

    with pd.ExcelWriter(os.path.join(info['paths']['datapath'],'interpretation.xlsx')) as writer:
        params.to_excel(writer, columns=columns, sheet_name='results')

        txt = []
        for id, row in params.iterrows():
            txt.append('T={0:.1f}C'.format(row['temp']))
        
        params['annotate'] = 1        
        params['txt'] = txt
        params['offset_x'] = 0
        params['offset_y'] = 0
        
        columns = ['load', 'temp', 'eps100', 'epsf', 'epss', 'annotate', 'txt', 'offset_x', 'offset_y']
        #params.to_excel(writer, columns=columns, sheet_name='plotting')

        # only write the columns that actually exist
        existing = [c for c in columns if c in params.columns]
        missing = set(columns) - set(existing)
        if missing:
            # optional: warn the user which columns were skipped
            print(f"Warning: skipping missing plotting columns: {sorted(missing)}")
        params.to_excel(writer, columns=existing, sheet_name='plotting')


# ---------------------------------------------------------------
# Plot consolidation curve
# ---------------------------------------------------------------    
if plot_full_timeseries:    
    f = pyoedometer.plot_full_overview(lvdt_dat, pt100T, hoboT, history, sample_info, plot_temp=plot_full_timeseries_temp)

    if os.path.exists(latex_info['interpretation_file']):
        params = pd.read_excel(latex_info['interpretation_file'], sheet_name='plotting')
        
        if plot_full_timeseries_temp:
            pass
            # Tyl = [np.floor(params['temp'].min()-1), np.ceil(params['temp'].max()+1)]
            # f.axes[1].set_ylim(Tyl)
        
            # Eyl = [np.ceil(params['epsf'].max()+1), np.floor(params['epsf'].min()-1)]
            # f.axes[0].set_ylim(Eyl)
    
    f.savefig(os.path.join(info['paths']['savefigpath'], 'lvdt_strain_full.png'), dpi=200,
              bbox_inches='tight')
    if close_figs_after_save:
        plt.close(f)    
    
# ---------------------------------------------------------------
# Plot consolidation curve
# ---------------------------------------------------------------    

if os.path.exists(info['paths']['interpretation_file']):
    params = pd.read_excel(info['paths']['interpretation_file'], sheet_name='plotting')
    
    f = plot_consolidation(params)
    
    f.savefig(os.path.join(info['paths']['savefigpath'], 'consolidation_curve.png'), 
          dpi=300, bbox_inches = 'tight')
          
    yl = f.axes[0].get_ylim()
    f.axes[0].set_ylim([17.5, yl[1]])
    
    f.savefig(os.path.join(info['paths']['savefigpath'], 'consolidation_curve_yl.png'), 
              dpi=300, bbox_inches = 'tight')
              
    if close_figs_after_save:
        plt.close(f)

# ---------------------------------------------------------------
# Produce LaTeX file
# ---------------------------------------------------------------    

pyoedio.produce_latex_file(info)


# ---------------------------------------------------------------
# Completed processing
# ---------------------------------------------------------------    

print(" ")
print('Processing completed.')
print('Results saved to: {0}'.format(info['paths']['datapath']))
print('Figures saved to: {0}'.format(info['paths']['savefigpath']))
print('LaTeX file saved to: {0}'.format(latex_info['filename']))
print(" ")
print("If figures are not shown, you can run issue the command `plt.show(block=False)` to display them.")