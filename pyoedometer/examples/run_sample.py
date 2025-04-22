from importlib import reload
import os
import matplotlib
matplotlib.use('Qt5Agg') # Make sure that we are using QT5
import matplotlib.pyplot as plt

import pandas as pd

import pyoedometer
import pyoedometer.pyoedio as pyoedio
reload(pyoedometer)
reload(pyoedio)

from pyoedometer import *



# Read information from files
info = read_config('./config_sample.yml')

config = info['config']
sample_info = info['sample']
history = info['history']
interpret = info['interpret']
latex_info = info['latex']

lvdt_dat = read_lvdt_data(info['lvdt'], info['sample'])
pt100T = read_pt100_data(info['pt100'])
hoboT = read_hobo_data(info['hobo'])

steps_to_interpret = []
steps_to_plot =  [] 
steps_to_overview = []
steps_to_save = []

# Specify which steps to process or plot
#steps_to_interpret = []
#steps_to_interpret = [13,14,15,16] 
steps_to_interpret = [d['step'] for d in interpret]

#steps_to_plot = []
#steps_to_plot =  [14,14.1,15] 
#steps_to_plot = [13,14,15,16]
steps_to_plot = [d['step'] for d in history]

#steps_to_overview = [] # range(26,37,1) # []
steps_to_overview = [15,16,17,18,19,20,21,22,23]
#steps_to_overview = [d['step'] for d in history]    

#steps_to_save = []
#steps_to_save = [d['step'] for d in history]

plot_full_timeseries = True
close_figs_after_save = False


# create all output dirs specified in the config section of the configuration file.
pyoedio.create_dirs([v for k,v in config.items() if 'path' in k])


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
params_list = []

for step_id, hist in enumerate(history):
    
    hist['name'] = info['sample']['name']
    
    if (hist['step'] in steps_to_overview or 
        hist['step'] in steps_to_plot or 
        hist['step'] in steps_to_interpret or
        hist['step'] in steps_to_save):
        # select the data
        [lvdt_step, pt100_step, hobo_step] = select_data([lvdt_dat, pt100T, hoboT], hist)
        lvdt_step = add_minutes(lvdt_step)
        pt100_step = add_minutes(pt100_step)
        hobo_step = add_minutes(hobo_step)
    else:
        continue
    
    print('Processing {0:.0f} kPa, {1:.0f} C'.format(hist['load'], hist['temp']))
    
    params = {'step': hist['step'], 'load': hist['load'], 'temp': hist['temp'],
              'h0': sample_info['h0']}
    
    if hist['step'] in steps_to_save:
        pyoedio.save_step(config, sample_info, hist, lvdt_step, pt100_step, hobo_step)
    
    if hist['step'] in steps_to_interpret:
        # do interpretation and store interpreted 
        # parameters in dictionary params
        
        h0 = sample_info['h0']   # mm     What thickness should be used here?
    
        int_id = get_key_index(interpret, 'step', hist['step'])
        if int_id is not None:
            
            params.update(basic_interpretation(lvdt_step, interpret[int_id], history, params))
            
            if 'temp' in interpret[int_id]:
                params['temp'] = get_step_temp(pt100_step, interpret[int_id]['temp'])
    
            if 'timec' in interpret[int_id]: 
                timec_info = interpret[int_id]['timec']
            
                if timec_info['type'] == 'iso_sqrt':
                    params = interpret_iso17892_5(lvdt_step['minutes'].values, lvdt_step['eps'].values,
                                                    interpret[int_id], history, params)
                else:
                    raise ValueError('Unknown interpretation method!')
                
                params = interpret_k0(params)

    if hist['step'] in steps_to_overview:
        f = plot_step_overview_hobo2(lvdt_dat, pt100T, hoboT, hist)
        ovname = '{0:02.0f}_{1}_raw_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                info['sample']['name'].replace(' ','-'),
                                                                hist['load'], hist['temp'])
        f.savefig(os.path.join(config['savefigpath'], ovname), dpi=200)
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
        f.savefig(os.path.join(config['savefigpath'], tcname), dpi=200)
        if close_figs_after_save:
            plt.close(f)
    
    params_list.append(params)


# ---------------------------------------------------------------
# SAVE time curve interpretations to excel file
# ---------------------------------------------------------------    
if len(steps_to_interpret) > 0:
    params = pd.DataFrame.from_dict(params_list)
    params = params.set_index('step')

    columns = params.columns
    order = ['load', 'temp', 'eps0', 'eps50', 'eps90', 'eps100', 'epsf', 'epss', 
             'Cv', 'k0', 't50', 't90', 't100', 'sqrt_a', 'sqrt_b', 'log_a', 'log_b', 
             't_start', 't_end']

    order.extend([col for col in columns if col not in order])
    columns = [col for col in order if col in columns]

    writer = pd.ExcelWriter(os.path.join(config['datapath'],'interpretation.xlsx'))
    params.to_excel(writer, columns=columns, sheet_name='results')

    txt = []
    for id, row in params.iterrows():
        txt.append('T={0:.1f}C'.format(row['temp']))
    
    params['annotate'] = 1        
    params['txt'] = txt
    params['offset_x'] = 0
    params['offset_y'] = 0
    
    columns = ['load', 'temp', 'eps100', 'epsf', 'epss', 'annotate', 'txt', 'offset_x', 'offset_y']
    params.to_excel(writer, columns=columns, sheet_name='plotting')
    writer.save()

# ---------------------------------------------------------------
# Plot consolidation curve
# ---------------------------------------------------------------    
if plot_full_timeseries:    
    f = pyoedometer.plot_full_overview(lvdt_dat, pt100T, hoboT, history, sample_info)    
    f.savefig(os.path.join(config['savefigpath'], 'lvdt_strain_full.png'), dpi=200,
              bbox_inches='tight')
    if close_figs_after_save:
        plt.close(f)    
    
# ---------------------------------------------------------------
# Plot consolidation curve
# ---------------------------------------------------------------    

if os.path.exists(latex_info['interpretation_file']):
    params = pd.read_excel(latex_info['interpretation_file'], sheetname='plotting')
    
    f = plot_consolidation(params)

    f.savefig(os.path.join(config['savefigpath'], 'consolidation_curve.png'), 
              dpi=200, bbox_inches = 'tight')
    if close_figs_after_save:
        plt.close(f)

# ---------------------------------------------------------------
# Produce LaTeX file
# ---------------------------------------------------------------    

pyoedio.produce_latex_file(info)
    
    

    


