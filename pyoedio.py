import os
import posixpath
import pandas as pd
import numpy as np

# use double curly braces {{}} for any 'LaTeX' curly braces
# and single curly braces {} for any python formatting 
swfig = r"""
\begin{{sidewaysfigure}}{0}
	\centering
	\includegraphics[width=\linewidth]{{{1}}}{2}{3}
\end{{sidewaysfigure}}"""

#\caption{{{1}}}
#\label{{fig:{2}}}

normalfig = r"""
\begin{{figure}}{0}
	\centering
	\includegraphics[width=\linewidth]{{{1}}}{2}{3}
\end{{figure}}"""

def _process_file(f, func, **kwargs):
    return func(f, **kwargs)

def with_file(f, func, **kwargs):
    if hasattr(f, 'read'):
        return func(f, **kwargs)
    else:
        mode = kwargs.get('filemode', 'a')
        with open(f, mode) as f:
            return func(f, **kwargs)

def create_dirs(dir_list):
    for d in dir_list:
        if not os.path.exists(d):
            os.makedirs(d)
            
def save_step(config, sample_info, hist, lvdt_step=None, pt100_step=None, hobo_step=None):
    fname = '{0:02.0f}_{1}_{2:g}kPa_{3:g}C.xlsx'.format(hist['step'], 
                                                        sample_info['name'].replace(' ','-'),
											            hist['load'], hist['temp'])
    
    print('Saving: {0}'.format(fname))
    
    writer = pd.ExcelWriter(os.path.join(config['savedatapath'],fname))
    
    row = 0
    for d in [sample_info, hist]:
        df = pd.DataFrame.from_dict(d, orient='index')
        df.to_excel(writer, sheet_name='info', header=False, startrow=row)
        row += len(df)+2
    
    if lvdt_step is not None:
        lvdt_step.to_excel(writer, sheet_name='LVDT_data')
        
    if pt100_step is not None:
        pt100_step.to_excel(writer, sheet_name='sample_temperature')
    
    if hobo_step is not None:
        hobo_step.to_excel(writer, sheet_name='room_temperatures')
    
    writer.save()

    
            
def append_sideways_figure(file, figfilepath, caption=None, label=None, loc='[htp]'):
    """
    Parameters
    ----------
    file: open file object
        text file to append figure to
    """
    txt = swfig
    
    if caption is not None:
        cap = '\n\caption{{{0}}}'.format(caption)
    else:
        cap = ''
    
    if label is not None:
        lab = '\n\label{{fig:{0}}}'.format(label)
    else:
        lab = ''
    
    file.write(txt.format(loc, figfilepath, cap, lab)+'\n\n')


def append_normal_figure(file, figfilepath, caption=None, label=None,loc='[htp]'):
    """
    Parameters
    ----------
    file: open file object
        text file to append figure to
    """
    txt = normalfig
    
    if caption is not None:
        cap = '\n\caption{{{0}}}'.format(caption)
    else:
        cap = ''
    
    if label is not None:
        lab = '\n\label{{fig:{0}}}'.format(label)
    else:
        lab = ''
        
    file.write(txt.format(loc, figfilepath, cap, lab)+'\n\n')


    
def append_table(file, df, na_rep='', centering=False, **kwargs):
    
    file.write('\\begin{table}\n')
    
    if centering:
        file.write('\\centering\n')
    
    txt = df.to_latex(na_rep=na_rep, **kwargs)
    file.write(txt)

    file.write('\\end{table}')

    file.write('\n\n')
    


def append_newpage(file):
    file.write('\\clearpage\n\n')

def append_chapter_title(file, title, label=None):
    txt = '\chapter{{{0}}}'.format(title)
    if label is not None:
        txt += '\label{{{0}}}'.format(label)
    
    file.write(txt+'\n\n')
    
    
def produce_latex_file(info):
    config = info['config']
    sample_info = info['sample']
    history = info['history']
    interpret = info['interpret']
    latex_info = info['latex']

    if not latex_info['produce_latex_file']:
        return
        
    steps = sorted([d['step'] for d in history])
    
    with open(latex_info['filename'], 'w') as f:
        
        append_chapter_title(f, 'Sample "{0}"'.format(sample_info['name']))

        tcname = posixpath.join('.', latex_info['figspath'], 'consolidation_curve.png')
        append_normal_figure(f, tcname)
        
        sifo = [['Name', sample_info['name']],
                ['Depth', sample_info['depth']],
                ['Initial height', '{0} mm'.format(sample_info['h0'])],
                ['Diameter', '{0:.1f} mm'.format(sample_info['diameter'])],
                ['Start date', history[0]['date']],
                ['End date', history[-1]['date2']]
                ]
        
        df = pd.DataFrame(sifo, columns=['Parameter','Value'])
        append_table(f, df, centering=True, index=False)
        
        append_newpage(f)

        if os.path.exists(latex_info['interpretation_file']):
            results = pd.read_excel(latex_info['interpretation_file'], sheetname='results')
            
            cols=['step', 'load', 'temp', 'epsf', 'eps100', 'epss', 'Cv']
            dropcols = [col for col in results.columns if col not in cols]
            results.drop(labels=dropcols, axis=1, inplace=True)                
            
            results.rename(columns = {'step':'Step', 'load':'$\\sigma$ [kPa]', 
                                    'temp': 'T [degC]', 'epsf': '$\\varepsilon_f$ [\%]', 
                                    'eps100':'$\\varepsilon_{100}$ [\%]',
                                    'epss':'$C_{\\alpha}$ [\%/lct]',
                                    'Cv':'$c_v$ [$m^2/s$]',
                                    }, inplace = True)
                                    
            def f3(x):
                if np.isnan(x):
                    return ''
                else:
                    return '{0:.3f}'.format(x)
                
            formatters = {'$\\sigma$ [kPa]': lambda x: '{0:.1f}'.format(x),
                        'T [degC]': lambda x: '{0:.1f}'.format(x),
                        '$\\varepsilon_f$ [\%]': lambda x: f3(x),
                        '$\\varepsilon_{100}$ [\%]': lambda x: f3(x),
                        '$C_{\alpha}$ [\%/lct]': lambda x: f3(x),
                        '$c_v$ [$m^2/s$]': lambda x: f3(x)
                        }
                        
            append_table(f, results, centering=True, 
                        na_rep='', index=False, escape=False, formatters=formatters)
        
        else:
            f.write('Interpretation data is missing...\n')
        
        append_newpage(f)

        fname = posixpath.join('.', latex_info['figspath'], 'lvdt_strain_full.png')
        append_sideways_figure(f, fname)

        append_newpage(f)        
        
        for hist in info['history']:
            ovname = '{0:02.0f}_{1}_raw_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                    info['sample']['name'].replace(' ','-'),
                                                                    hist['load'], hist['temp'])
            ovname = posixpath.join('.', latex_info['figspath'], ovname)
            append_sideways_figure(f, ovname)

            tcname = '{0:02.0f}_{1}_timec_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                      info['sample']['name'].replace(' ','-'),
                                                                      hist['load'], hist['temp'])
            tcname = posixpath.join('.', latex_info['figspath'], tcname)
            append_normal_figure(f, tcname)
            
            append_newpage(f)
    
    

