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
                                                        sample_info['name'].replace(' ','_'),
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


    
def append_table(file, df, na_rep='', small=True, sideways=False, centering=False, loc='[htp]', **kwargs):
    
#    if sideways:
#        file.write('\\begin{{sidewaystable}}{0}\n'.format(loc))
#    else:
#        file.write('\\begin{{table}}{0}\n'.format(loc))
    
    if sideways:
        file.write('\\begin{turn}{90}\n')
        file.write('\\begin{minipage}{0.9\\textheight}\n')
        file.write('\\begin{table}[H]\n')
    else:
        file.write('\\begin{{table}}{0}\n'.format(loc))
    
    if centering:
        file.write('\\centering\n')
        
    if small:
        file.write('\\small\n')
    
    txt = df.to_latex(na_rep=na_rep, **kwargs)
    file.write(txt)

    file.write('\\end{table}')

    if sideways:
        file.write('\\end{minipage}\n')
        file.write('\\end{turn}\n')
    
#    if sideways:
#        file.write('\\end{sidewaystable}')
#    else:
#        file.write('\\end{table}')

    file.write('\n\n')

    
def append_section_heading(file, secname, label=None, indent=0):
    if label is None:
        file.write('{0}\\section{{{1}}}\n\n'.format( ' '*indent,secname))
    else:
        file.write('{0}\\section{{{1}}}\label{{{2}}}\n\n'.format( ' '*indent,secname,label))
    

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
        
        append_chapter_title(f, 'Sample "{0}"'.format(sample_info['name']), 
                             label='app:{0}'.format(sample_info['nickname'].lower().replace(' ','_')))

        append_section_heading(f,'Consolidation curve')
        
        tcname = posixpath.join('.', latex_info['figspath'], 'consolidation_curve.png')
        append_normal_figure(f, tcname, loc='[hp!]')
        
        sifo = [['Name', sample_info['name']],
                ['Depth', sample_info['depth']],
                #['Initial height', '{0} mm'.format(sample_info['h0'])],
                #['Diameter', '{0:.1f} mm'.format(sample_info['diameter'])],
                ['Start date', history[0]['date']],
                ['End date', history[-1]['date2']]
                ]
        
        df = pd.DataFrame(sifo, columns=['Parameter','Value'])
        append_table(f, df, centering=True, index=False)
        
        append_newpage(f)
        
        cfname, cfext = os.path.splitext(os.path.basename(latex_info['filename']))
        cfname = cfname + '_classification.tex'
        
        f.write('\\IfFileExists{{./{0}}}{{\n'.format(cfname))
        append_section_heading(f, 'Classification parameters', indent=4)
        f.write('    \\input{{./{0}}}\n'.format(cfname))
        f.write('    \\clearpage\n')
        f.write('}{}\n\n')
        
        append_section_heading(f, 'Overview of load steps and interpreted results')
        
        if os.path.exists(latex_info['interpretation_file']):
            results = pd.read_excel(latex_info['interpretation_file'], sheetname='results')
            
            def fixed(x, sig=3):
                txt = '{{0:.{0:d}f'.format(sig) + '}'
                if np.isnan(x):
                    return ''
                else:
                    return txt.format(x)
                
            def scientific(x, sig=3):
                txt = '\\num{{{0:.'+'{0:d}e'.format(sig) + '}}}'
                if np.isnan(x):
                    return ''
                else:
                    return txt.format(x)
                        
            formatters = {'step':   lambda x: fixed(x, 0),
                          'load':   lambda x: fixed(x, 1),
                          'temp':   lambda x: fixed(x, 1),
                          'eps0':   lambda x: fixed(x, 3),
                          'eps50':  lambda x: fixed(x, 3),
                          'eps90':  lambda x: fixed(x, 3),
                          'eps100': lambda x: fixed(x, 3),
                          'epsf':   lambda x: fixed(x, 3),
                          't50':    lambda x: fixed(x, 3),
                          't90':    lambda x: fixed(x, 3),
                          't100':   lambda x: fixed(x, 3),
                          'epss':   lambda x: fixed(x, 3),
                          'Cv':     lambda x: scientific(x, 3),
                          'K':      lambda x: fixed(x, 0),
                          'k0':     lambda x: scientific(x, 3)
                         }
            

            columns = ['step',
                       'load',
                       'temp', 
                       'eps0',
                       'eps50', 
                       'eps100',
                       'epsf',
                       'epss', 
                       'Cv',   
                       'K',    
                       'k0']
            
            dropcols = [col for col in results.columns if col not in columns]
            df = results.drop(labels=dropcols, axis=1, inplace=False)                
            cols = [col for col in columns if col in df.columns]
            
            # Since we are using siunitx to typeset the columns, 
            # non-numerical cell contents must be braced...
            headers = ['{Step}', 
                       '{$\\sigma$ [\si{kPa}]}', 
                       '{$T$ [\si{\celsius}]}',
                       '{$\\varepsilon_{0}$ [\%]}',
                       '{$\\varepsilon_{50}$ [\%]}',
                       '{$\\varepsilon_{100}$ [\%]}',
                       '{$\\varepsilon_{f}$ [\%]}',
                       '{$C_{\\alpha}$ [\%/lct]}',
                       '{$c_v$ [\si{m^2/s}]}',
                       '{$K$ [\si{kPa}]}',
                       '{$k_0$ [\si{m/s}]}']                       

                       
            heads = [headers[id] for id, col in enumerate(columns) if col in df.columns]
                        
            append_table(f, df, centering=True, sideways=True, loc='[H]', columns=cols, header=heads, 
                        na_rep='', index=False, escape=False, formatters=formatters,
                        column_format='c'*len(heads))

            
            
#            columns = ['step',
#                       'load',
#                       'temp', 
#                       'eps0',
#                       'eps50', 
#                       'eps90',
#                       'eps100',
#                       'epsf',
#                       't50',
#                       't90',  
#                       't100']
#            
#            dropcols = [col for col in results.columns if col not in columns]
#            df = results.drop(labels=dropcols, axis=1, inplace=False)                
#            cols = [col for col in columns if col in df.columns]
#            
#            headers = ['Step', 
#                       '$\\sigma$ [kPa]', 
#                       '$T$ [\si{\celsius}]',
#                       '$\\varepsilon_{0}$ [\%]',
#                       '$\\varepsilon_{50}$ [\%]',
#                       '$\\varepsilon_{90}$ [\%]',
#                       '$\\varepsilon_{100}$ [\%]',
#                       '$\\varepsilon_{f}$ [\%]',
#                       '$t_{50}$ [min]',
#                       '$t_{90}$ [min]',
#                       '$t_{100}$ [min]']
#                       
#            heads = [headers[id] for id, col in enumerate(columns) if col in df.columns]
#                        
#            append_table(f, df, centering=True, sideways=True, loc='[hp!]', columns=cols, header=heads, 
#                        na_rep='', index=False, escape=False, formatters=formatters)
#
#            append_newpage(f)
#                        
#            columns = ['step',
#                       'load',
#                       'temp', 
#                       'epss', 
#                       'Cv',   
#                       'K',    
#                       'k0']
#
#            dropcols = [col for col in results.columns if col not in columns]
#            df = results.drop(labels=dropcols, axis=1, inplace=False)                
#            cols = [col for col in columns if col in df.columns]
#                       
#            headers = ['Step', 
#                       '$\\sigma$ [kPa]', 
#                       '$T$ [\si{\celsius}]',
#                       '$C_{\\alpha}$ [\%/lct]',
#                       '$c_v$ [$m^2/s$]',
#                       '$K$ [$kPa$]',
#                       '$k_0$ [$m/s$]']
#                       
#            heads = [headers[id] for id, col in enumerate(columns) if col in df.columns]
#                        
#            append_table(f, df, centering=True, loc='[hp!]', columns=cols, header=heads, 
#                        na_rep='', index=False, escape=False, formatters=formatters)

                        
        else:
            f.write('Interpretation data is missing...\n')
        
        append_newpage(f)

        fname = posixpath.join('.', latex_info['figspath'], 'lvdt_strain_full.png')
        append_sideways_figure(f, fname)

        append_newpage(f)        
        
        f.write('\\iftimecurves\n')
        f.write('\\fakesection{Time curves}\n\n')
        
        for hist in info['history']:
            ovname = '{0:02.0f}_{1}_raw_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                    info['sample']['name'].replace(' ','_'),
                                                                    hist['load'], hist['temp'])
            ovname = posixpath.join('.', latex_info['figspath'], ovname)
            append_sideways_figure(f, ovname)

            tcname = '{0:02.0f}_{1}_timec_{2:g}kPa_{3:g}C.png'.format(hist['step'], 
                                                                      info['sample']['name'].replace(' ','_'),
                                                                      hist['load'], hist['temp'])
            tcname = posixpath.join('.', latex_info['figspath'], tcname)
            append_normal_figure(f, tcname)
            
            append_newpage(f)
    
        f.write('\\fi\n')

