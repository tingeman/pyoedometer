import os
from subprocess import *

scripts = ['python ./sample_1.py',
           'python ./sample_2.py',
           'python ./sample_3.py',
           'python ./sample_4.py',
           'python ./sample_5.py',
           'python ./sample_6.py']

#scripts = ['python ./sample_5.py',
#           'python ./sample_6.py']


# 1) run all samples for all steps
# 2) introduce step 14.1 for samples 3+4
# 3) run samples 3+4 doing only interpret for all steps including 14.1
# 4) update interpret_sample_x for all samples 1-6
# 5) remove step 14.1 for sample 3+4
# 6) run interpret for all samples, in order to write correct latex tables.



for script in scripts:
    print('Running {0}\n'.format(script))
    #p = Popen(script, shell=True, stdin=PIPE, stdout=PIPE)
    #output = p.communicate()
    #print(output[0])
    os.system(script)