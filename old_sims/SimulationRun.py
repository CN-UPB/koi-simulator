import os
import sys


RUNS  = range(1,2);

for n in RUNS:
        print(' ######## STARTING SIMULATION  ' + str(n) + ' ########');
        f = open('dynParam.ini', 'w')
        f.write('seed-set = ' + str(n+10) + '\n');
        f.close();
        os.system('mkdir -p results');
        os.system('./abstractLTEChannelModel');
        os.system('mv dynParam.ini results/');
        os.system('rm -rf results' + str(n));
        os.system('mv results results_chan_' + str(n));

        print(' ######## END SIMULATION  ' + str(n) + ' ########');

