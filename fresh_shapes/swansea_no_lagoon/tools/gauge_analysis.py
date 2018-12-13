import matplotlib
matplotlib.use('Agg')
import h5py
import uptide
import datetime
#from matplotlib.pylab import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

det_file = '../outputs/diagnostic_detectors.hdf5'
output_file = '../processed/'
df = h5py.File(det_file, 'r+')

data_constituents = ['Q1', 'O1', 'P1', 'S1', 'K1', '2N2', 'MU2', 'N2', 'NU2', 'M2', 'L2', 'T2', 'S2', 'K2', 'M4', 'Z0'] # the ones present in the data
thetis_constituents = ['Q1', 'O1', 'P1', 'K1', 'N2', 'M2', 'S2', 'K2', 'M4']# the ones you want to analyse for
 
tide = uptide.Tides(thetis_constituents)
tide.set_initial_time(datetime.datetime(2003,5,6,8,0))  # make sure this is the same as in the tidal forcing of your run

tidegauge_file = '../inputs/tide_constituents_hyphen_names.dat'
gauge_year = np.loadtxt(tidegauge_file, skiprows=22, usecols=(5,))
gauge_names = np.loadtxt(tidegauge_file, skiprows=22, usecols=(1,), dtype=str)
gauge_amps = np.loadtxt(tidegauge_file, skiprows=22, usecols=np.arange(6, 36, 2))
gauge_phases = np.loadtxt(tidegauge_file, skiprows=22, usecols=np.arange(7, 36, 2))

t = df['time']

dt=100.
spin = int(8*24*60*60/dt)
spin = 0

print(len(t),spin)

def plot_amplitudes(constituent):
    cno_data = (np.array(data_constituents)==constituent).argmax()
    #print(cno_data)
    cno_thetis = (np.array(thetis_constituents)==constituent).argmax()
    #print(cno_thetis)
    plt.figure(figsize=(3,4))
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    maxamp = 0.0
    rmsesum = []
    guagesum = 0
    for name, data in df.items():
        #print(name)
        name=name.split('_')[0]
        if name not in gauge_names:
            #print('Not plotting detector:', name)
            continue

        ind = (gauge_names==name).argmax()
        if gauge_year[ind] < 1975 :
            continue
        amp, pha = uptide.harmonic_analysis(tide, data[spin:,0], t[spin:])
        plt.plot(gauge_amps[ind, cno_data], amp[cno_thetis], '.')
        plt.annotate(name, (gauge_amps[ind, cno_data], amp[cno_thetis]), size=3)
        maxamp = max(maxamp, gauge_amps[ind, cno_data], amp[cno_thetis])
        error = gauge_amps[ind, cno_data] - amp[cno_thetis]
        #print(len(gauge_amps[ind, cno_data]), type(gauge_amps[ind, cno_data]))
        #print(len(amp[cno_thetis]), type(amp[cno_thetis]))
        l2error = (gauge_amps[ind, cno_data] - amp[cno_thetis])**2
        #df.create_dataset('error1', data = error)
        #df.create_dataset('l2error', data = l2error)
        rmsesum.append(l2error)
        guagesum = guagesum + gauge_amps[ind, cno_data]
    rmse = sqrt(sum(rmsesum)/len(rmsesum))
    guageavg = guagesum/len(rmsesum)
    nrmse = rmse/guageavg
    print('{}-constituent rmse = %f'.format(constituent) % rmse)
    print('{}-constituent nrmse = %f'.format(constituent) % nrmse)
    plt.plot([0, 5], [0, 5], 'k')
    plt.xlabel('Gauge amplitude (m)')
    plt.xlim([0,5])
    plt.ylim([0,5])
    plt.ylabel('Thetis amplitude (m)')
    plt.title('{}-constituent'.format(constituent))
    plt.savefig(output_file+'{}-amplitude.png'.format(constituent), dpi = 300)
    return rmse, nrmse

# rmse = {'M2':'','O1':'','K1':'','M4':'','S2':''}
# nrmse = {'M2':'','O1':'','K1':'','M4':'','S2':''}

# [rmse['M2'], nrmse['M2']]=plot_amplitudes('M2')
# [rmse['O1'], nrmse['O1']]=plot_amplitudes('O1')
# [rmse['K1'], nrmse['K1']]=plot_amplitudes('K1')
# [rmse['M4'], nrmse['M4']]=plot_amplitudes('M4')
# [rmse['S2'], nrmse['S2']]=plot_amplitudes('S2')

# print(rmse)
# print(nrmse)
#hf.close()
#show()

rmse = {'M2':''}
nrmse = {'M2':''}

text_file_rmse = open("rmse.txt", "w")
text_file_nrmse = open("nrmse.txt", "w")


for tc in thetis_constituents:
    #rmse.append(tc,)
    [tc_rmse, tc_nrmse] = plot_amplitudes(tc)
    # print(tc)
    # print(type(tc))
    # print(tc_rmse)
    # print(type(str(tc_rmse)))
    # print(len(str(tc_rmse)))
    rmse[tc] = str(tc_rmse)
    nrmse[tc] = str(tc_nrmse)
    #rmse.update({tc,str(tc_rmse)})
    #nrmse.update({tc,str(tc_nrmse)})
    text_file_rmse.write(tc + ' : ' + str(tc_rmse) + '\n')
    text_file_nrmse.write(tc + ' : ' + str(tc_nrmse) + '\n')

