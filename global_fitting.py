import os, sys
from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.signal import residue
import seaborn as sns
import random as rnd
import csv
from scipy import linalg as la
import matplotlib as mpl
#mpl.use('Agg')
import scipy.io as sio
from scipy.linalg import hankel
from scipy.linalg import pinv
from scipy.linalg import eig
from scipy.optimize import minimize                                            
import scipy.constants as sc
import time
import os
import sys


from scipy.signal import chirp, find_peaks, peak_widths
from scipy.interpolate import interp1d 


def read_txt(t,dt,Nmod,Nmodtext,ff,r,DOS_space,grid,version,wvlen, i):
    global delta_w
    dirName='/Users/erglisa/Desktop/SHAHEEN/WVLEN=' + str(wvlen)+'/MOD_FF=' + \
            str(ff) + '_r=' + str(r) + '_Nmod=' + Nmodtext +'_t=' + str(t-1) + \
            '_grid=' + grid + '_DOS=' + DOS_space

    file_path=dirName+'/MOD_FF='+str(ff)+'_r='+str(r)+'_Nmod=' + \
              Nmodtext +'_t='+str(t-1)+'_'+str(i)+'.txt'

    lama = np.loadtxt(file_path, delimiter=',', unpack=True)[0, :]
    norm_data = np.loadtxt(file_path, delimiter=',', unpack=True)[3, :]
    norm_data=norm_data[np.where(lama>=0.2)]
    lama=lama[np.where(lama>=0.2)]
    f=interp1d(lama, norm_data, kind='cubic')
    # plt.plot(lama, norm_data, '.', color='tab:blue')
    inter_lama = np.linspace(np.min(lama), np.max(lama), len(lama)*100) #interpolate lama
    # plt.plot(inter_lama, f(inter_lama), color='tab:red')
    # plt.show()
    

    norm_data=f(inter_lama)
    lama=inter_lama
    peaks, _ = find_peaks(norm_data, height=np.max(norm_data)) #location of max peak
    results_half = peak_widths(norm_data, peaks, rel_height=0.5) #width of peaks
    print(np.around(results_half[0]*2))

    max_peak = np.max(norm_data)
    index = np.where(norm_data==max_peak)
    peak_1 = index-np.around(results_half[0])
    peak_2 = index+np.around(results_half[0])
    lama_1=lama[int(peak_1)]
    lama_2=lama[int(peak_2)]

    print('lama1 is', lama_1)
    print('lama2 is', lama_2)
    print(results_half)
    # plt.plot(lama,norm_data)
    # plt.plot(lama[peaks], norm_data[peaks], "x")
    # plt.hlines(*results_half[1:], color="blue")
    # plt.hlines(results_half[1],lama[int(results_half[2])], lama[int(results_half[3])], color="red")
    print('lama1_interpolated is ', lama[int(results_half[2])])
    print('lama2_interpolated is ', lama[int(results_half[3])])
    delta_w = lama[int(results_half[3])] - lama[int(results_half[2])]
    print('delta w is', delta_w)

    # plt.show()

    index = np.where(norm_data==max_peak)
    print(index)
    # plt.plot(lama[int(peak_1)],norm_data[int(peak_1)], "x")
    # plt.plot(lama[int(peak_2)],norm_data[int(peak_2)], "x")
    # plt.plot(lama, norm_data)
    # plt.show()

# i=0.27
# read_txt(100001, 5.89664e-12 ,50*50 ,'50x50', i, 0.04, '0.4x1.6_0.4_1.6', '400x400', i, 0.5, 17) 

delta_w_list = []
delta_w_mean_list = []
# i=0.15
files = [i+1 for i in range(90)]
filling_fraction = [0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27, 0.3]
# filling_fraction = [0.03, 0.06, 0.09, 0.12]
print(files)
for i in filling_fraction:
    delta_w_list = []
    for version in range(90):
        read_txt(100001, 5.89664e-12 ,50*50 ,'50x50', i, 0.04, '0.4x1.6_0.4_1.6', '400x400', i, 0.5, version+1) 
        delta_w_list.append(delta_w)

    delta_w_mean=np.mean(delta_w_list)
    delta_w_mean_list.append(delta_w_mean)

print(delta_w_mean_list)
plt.plot(filling_fraction, delta_w_mean_list, '.')
plt.xlabel('filling fraction')
plt.ylabel('delta w')
plt.suptitle('peak width with increasing ff')
plt.show()
print('delta_w_mean is ', delta_w_mean)

