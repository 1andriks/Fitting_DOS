import os, sys
from matplotlib import pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from scipy.signal import residue
import seaborn as sns
import random as rnd
import csv
np.set_printoptions(threshold=sys.maxsize)

resdir = 'res2/'

# import geometryproc
# from geometryproc import Geometry
# Geometry(resdir+'geometry.in', xlim=(0.0,0.3), zlim = (0.0, 0.3)).draw(show=False)


import fdtd_utils
from fdtd_utils import tdcmt
from fdtd_utils import myfun
from fdtd_utils import frac_fit
from fdtd_utils import exp_fit
#spectr=tdcmt(resdir+'src.bin',resdir+'ful.bin',1.17541e-12,50000)

fmax = 5
#spectr=tdcmt('/Users/erglisa/Desktop/source_ff=0.27/modes.bin', '/Users/erglisa/Desktop/cuboid_ff=0.27_17/modes.bin',5.89664e-12,100000, 2500)
spectr=tdcmt('/Users/erglisa/Desktop/source/modes.bin', '/Users/erglisa/Desktop/cuboid/modes.bin',5.89664e-12,200000, 25)
nmodes = spectr.smode.shape[1]
# randomly choose 6 modes to plot
np.random.seed(1)

modes = np.random.choice(np.arange(nmodes),size=1, replace=False)
modes.sort()

#plt.figure(figsize=(5,3), dpi = 200)
plt.xlim(0, fmax)
# plt.ylim(0, 1.8)
plt.xlabel('$1/\lambda$')
#plt.title('modes %s' % str(modes))
FTED = np.array([])


#while error>=tolerance:
# for i in range(nmodes):


i=0
nmodes=25

COEFF= []
coeff_matrix=np.zeros((20,2))

for i in range(nmodes):
    poles =2
    th=0.5
    spectr.fitone(i,th,fmax, N=poles, plot=False)
    # print('nonzero indexes of fted_full are ', np.nonzero(spectr.fted_full))
    limit = np.where(spectr.f<=fmax) 
    spectr.smode=spectr.smode[limit]
    #print('size of smode is', spectr.smode.shape)
    masking = spectr.fted_full !=0
    
    error = np.sum((np.abs(spectr.smode[:,i])**2-np.abs(spectr.fted_full)**2)**2)
    #error = np.sum((np.abs(spectr.smode[masking,i])**2-np.abs(spectr.fted_full[masking])**2)**2)
        #print('error at ' + str(i) + 'th mode is', error)
    # print('spectr.mode is', spectr.smode[masking,i])
    # print('spectr.mode length  is', len(spectr.smode[:,i]))
    # print('spectr.mode length with masking is', len(spectr.smode[masking,i]))
    # print('ftd_full length is', len(spectr.fted_full[masking]))

    while error>30:
        #th=th+0.01
        poles=poles+1
        spectr.fitone(i,th,fmax, N=poles, plot=False)
        spectr.smode=spectr.smode[limit]
        error = np.sum((np.abs(spectr.smode[:,i])**2-np.abs(spectr.fted_full)**2)**2)
            #print('error at ' + str(i) + 'th mode is', error)

    plt.plot(spectr.f[limit],np.abs(spectr.smode[:,i])**2,'b',spectr.w,np.abs(spectr.fted)**2,'r')

    fted_T=spectr.fted_full
    FTED = np.column_stack((FTED,fted_T)) if FTED.size else fted_T
    #print('FTED is', FTED)
    print('len of a.a is', len(spectr.a.a))
    # coeff_matrix.resize((len(spectr.a.a),2), refcheck=False) ##either of these work
    coeff_matrix = np.resize(coeff_matrix, (len(spectr.a.a),2)) #either of these work
    # print(coeff_matrix)

    coeff_matrix[:,0]=spectr.a.a
    coeff_matrix[:,1]=spectr.a.b
    print('len of coeff matrix is', len(coeff_matrix))
    COEFF.append(coeff_matrix)
    print(len(spectr.a.a))
# print('COEFF is ', COEFF)
with open('coefficients.npz', 'wb') as f:
    np.save(f, COEFF)

#     plt.show()

# # # print('fted_T shape is', fted_T.shape)
# # # print('fted shape is', spectr.fted.shape)



# fted_DOS=np.sum(np.abs(FTED)**2, axis=1)
# # print(fted_DOS)
# # print('fted_Dos is', fted_DOS)
# plt.figure(figsize=(18,8))
# plt.plot(spectr.w_full,fted_DOS+0.2,'tab:red', '-')
with open('coeff.txt', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(spectr.a.a, spectr.a.b))
f.close()

# plt.xlim(0,5)
# # plt.ylim(0,2)
# plt.plot(spectr.f, spectr.dos, 'tab:blue')

# plt.show()
# # spectr.fitone(24,0.5,fmax, N=poles)





# fig3, ax3 = plt.subplots(2, 1)
# fig3.suptitle('ff=', fontsize=16, weight='bold')   #fontname="Times New Roman"#, style='italic' )
# ax3[0].plot(spectr.f, spectr.dos, linewidth=2, label='Silicon')
# ax3[0].plot(spectr.f, spectr.pdss, label='Random Structure')
# ax3[0].set(title='Density of states ', xlabel="1/lambda", ylabel="Amplitude")
# ax3[0].legend()
# ax3[1].plot(spectr.f, spectr.dos)
# ax3[1].set(title='Density of states (Random Structure)', xlabel="1/lambda", ylabel="Norm. Amplitude")
# fig3.tight_layout()
# fig3.subplots_adjust(top=0.86)

# # plt.savefig(self.dirName+'/MOD_FF='+str(self.ff)+'_r='+str(self.r) + '_Nmod='+self.Nmodtext+'_t='+str(self.t-1)+'_'+str(i)+'.png', dpi=500)
# plt.show()
# plt.close()    


#print (np.load("FT.npy"))