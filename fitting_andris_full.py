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
import pyfftw
import multiprocessing
from scipy import fftpack

from scipy.signal import chirp, find_peaks, peak_widths


def read_txt(t,dt,Nmod,Nmodtext,ff,r,DOS_space,grid,version,wvlen, i):
    global lama_1
    global lama_2
    dirName='/Users/erglisa/Desktop/SHAHEEN/WVLEN=' + str(wvlen)+'/MOD_FF=' + \
    str(ff) + '_r=' + str(r) + '_Nmod=' + Nmodtext +'_t=' + str(t-1) + \
    '_grid=' + grid + '_DOS=' + DOS_space

    file_path=dirName+'/MOD_FF='+str(ff)+'_r='+str(r)+'_Nmod=' + \
    Nmodtext +'_t='+str(t-1)+'_'+str(i)+'.txt'

    lama = np.loadtxt(file_path, delimiter=',', unpack=True)[0, :]
    norm_data = np.loadtxt(file_path, delimiter=',', unpack=True)[3, :]
    norm_data=norm_data[np.where(lama>=1.2)]
    lama=lama[np.where(lama>=1.2)]
    peaks, _ = find_peaks(norm_data, height=np.max(norm_data)) #location of max peak
    results_half = peak_widths(norm_data, peaks, rel_height=0.5) #width of peaks
    print(np.around(results_half[0]*2))

    max_peak = np.max(norm_data)
    index = np.where(norm_data==max_peak)
    peak_1 = index-np.around(results_half[0]*2)
    peak_2 = index+np.around(results_half[0]*2)
    
    lama_1=lama[int(peak_1)]
    lama_2=lama[int(peak_2)]
    print('lama1 is', lama_1)
    print('lama2 is', lama_2)



    plt.plot(norm_data)
    plt.plot(peaks, norm_data[peaks], "x")
    plt.hlines(*results_half[1:], color="red")
    plt.show()


    index = np.where(norm_data==max_peak)
    print(index)
    plt.plot(lama[int(peak_1)],norm_data[int(peak_1)], "x")
    plt.plot(lama[int(peak_2)],norm_data[int(peak_2)], "x")
    plt.plot(lama, norm_data)
    plt.show()
####### TDCMT CLASS
class tdcmt(object):

    #private
    # fit the spectrum of mode id from th*max to max below fmax with minimum number of poles N
    def fitone(self,id,th,fmax, N = 3, plot=False):
        
        # self.smode = np.load("FT.npy")
        self.smode=self.smode_original # define smode_original otherwise in gives error about boolean, see fitting_andris.py for loop limiting
        # print('len of self.mode is ', len(self.smode))
        
        ff=self.f<=fmax
        MAX=np.max(np.abs(self.smode[ff,id]))
        msk=(np.abs(self.smode[:,id])>=th*MAX)*(self.f<=fmax)#(self.f>=f0)*(self.f<=f1)
        self.w=self.f[msk]

        fun=self.smode[msk,id]

        a=frac_fit()
        a.fit0(self.w,fun,np.real(fun[0]),N)
        self.fted=a.out(self.w)


# lim=np.where(self.f<=fmax)
# self.w=np.zeros(len(ff[lim]))
# print('length of ff is', len(ff))
# print('self.w is',self.w)
# print('size of w is', len(self.w))
# ind_pos=np.where((np.abs(self.smode[:,id])>=th*max)*(self.f<=fmax))
# indexes = ind_pos[0].tolist()
# print('index is', indexes)
# replacements = self.f[indexes]
# for index, replacement in zip(indexes, replacements):
#     self.w[index]=replacement
# print(ind_pos)
# print('new self.w is',self.w)
# print('size of w is', len(self.w))

        lim=np.where(self.f<=fmax)
        self.w_full=self.f[lim]
        # print('length of ff is', len(ff))
        # print('self.w is',self.w)
        # print('size of w is', len(self.w))
        ind_pos=np.where((np.abs(self.smode[:,id])>=th*MAX)*(self.f<=fmax))
        
        indexes = ind_pos[0].tolist()
        # print('index is', indexes)
        #replacements = self.f[indexes]
        self.fted_full=np.zeros(len(ff[lim]), dtype=np.complex)
        for index, replacement in zip(indexes, self.fted):
            self.fted_full[index]=replacement

        # print(ind_pos)
        # print('new self.w is',self.w)
        # print('size of w is', len(self.w))
        #print('fted_full is ', self.fted_full)

        plt.figure(1)
        self.a = a
        #print('self.fted is', self.fted)
        if plot:
            #plt.clf()
            plt.plot(self.f,np.abs(self.smode[:,id])**2,'b',self.w,np.abs(self.fted)**2,'r')
        



        
    
    #public
    def __init__(self,src,ful,dt,tsteps, Nmod):
        # src source tdcmt file, ful full tdcmt file, tsteps fdtd simulation steps
        # unpack data
        self.src=src
        self.ful=ful
        self.dt=dt
        self.tsteps=tsteps
        self.Nmod=Nmod
    

    def fft(self):
        src = self.src
        ful = self.ful
        dt = self.dt
        tsteps=self.tsteps
        Nmod=self.Nmod
                
        with open(self.src, 'rb') as s:
            src_data = np.fromfile(s, dtype=np.float)
            sr = np.reshape(src_data, [tsteps+1, Nmod])
        
        with open(self.ful, 'rb') as c:
            cbd_data = np.fromfile(c, dtype=np.float)
            fu = np.reshape(cbd_data, [tsteps+1, Nmod])
            

        start = time.time()
        pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
        pdss = pyfftw.interfaces.scipy_fftpack.fft(sr, axis=0)
        end = time.time()
        print(end-start)
        
        start = time.time()
        pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
        tmpf = pyfftw.interfaces.scipy_fftpack.fft(fu, axis=0)
        end = time.time()
        print(end-start)
 



        #normalize data
        pdss=np.abs(pdss)**2
        pdss=np.sum(pdss,axis=1) #pds spectrum of the source
        nrm=1./np.sqrt(pdss)
        f=np.fft.fftfreq(tsteps+1,dt)/sc.c #frequencies 
        for i in range(tmpf.shape[1]):
            mask = pdss < 0.02e5
            # mask2 = np.where(f>4.6)
            # mask3 = np.where(f<1)
            mask2 = np.where(f>lama_2)
            mask3 = np.where(f<lama_1)
            tmpf[:,i]=tmpf[:,i]*nrm
            tmpf[mask, i] = 0
            tmpf[mask2, i]=0
            tmpf[mask3, i]=0

        #compute quantities of interest
        dim=int(tsteps/2) #half frequency


        self.f=f[0:dim] #positive frequency
        pdsf=np.abs(tmpf)**2 #pds full
        dos=np.sum(pdsf,axis=1)
        self.dos=dos[0:dim] #dos
        self.pdss=pdss[0:dim] #pds of source
        self.smode=tmpf[0:dim,:] #normalized mode spectrum
        self.smode_original=tmpf[0:dim,:] #define new variable to use in fitone, otherwise gives boolean error
        
        np.save('fft_data.npy',self.smode)



class frac_fit(object):
    #get an odd line
    #returns (-1)n-1*w^2n-2 in column format
    def getline_odd(self,w,n):
        tmp=w.copy()
        for i in range(len(tmp)):
            tmp[i]=((-1)**(n+1))*(tmp[i]**(2*n-1))
        return tmp

    #get an even line
    #returns (-1)n-1*w^2n-2 in column format
    def getline_even(self,w,n):
        tmp=w.copy()
        for i in range(len(w)):
            tmp[i]=((-1)**(n))*(tmp[i]**(2*n))
        return tmp


    # fit complex data with (a_0+sum_n=1^N an(iw)^n)/(1+sum_n=1^N bn(iw)^n)
    # at M frequencies M=len(w)=len(data)
    # with given behavior at w=0 a_0=real(data(0))
    def fit0(self,w,data,data_0,N):
        M2=len(w) #2*number of frequency points
        #assemble linear system
        #RHS
        f=np.zeros((2*M2,1))
        self.a0=data_0
        for i in range(M2):
            f[i]=np.real(data[i])-self.a0
            f[M2+i]=np.imag(data[i])
        #B matrix
        B=np.zeros((2*M2,4*N))
        # assemble B 
        for n in np.arange(1,N+1): #1,...,N
            # assemble top
            eve=self.getline_even(w,n)
            odd=self.getline_odd(w,n)
            B[0:M2,2*n-1]=eve
            B[0:M2,2*N+2*n-1]=-np.real(data)*eve
            B[0:M2,2*N+2*n-2]=np.imag(data)*odd
            # assemble bottom
            B[M2::,2*n-2]=odd
            B[M2::,2*N+2*n-1]=-np.imag(data)*eve
            B[M2::,2*N+2*n-2]=-np.real(data)*odd
        # solve system        
        self.f=f
        self.B=B
        #self.sol=np.linalg.lstsq(B,f,rcond=-1)
        # weighted lsq
        W=np.diag(f.flatten()**2)
        Bw=np.matmul(W,B)
        fw=np.matmul(W,f)
        self.sol=np.linalg.lstsq(Bw,fw,rcond=-1)
        
    # get output of the fit on frequency w
    def out(self,w):
        #construct poly of numerator/denominator
        sol=self.sol[0].flatten()
        dim=int(len(sol)/2) #number of coefficients at numerator or denominator
        a=np.zeros(dim+1) #plus 0 order term
        b=np.zeros(dim+1) #plus 0 order term
        a[0:-1]=sol[dim-1::-1]
        b[0:-1]=sol[-1:dim-1:-1]
        a[-1]=self.a0
        b[-1]=1
        self.a=a
        self.b=b
        # print('sol is', sol)
        # print('plynomial is', np.polyval(a,1j*w)/np.polyval(b,1j*w))
        return np.polyval(a,1j*w)/np.polyval(b,1j*w)









np.set_printoptions(threshold=sys.maxsize)
fmax = 5
N_mod = 2500
i=0.03
#spectr=tdcmt('/Users/erglisa/Desktop/source_ff=0.27/modes.bin', '/Users/erglisa/Desktop/cuboid_ff=0.27_17/modes.bin',5.89664e-12,100000, 2500)
# read_txt(100001, 5.89664e-12 ,50*50 ,'50x50', i, 0.04, '0.4x1.6_0.4_1.6', '400x400', i, 0.5, 17)
# spectr=tdcmt('/Users/erglisa/Desktop/source_ff=0.27/modes.bin', 
# '/Users/erglisa/Desktop/cuboid_ff=0.27_17/modes.bin',5.89664e-12,100000, N_mod)

read_txt(100001, 5.89664e-12 ,50*50 ,'50x50', i, 0.04, '0.4x1.6_0.4_1.6', '400x400', i, 0.5, 1)  
spectr=tdcmt('/Users/erglisa/Desktop/source_ff=0.03/modes.bin', 
'/Users/erglisa/Desktop/cuboid_ff=0.03_1/modes.bin',5.89664e-12,100000, N_mod)

spectr.fft()



nmodes = spectr.smode.shape[1]
# randomly choose 6 modes to plot
np.random.seed(1)

modes = np.random.choice(np.arange(nmodes),size=1, replace=False)
modes.sort()

#plt.figure(figsize=(5,3), dpi = 200)
    # plt.xlim(0, fmax) 
# plt.ylim(0, 1.8)
plt.xlabel('$1/\lambda$')
#plt.title('modes %s' % str(modes))
FTED = np.array([])



i=0
nmodes=2500

COEFF= []
coeff_matrix=np.zeros((20,2))
w_list = []

for i in range(nmodes):
    error=[]
    poles =6
    # th=0.2+i/nmodes*0.4
    th=0.5#-i/nmodes*0.15
    spectr.fitone(i,th,fmax, N=poles, plot=False)

    # #print('nonzero indexes of fted_full are ', np.nonzero(spectr.fted_full))
    limit = np.where(spectr.f<=fmax) 
    spectr.smode=spectr.smode[limit]
    #print('size of smode is', spectr.smode.shape)
    masking = spectr.fted_full !=0

    #error = np.sum((np.abs(spectr.smode[:,i])**2-np.abs(spectr.fted_full)**2)**2)

    error.append(np.sum((np.abs(spectr.smode[masking,i])**2-np.abs(spectr.fted_full[masking])**2)**2))
    error1 = (np.abs(spectr.smode[masking,i])**2-np.abs(spectr.fted_full[masking])**2)**2
            # print('==============\n\n')
            # print('error at ' + str(i) + 'th mode is', error)

    # print('spectr.mode is', spectr.smode[masking,i])
    # print('spectr.mode length  is', len(spectr.smode[:,i]))
    # print('spectr.mode length with masking is', len(spectr.smode[masking,i]))
    # print('ftd_full length is', len(spectr.fted_full[masking]))


    while error[-1]>50:
        poles=poles+1
        print('\n')
        print('poles:', poles, 'th:', th)
        spectr.fitone(i,th,fmax, N=poles, plot=False)
        spectr.smode=spectr.smode[limit]
        #error = np.sum((np.abs(spectr.smode[:,i])**2-np.abs(spectr.fted_full)**2)**2)
        error.append(np.sum((np.abs(spectr.smode[masking,i])**2-np.abs(spectr.fted_full[masking])**2)**2))
        error1 = (np.abs(spectr.smode[masking,i])**2-np.abs(spectr.fted_full[masking])**2)**2
        print('error at ' + str(i) + 'th mode is', error)
        if poles ==10:
            pole_index = error.index(min(error))
            poles= pole_index+6 #6 is initial number of poles
            spectr.fitone(i,th,fmax, N=poles, plot=False)
            spectr.smode=spectr.smode[limit]
                    # print('\n')
                    # print('pole_index is', pole_index)
                    # print('END_poles are', poles)
            break
    print('mode is', i)    




    w_list.append(spectr.w)
    plt.plot(spectr.f[limit],np.abs(spectr.smode[:,i])**2,'b',spectr.w,np.abs(spectr.fted)**2,'r')

    fted_T=spectr.fted_full
    FTED = np.column_stack((FTED,fted_T)) if FTED.size else fted_T
        # print('# of poles:', poles)
        # print('len of a.a is', len(spectr.a.a))
    # coeff_matrix.resize((len(spectr.a.a),2), refcheck=False) ##either of these work
    coeff_matrix = np.resize(coeff_matrix, (len(spectr.a.a),2)) #either of these work
    # print(coeff_matrix)

    coeff_matrix[:,0]=spectr.a.a
    coeff_matrix[:,1]=spectr.a.b
    # print('len of coeff matrix is', len(coeff_matrix))
    COEFF.append(coeff_matrix)
    # print(len(spectr.a.a))
    # print('COEFF is ', COEFF)

    # plt.figure()
    # plt.plot(error1, '.')

    # plt.show()

# print('wlist is', w_list)

with open('coefficients.npz', 'wb') as f:
    np.save(f, COEFF)
with open('w.npz', 'wb') as h:
    np.save(h, w_list)
with open('w_full.npz', 'wb') as g:
    np.save(g, spectr.w_full)



# plt.show()

# # # print('fted_T shape is', fted_T.shape)
# # # print('fted shape is', spectr.fted.shape)



fted_DOS=np.sum(np.abs(FTED)**2, axis=1)
under=np.where(spectr.f<=5)
total_error=fted_DOS+0.2-spectr.dos[under]
print('total_error is', np.sum(total_error**2))


plt.figure(figsize=(18,8))
plt.plot(spectr.w_full,fted_DOS+0.2, '-',color='tab:red', label='fitted')
plt.plot(spectr.f, spectr.dos, 'tab:blue', label='Simulation')
plt.plot(spectr.w_full,total_error, '.', color='tab:green', label='Error')
plt.xlim(0,5)
# plt.ylim(0,2)
plt.legend()
plt.show()

fig3, ax3 = plt.subplots(2, 1)
fig3.suptitle('DOS - simulation and fitted', fontsize=16, weight='bold')   #fontname="Times New Roman"#, style='italic' )
ax3[0].plot(spectr.w_full,fted_DOS+0.2, '-',color='tab:red', label='fitted')
ax3[0].plot(spectr.f, spectr.dos, 'tab:blue', label='Simulation')
ax3[0].set(title='Density of states ', xlabel="1/lambda", ylabel="Amplitude")
ax3[0].legend()
ax3[1].plot(spectr.w_full,total_error, '.', color='tab:green', label='Error')
ax3[1].set(title='Error = (simulation-fitted)', xlabel="1/lambda", ylabel="Error")
fig3.tight_layout()
fig3.subplots_adjust(top=0.86)
ax3[0].set_xlim(0,5)
ax3[1].set_xlim(0,5)

# plt.savefig(self.dirName+'/MOD_FF='+str(self.ff)+'_r='+str(self.r) + '_Nmod='+self.Nmodtext+'_t='+str(self.t-1)+'_'+str(i)+'.png', dpi=500)
plt.show()
plt.close()    




#print (np.load("FT.npy"))

with open('coeff.txt', 'w') as f:
    writer = csv.writer(f)
    writer.writerows(zip(spectr.a.a, spectr.a.b))
f.close()












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

