import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy import fftpack

dt=3.93109e-12
steps=35001
modes=4
data=np.fromfile('/Users/rhlin/desktop/res-good-1d/td.bin', dtype=np.float64)
data=np.reshape(data,(steps,modes))
x=np.arange(0,steps,1)
plt.figure()
plt.plot(x,data[:,0],x,data[:,1],x,data[:,2],x,data[:,3])
plt.xlim(0,2000)
plt.show()


freq_sum=np.zeros((steps-1)/2)
freq_sum_geom=np.zeros((steps-1)/2)

for x in range(0,modes):

    sig=data[:,x]
    sample_freq=fftpack.fftfreq(sig.size, d=dt)
    sig_fft=fftpack.fft(sig)
    pidxs=np.where(sample_freq>0)
    freqs=sample_freq[pidxs]
    power=np.abs(sig_fft)[pidxs]
    freq_sum=freq_sum+power**2
    

fig=plt.figure()
plt.plot(freqs/3e8,freq_sum)
fig.suptitle('spectrum of source')
plt.xlim(0,4)

plt.show()

############################process results with geometry
data_geom=np.fromfile('/Users/rhlin/desktop/res/td.bin', dtype=np.float64)
data_geom=np.reshape(data_geom,(steps,modes))

for x in range(0,modes):
    
    sig=data_geom[:,x]
    sample_freq=fftpack.fftfreq(sig.size, d=dt)
    sig_fft=fftpack.fft(sig)
    pidxs=np.where(sample_freq>0)
    freqs=sample_freq[pidxs]
    power_geom=np.abs(sig_fft)[pidxs]
    freq_sum_geom=freq_sum_geom+power_geom**2

fig2=plt.figure()
plt.plot(freqs/3e8,freq_sum_geom)
fig2.suptitle('spectrum of response')
plt.xlim(0,4)
plt.show()

norm_DOS=freq_sum_geom/freq_sum
fig3=plt.figure()
plt.plot(freqs/3e8,norm_DOS)
fig3.suptitle('normalized DOS')
plt.xlim(0,4)
plt.ylim(0,7.5)
#plt.show()

############################

dt=4.698511e-12
steps=30001
modes=4

data=np.fromfile('/Users/rhlin/desktop/res2/src.bin', dtype=np.float64)
data=np.reshape(data,(steps,modes))
x=np.arange(0,steps,1)

freq_sum_src2=np.zeros((steps-1)/2)

for x in range(0,modes):
    
    sig=data[:,x]
    sample_freq=fftpack.fftfreq(sig.size, d=dt)
    sig_fft=fftpack.fft(sig)
    pidxs=np.where(sample_freq>0)
    freqs_src=sample_freq[pidxs]
    power=np.abs(sig_fft)[pidxs]
    freq_sum_src2=freq_sum_src2+power**2

data=np.fromfile('/Users/rhlin/desktop/res2/ful.bin', dtype=np.float64)
data=np.reshape(data,(steps,modes))
x=np.arange(0,steps,1)

freq_sum_geom2=np.zeros((steps-1)/2)

for x in range(0,modes):
    
    sig=data[:,x]
    sample_freq=fftpack.fftfreq(sig.size, d=dt)
    sig_fft=fftpack.fft(sig)
    pidxs=np.where(sample_freq>0)
    freqs2=sample_freq[pidxs]
    power=np.abs(sig_fft)[pidxs]
    freq_sum_geom2=freq_sum_geom2+power**2


fig4=plt.figure()
plt.plot(freqs/3e8,norm_DOS, freqs2/3e8,freq_sum_geom2/freq_sum_src2)
plt.scatter(1,1.9)
fig4.suptitle('DOS')
plt.xlim(0,4)
plt.ylim(0,10)

plt.show()



