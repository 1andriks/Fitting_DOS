import numpy as np
#from scipy import linalg as la
#from matplotlib import pyplot as plt
import math
import pickle
#%matplotlib
import time
import os
import sys #system library
#import pso
    

################ CMD PARAMS
resdir=sys.argv[1] #res directory

################ GENERAL PARAMETERS
wx=0.3
wz=0.3
dx=0.001 #resolution size

################ OBJECTIVE FUNCTION
def myfun(x):

    # x=[wx,wy]
    dz=x[0]

    ##### SET UP dimx and dimy and dimz
    dimx=round(wx/dx)+1
    dimz=round(wz/dx)+1
    upml=20

    ##### SET UP AN EMPTY GEOMETRY FILE    
    f=open('geometry.in', 'w')
    f.write('#air\n@ 0 1.\n#background\nsphere 0 0. 0. 10.\n')
    f.close()

    #### SET NUMBER OF MODES
    nmodes=np.zeros(2,dtype=int)
    nmodes[0]=1
    nmodes[1]=6
    np.savetxt('nmodes.txt',nmodes,fmt='%d')

    # ##### SET UP TFSF 
    a=np.zeros(4)
    a[0]=0.
    a[1]=wx-dx/2.
    a[2]=0.05 #upml*dx+20*dx        #z0
    a[3]=wz - a[2]
    #np.savetxt('tfsf.txt',a)        

    #### DEFINE MY CUBOID
    cbd=np.zeros(4)
    cbd[0]=0.
    cbd[1]=wx-dx/2.
    cbd[2]=wz/2.-dz/2.
    cbd[3]=wz/2.+dz/2.
    np.savetxt('cbd.txt',cbd)        

    ##### SIMULATE
    strc=('mpiexec -n 24 ./nano input.in tfsf_box %f %f %f %f total_steps 50000 io_screen 10000 waveform spl+0.25+0 box_size %f %f number_of_cells %i %i upml_size 0 %i erg_io 30000 base_dir '+resdir) % (a[0],a[1],a[2],a[3],wx,wz,dimx,dimz,upml)
    os.system(strc)
    os.rename(resdir+'/td.bin',resdir+'/src.bin')

    ##### CREATE THE GEOMETRY
    f=open('geometry.in', 'a')
    # material
    f.write('#silicon at 750nm\n@ 2 13.91\n')
    #set up cuboid
    f.write('cuboid 2 %f %f %f %f\n' % (cbd[0],cbd[2],cbd[1],cbd[3]))
    f.close()                

    ####### SIMULATE
    os.system(strc)       
    os.rename(resdir+'/td.bin',resdir+'/ful.bin')


############### BOUNDARY

def resolution(x):
    for i in np.arange(0,x.size):
        dim=np.round(x[i]/dx)
        x[i]=dx*dim

x=np.array([0.1]);
resolution(x) #tranform x in multiple of the resolution
b=myfun(x)


