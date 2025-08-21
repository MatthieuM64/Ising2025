#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import sys
import time
import multiprocessing
import gc

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

#Physical parameters
beta=1.
h=0
LX=512
LY=512
init=0
ran=0

#Time intervals
DT=25
tmax=100000

#Movie creation
dpi=180
NCPU=4
multi=True
movie=False

#Read parameters in command line
for arg in sys.argv[1:]:
	if "-beta=" in arg:
		beta=float(arg[6:])
	elif "-h=" in arg:
		h=float(arg[3:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-ran=" in arg:
		ran=int(arg[5:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie" in arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)

#Multi-threading		
if NCPU==1:
	multi=False
elif NCPU>1:
	multi=True
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

#Global parameters
Nsites=LX*LY
nbytes=(Nsites+7)//8
colors_map1=[(1,0,0),(0,0,1)]
cmap1=LinearSegmentedColormap.from_list('my_list1', colors_map1, N=256)

#Create one snapshot (one frame of the movie)
def Snapshot(i):
	if not os.path.isfile('snapshots/figure_Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,h,LX,LY,init,ran,i)):
		t=i*DT

		fig=plt.figure(figsize=(6,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.11,right=0.96,bottom=0.06,top=0.91,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		
		#Import data from binary file
		data=np.fromfile('data_Ising_dynamics2d/Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.bin'%(beta,h,LX,LY,init,ran,t),dtype=np.uint8,count=nbytes)
		STATE=np.unpackbits(data).reshape(LY,LX)
		mag=1-2*np.mean(STATE)
		
		x=np.linspace(0,LX,LX)
		y=np.linspace(0,LY,LY)
		plt.pcolormesh(x,y,STATE,vmin=0,vmax=1,rasterized=True,cmap=cmap1)
		
		#plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0,1.03*LY,'$t=%d$'%(t),ha='left',va='center',fontsize=20)
		plt.text(LX,1.03*LY,'$\\beta=%.8g$, $h=%.8g$, $m=%.4f$'%(beta,h,mag),ha='right',va='center',fontsize=20)
		
		plt.savefig('snapshots/figure_Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(beta,h,LX,LY,init,ran,i),dpi=dpi)
		plt.close()
		
		print('-snap=%d/%d -t=%d -mag=%.8g -tcpu=%d'%(i+1,Nsnap,t,mag,time.time()-clock))
		del fig,STATE,data

#Find the datafile and create the corresponding snapshots
os.system('mkdir -p snapshots/')

i=0
ARG=[]
while os.path.isfile('data_Ising_dynamics2d/Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.bin'%(beta,h,LX,LY,init,ran,i*DT)) and i*DT<=tmax:
	ARG.append(i)
	i+=1	
Nsnap=len(ARG)
	
print('%d Snapshots'%len(ARG))
if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

#Create the movie
if movie:
	os.system('mkdir -p movies')	
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -crf 30 -s %dx%d movies/movie_Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(beta,h,LX,LY,init,ran,6*dpi,6*dpi,beta,h,LX,LY,init,ran))
	if Nsnap>tmax/DT:
		os.system('rm snapshots/figure_Ising_state_beta=%.8g_h=%.8g_LX=%d_LY=%d_init=%d_ran=%d_*.png'%(beta,h,LX,LY,init,ran))

print('OK - time=%d sec'%(time.time()-clock))
