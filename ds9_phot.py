install_path='/opt/local/bin/'
import sys
sys.path.append(install_path)
import numpy as np
#from astropy.io import fits
#from astropy.modeling import models,fitting
import os
from scipy.optimize import curve_fit,minimize
import pyds9
import matplotlib as mpl
mpl.use('tkagg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Tkinter as tk
import math
from photutils import CircularAperture,CircularAnnulus,aperture_photometry
from photutils.background import Background2D,MMMBackground

#to do:
#1: Add interactive way to change radial profile display, perhaps via menu: need to add Moffat profile, that is default for iraf imexam -r. need to be able to compare.
#3: Once centroid is determined, add interactive way to do aperture photometry, perhaps with photutils: http://photutils.readthedocs.io/en/stable/

#package management:
#in /opt/local/bin, installed matplotlib, numpy, scipy, pyfits, pyds9
#tested with most recent versions of all packages, and works

MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5

def zscale (image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None):
    """Implement IRAF zscale algorithm
    nsamples=1000 and contrast=0.25 are the IRAF display task defaults
    bpmask and zmask not implemented yet
    image is a 2-d numpy array
    returns (z1, z2)
    """

    # Sample the image
    samples = zsc_sample (image, nsamples, bpmask, zmask)
    npix = len(samples)
    samples.sort()
    zmin = samples[0]
    zmax = samples[-1]
    # For a zero-indexed array
    center_pixel = (npix - 1) / 2
    if npix%2 == 1:
        median = samples[center_pixel]
    else:
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1])

    #
    # Fit a line to the sorted array of samples
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    ngrow = max (1, int (npix * 0.01))
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, KREJ, ngrow,
                                             MAX_ITERATIONS)

    if ngoodpix < minpix:
        z1 = zmin
        z2 = zmax
    else:
        if contrast > 0: zslope = zslope / contrast
        z1 = max (zmin, median - (center_pixel - 1) * zslope)
        z2 = min (zmax, median + (npix - center_pixel) * zslope)
    return z1, z2

def zsc_sample (image, maxpix, bpmask=None, zmask=None):
    
    # Figure out which pixels to use for the zscale algorithm
    # Returns the 1-d array samples
    # Don't worry about the bad pixel mask or zmask for the moment
    # Sample in a square grid, and return the first maxpix in the sample
    nc = image.shape[0]
    nl = image.shape[1]
    stride = max (1.0, math.sqrt((nc - 1) * (nl - 1) / float(maxpix)))
    stride = int (stride)
    samples = image[::stride,::stride].flatten()
    return samples[:maxpix]
    
def zsc_fit_line (samples, npix, krej, ngrow, maxiter):

    #
    # First re-map indices from -1.0 to 1.0
    xscale = 2.0 / (npix - 1)
    xnorm = np.arange(npix)
    xnorm = xnorm * xscale - 1.0

    ngoodpix = npix
    minpix = max (MIN_NPIXELS, int (npix*MAX_REJECT))
    last_ngoodpix = npix + 1

    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad
    badpix = np.zeros(npix, dtype="int32")

    #
    #  Iterate

    for niter in range(maxiter):

        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
            break
        
        # Accumulate sums to calculate straight line fit
        goodpixels = np.where(badpix == GOOD_PIXEL)
        sumx = xnorm[goodpixels].sum()
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum()
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum()
        sumy = samples[goodpixels].sum()
        sum = len(goodpixels[0])

        delta = sum * sumxx - sumx * sumx
        # Slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta
        slope = (sum * sumxy - sumx * sumy) / delta
        
        # Subtract fitted line from the data array
        fitted = xnorm*slope + intercept
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        ngoodpix, mean, sigma = zsc_compute_sigma (flat, badpix, npix)

        threshold = sigma * krej

        # Detect and reject pixels further than k*sigma from the fitted line
        lcut = -threshold
        hcut = threshold
        below = np.where(flat < lcut)
        above = np.where(flat > hcut)

        badpix[below] = BAD_PIXEL
        badpix[above] = BAD_PIXEL
        
        # Convolve with a kernel of length ngrow
        kernel = np.ones(ngrow,dtype="int32")
        badpix = np.convolve(badpix, kernel, mode='same')

        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0])
        
        niter += 1

    # Transform the line coefficients back to the X range [0:npix-1]
    zstart = intercept - slope
    zslope = slope * xscale

    return ngoodpix, zstart, zslope

def zsc_compute_sigma (flat, badpix, npix):

    # Compute the rms deviation from the mean of a flattened array.
    # Ignore rejected pixels

    # Accumulate sum and sum of squares
    goodpixels = np.where(badpix == GOOD_PIXEL)
    sumz = flat[goodpixels].sum()
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum()
    ngoodpix = len(goodpixels[0])
    if ngoodpix == 0:
        mean = None
        sigma = None
    elif ngoodpix == 1:
        mean = sumz
        sigma = None
    else:
        mean = sumz / ngoodpix
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1))
        if temp < 0:
            sigma = 0.0
        else:
            sigma = math.sqrt (temp)

    return ngoodpix, mean, sigma

'''end zscale routines
# The above functions are available at the url
# https://github.com/spacetelescope/stsci.numdisplay/blob/master/lib/stsci/numdisplay/zscale.py
# and are subject to the following license:
#
# Copyright (C) 2005 Association of Universities for Research in Astronomy (AURA)

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
'''
###############################################################################
###############################################################################

#try new implementation
def gauss1d(x, amplitude, xo, sigma_x, offset):
    return amplitude*np.exp(-((x-xo)**2.)/(2.*sigma_x**2.))+offset

def gauss2d((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def moffat1d(x,amplitude,xo,alpha,beta,offset):
    #fwhm is 2*gamma
    return (amplitude*(1.+((x-xo)/alpha)**2.)**-beta)+offset

def moffat2d((x,y),amplitude,xo,yo,gamma_x,gamma_y,theta,alpha,offset):
    #amplitude,xo,yo,gamma_x,gamma_y,theta,alpha,offset=float(amplitude),float(xo),float(yo),float(gamma_x),float(gamma_y),float(theta),float(alpha),float(offset)
    g=amplitude*(1.+(((x-xo)/gamma_x)**2.+((y-yo)/gamma_y)**2.)**-alpha)+offset
    return g.ravel()

###############################################################################
###############################################################################

# def radial_profile(data, center):
#     y,x = np.indices((data.shape)) # first determine radii of all pixels
#     r = np.sqrt((x-center[0])**2+(y-center[1])**2)
#     ind = np.argsort(r.flat) # get sorted indices
#     sr = r.flat[ind] # sorted radii
#     sim = data.flat[ind] # image values sorted by radii
#     ri = sr.astype(np.int32) # integer part of radii (bin size = 1)
#     # determining distance between changes
#     deltar = ri[1:] - ri[:-1] # assume all radii represented
#     rind = np.where(deltar)[0] # location of changed radius
#     nr = rind[1:] - rind[:-1] # number in radius bin
#     csim = np.cumsum(sim, dtype=np.float64) # cumulative sum to figure out sums for each radii bin
#     tbin = csim[rind[1:]] - csim[rind[:-1]] # sum for image values in radius bins
#     radialprofile = tbin/nr # the answer
#     return radialprofile


def radial_profile(data, center):
    y,x = np.indices((data.shape)) # first determine radii of all pixels
    r = np.sqrt((x-center[0])**2+(y-center[1])**2)
    ind = np.argsort(r.flat) # get sorted indices
    sr = r.flat[ind] # sorted radii
    sim = data.flat[ind] # image values sorted by radii
    return sr,sim

def radial_profile_gauss1d_chisquare((x,y,amplitude,sigma_x,sky),image):
    # x=op_vals[0]
    # y=op_vals[1]
    # amplitude=op_vals[2]
    # if amplitude<0 or sigma_x<0 or sky<0:return 50
    radial=radial_profile(image,(x,y))
    # print 0.5*np.sum((radial[1]-gauss1d(radial[0],amplitude,0,sigma_x,sky))**2.)
    return 0.5*np.sum((radial[1]-gauss1d(radial[0],amplitude,0,sigma_x,sky))**2.)
    
    

###############################################################################
###############################################################################

#def onmove(event):
        #if event.inaxes:
		#just update the labeltext
		#labeltext.set('X: '+str(int(event.ydata)+x-new_image_halfwidth))
		#label2text.set('Y: '+str(int(event.xdata)+y-new_image_halfwidth))
		#label3text.set('Value: '+str(d.get('data physical '+str(int(event.ydata)+x-new_image_halfwidth)+' '+str(int(event.xdata)+y-new_image_halfwidth)+' 1 1 yes')))

def onpress(event):
        if event.key=='p':
            	filename=d.get('file')
		#hdu=fits.open(filename)
		hdu=d.get_pyfits()
		data=hdu[0].data
		x=d.get('crosshair image')
		x,y=x.split()
		x,y=int(float(x)),int(float(y))

		new_image_halfwidth=int(entry1.get())
		aperture_radius=int(entry2.get())
		inner_sky_radius=int(entry3.get())
		outer_sky_radius=int(entry4.get())

		new_image=data[y-new_image_halfwidth:y+new_image_halfwidth,x-new_image_halfwidth:x+new_image_halfwidth]

		x_grid=np.arange(x-new_image_halfwidth,x+new_image_halfwidth,1)
		y_grid=np.arange(y-new_image_halfwidth,y+new_image_halfwidth,1)
		x_grid, y_grid = np.meshgrid(x_grid, y_grid)
		guess=(np.amax(new_image)-new_image[0,0],x,y,3,3,0,new_image[0,0])
		popt,pcov=curve_fit(gauss2d,(x_grid,y_grid),new_image.ravel(),p0=guess)
		#print popt
		x=int(popt[1])
		y=int(popt[2])
                labeltext.set('X: '+str(x))
                label2text.set('Y: '+str(y))


		new_image=data[y-new_image_halfwidth:y+new_image_halfwidth,x-new_image_halfwidth:x+new_image_halfwidth]
        	# for artist in ax1.get_children():
            	# 	if hasattr(artist,'get_label') and artist.get_label()=='centroid':
                # 		artist.remove()
		ax1.clear()
        	ax1.matshow(np.flip(new_image,axis=0),cmap='gray',origin='upper',clim=zscale(new_image),zorder=0)
    		ax1.scatter([new_image_halfwidth+1],[new_image_halfwidth-1],marker='+',s=120,c='k',zorder=1)#ax1.scatter([popt[1]-x+new_image_halfwidth],[popt[2]-y+new_image_halfwidth],marker='+',s=120,c='k',zorder=1)
		aperture_circle=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=aperture_radius,linewidth=3,color='hotpink',fill=False,lw=3,zorder=2)
		ax1.add_patch(aperture_circle)
		inner_sky=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=inner_sky_radius,linewidth=3,color='lime',fill=False,zorder=2)
		ax1.add_patch(inner_sky)
		outer_sky=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=outer_sky_radius,linewidth=3,color='red',fill=False,zorder=2)
		ax1.add_patch(outer_sky)
		canvas.draw()

		#update the radial plot
		ax2.clear()
		radial=radial_profile(new_image,[new_image_halfwidth,new_image_halfwidth])
		#perform the aperture photometry
		#currently 2% different from IDL atv
		aperture=CircularAperture((x,y),r=aperture_radius)
		annulus=CircularAnnulus((x,y),r_in=inner_sky_radius,r_out=outer_sky_radius)
		phot_table=aperture_photometry(data,[aperture,annulus])
		#new background estimation
		bkg_mask=np.ma.masked_outside(radial[0].reshape(new_image.shape),inner_sky_radius,outer_sky_radius)
		bkg_map=Background2D(new_image,tuple(np.array(new_image.shape)/4),mask=bkg_mask.mask,exclude_mesh_method='all')
		bkg_map_med=MMMBackground().calc_background(bkg_map.data)#bkg_map_med=np.median(bkg_map.background)
		#print 'Map sky mean '+str(bkg_map_med)
		#bkg_mean=phot_table['aperture_sum_1']/annulus.area()
		#print 'Aperture sky mean '+str(bkg_mean)
		#phot_table['residual_aperture_sum']=phot_table['aperture_sum_0']-bkg_mean*aperture.area()
		phot_table['residual_aperture_sum']=phot_table['aperture_sum_0']-bkg_map_med*aperture.area()
		#print 'Map sky result: '+str(phot_table['aperture_sum_0']-bkg_map_med*aperture.area())
		#print "Aperture Photometry Result: "+str(phot_table['residual_aperture_sum'])
                label8text.set('Sky Value: '+str(int(bkg_map_med)))
		label9text.set('Aperture Counts: '+str(int(phot_table['residual_aperture_sum'][0])))
                label10text.set('Mag: '+str(-2.5*np.log10(int(phot_table['residual_aperture_sum'][0]))+25.)[:5])


                
		ax2.scatter(radial[0],radial[1])
                if var10.get()==1:
			ax2.plot(np.linspace(0,new_image_halfwidth,num=50),gauss1d(np.linspace(0,new_image_halfwidth,num=50),popt[0],0,np.mean([popt[3],popt[4]]),popt[6]),c='k',lw=2)
			ax2.text(0.5,0.93,'Gaussian FWHM: '+str(2.35482*np.mean([popt[3],popt[4]]))[:5],transform=ax2.transAxes,fontsize=int(15*scaling))
	    	if var11.get()==1:
                    	moffat1d_guess=(np.amax(new_image)-bkg_map_med,0,3,1,bkg_map_med)
		    	popt2,pcov2=curve_fit(moffat1d,radial[0],radial[1],p0=moffat1d_guess)
			ax2.plot(np.linspace(0,new_image_halfwidth,num=50),moffat1d(np.linspace(0,new_image_halfwidth,num=50),popt2[0],popt2[1],popt2[2],popt2[3],popt2[4]),c='r',lw=2)
                        ax2.text(0.5,0.85,'Moffat FWHM: '+str(2.0*popt2[2]*np.sqrt(2.0**(1./popt2[3])-1.))[:5],transform=ax2.transAxes,fontsize=int(15*scaling))
		ax2.grid(True,color='white',linestyle='-',linewidth=1)
		ax2.set_axisbelow(True)
		ax2.autoscale(False)
		ax2.set_xlim([0,new_image_halfwidth])
		ax2.set_xlabel('Radius (Pixels)')
		ax2.set_ylim([np.amin(radial[1]),np.amax(radial[1])])
		ax2.set_axis_bgcolor('0.85')
		ax2.axvline(aperture_radius,linewidth=2,color='hotpink')
		ax2.axvline(inner_sky_radius,linewidth=2,color='lime')
		ax2.axvline(outer_sky_radius,linewidth=2,color='red')
                ax2.axhline(bkg_map_med,linewidth=2,color='yellow')
		canvas2.draw()


def handler():
        open_opt=open(install_path+'phot_vars.py','r')
        linedata=open_opt.readlines()
        open_opt.close()
        linedata[0]='new_image_halfwidth_default='+entry1.get()+'\n'
        linedata[1]='aperture_radius='+entry2.get()+'\n'
        linedata[2]='inner_sky_radius='+entry3.get()+'\n'
        linedata[3]='outer_sky_radius='+entry4.get()+'\n'
	new_str=''
        for j in linedata:new_str+=j
        open_opt=open(install_path+'phot_vars.py','w')
        open_opt.write(new_str)
        open_opt.close()

def update():
	MyObject = type('MyObject', (object,), {})
	obj = MyObject()
	obj.key = 'p'
	handler()
	onpress(obj)


####################################main#######################################
#import the user defined display parameters
from phot_vars import *

## tkinter integration
root=tk.Tk()
#set title
root.title("Aperture Photometry")

scaling=0.75

label4text = tk.StringVar()
label4text.set('Half Width of Star Window')
label4 = tk.Label(root, textvariable=label4text, height=1)
label4.config(font=("Verdana", int(16*scaling)))
label4.grid(row=2,column=0,columnspan=1)

#set the width of the entry box column
#root.columnconfigure(2,width=

entry1=tk.Entry(root)
entry1.insert(0,str(new_image_halfwidth_default))
entry1.grid(row=2,column=2)
new_image_halfwidth=int(entry1.get())

label5text = tk.StringVar()
label5text.set('Aperture Radius (pix)')
label5=tk.Label(root, textvariable=label5text, height=1)
label5.config(font=("Verdana", int(16*scaling)))
label5.grid(row=3,column=0,columnspan=1)

entry2=tk.Entry(root)
entry2.insert(0,str(aperture_radius))
entry2.grid(row=3,column=2)
aperture_radius=int(entry2.get())

entry3=tk.Entry(root)
entry3.insert(0,str(inner_sky_radius))
entry3.grid(row=4,column=2)
inner_sky_radius=int(entry3.get())

label6text = tk.StringVar()
label6text.set('Inner Sky Radius (pix)')
label6=tk.Label(root, textvariable=label6text, height=1)
label6.config(font=("Verdana", int(16*scaling)))
label6.grid(row=4,column=0,columnspan=1)

entry4=tk.Entry(root)
entry4.insert(0,str(outer_sky_radius))
entry4.grid(row=5,column=2)
outer_sky_radius=int(entry4.get())

label7text = tk.StringVar()
label7text.set('Outer Sky Radius (pix)')
label7=tk.Label(root, textvariable=label7text, height=1)
label7.config(font=("Verdana", int(16*scaling)))
label7.grid(row=5,column=0,columnspan=1)

#label10text = tk.StringVar()
#label10text.set('Gaussian PSF')
#label10=tk.Label(root,textvariable=label10text,height=1)
#label10.config(font=("Verdana",12))
#label10.grid(row=2,column=3,columnspan=1)
var10 = tk.IntVar()
var10.set(1)
tk.Checkbutton(root, text="Gaussian PDF", font=("Verdana",int(16*scaling)), variable=var10,command=update).grid(row=2, column=3,columnspan=1)

var11 = tk.IntVar()
var11.set(1)
tk.Checkbutton(root, text="Moffat PDF", font=("Verdana",int(16*scaling)), variable=var11,command=update).grid(row=2, column=4,columnspan=1)


b = tk.Button(root, text="Update Parameters",command=update)#width=15
b.grid(row=6,column=1)

d=pyds9.DS9()
filename=d.get('file')
#hdu=fits.open(filename)
hdu=d.get_pyfits()
data=hdu[0].data

x=d.get('crosshair image')
x,y=x.split()
x,y=int(float(x)),int(float(y))


new_image=data[y-new_image_halfwidth:y+new_image_halfwidth,x-new_image_halfwidth:x+new_image_halfwidth]

x_grid=np.arange(x-new_image_halfwidth,x+new_image_halfwidth,1)
y_grid=np.arange(y-new_image_halfwidth,y+new_image_halfwidth,1)
x_grid, y_grid = np.meshgrid(x_grid, y_grid)
#guess=(np.amax(new_image)-new_image[0,0],x,y,3,3,0,new_image[0,0])
guess=(np.amax(new_image)-new_image[0,0],x_grid[new_image==np.amax(new_image)][0],y_grid[new_image==np.amax(new_image)][0],3,3,0,new_image[0,0])
#def gauss2d((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
popt,pcov=curve_fit(gauss2d,(x_grid,y_grid),new_image.ravel(),p0=guess)

#guess=(x,y,np.amax(new_image)-new_image[0,0],3,new_image[0,0])
#def radial_profile_gauss1d_chisquare((x,y,amplitude,sigma_x,sky),image):
#popt,pcov=curve_fit(gauss1d,x_grid,radial_profile(new_image,[new_image_halfwidth,new_image_halfwidth]),p0=guess)
#popt=minimize(radial_profile_gauss1d_chisquare,guess,args=(new_image),options={'maxiter':500},method = 'BFGS',bounds=((0,2*new_image_halfwidth),(0,2*new_image_halfwidth),(0,65534),(0,np.inf),(0,np.inf)))

#print popt
x=int(popt[1])
y=int(popt[2])

new_image=data[y-new_image_halfwidth:y+new_image_halfwidth,x-new_image_halfwidth:x+new_image_halfwidth]

#plot the star in the image, in zscale
fig=Figure(figsize=(int(6*scaling),int(6*scaling)))
ax1=fig.add_axes([0,0,1,1])
ax1.matshow(np.flip(new_image,axis=0),cmap='gray',origin='upper',clim=zscale(new_image),zorder=0)
ax1.autoscale(False)
ax1.scatter([new_image_halfwidth+1],[new_image_halfwidth-1],marker='+',s=120,c='k',zorder=1)#ax1.scatter([popt[1]-x+new_image_halfwidth],[popt[2]-y+new_image_halfwidth],marker='+',s=120,c='k',label='centroid',zorder=1)
aperture_circle=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=aperture_radius,linewidth=3,color='hotpink',fill=False,zorder=2)
ax1.add_patch(aperture_circle)
inner_sky=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=inner_sky_radius,linewidth=3,color='lime',fill=False,zorder=2)
ax1.add_patch(inner_sky)
outer_sky=plt.Circle((popt[1]-x+new_image_halfwidth,popt[2]-y+new_image_halfwidth),radius=outer_sky_radius,linewidth=3,color='red',fill=False,zorder=2)
ax1.add_patch(outer_sky)


radial=radial_profile(new_image,[new_image_halfwidth,new_image_halfwidth])
#perform the aperture photometry
#currently 2% different from IDL atv
aperture=CircularAperture((x,y),r=aperture_radius)
annulus=CircularAnnulus((x,y),r_in=inner_sky_radius,r_out=outer_sky_radius)
phot_table=aperture_photometry(data,[aperture,annulus])
#new background estimation
bkg_mask=np.ma.masked_outside(radial[0].reshape(new_image.shape),inner_sky_radius,outer_sky_radius)
bkg_map=Background2D(new_image,tuple(np.array(new_image.shape)/4),mask=bkg_mask.mask,exclude_mesh_method='all')
bkg_map_med=MMMBackground().calc_background(bkg_map.data)#bkg_map_med=np.median(bkg_map.background)
#print 'Map sky mean '+str(bkg_map_med)
#bkg_mean=phot_table['aperture_sum_1']/annulus.area()
#print 'Aperture sky mean '+str(bkg_mean)
#phot_table['residual_aperture_sum']=phot_table['aperture_sum_0']-bkg_mean*aperture.area()
phot_table['residual_aperture_sum']=phot_table['aperture_sum_0']-bkg_map_med*aperture.area()
#print 'Map sky result: '+str(phot_table['aperture_sum_0']-bkg_map_med*aperture.area())
#print "Aperture Photometry Result: "+str(phot_table['residual_aperture_sum'])
fig2=Figure(figsize=(int(6*scaling),int(6*scaling)))
fig2.set_facecolor('0.85')
ax2=fig2.add_axes([0.1,0.1,.9,.9])
ax2.grid(True,color='white',linestyle='-',linewidth=1)
ax2.set_axisbelow(True)
ax2.scatter(radial[0],radial[1])
ax2.plot(np.linspace(0,new_image_halfwidth,num=50),gauss1d(np.linspace(0,new_image_halfwidth,num=50),popt[0],0,np.mean([popt[3],popt[4]]),popt[6]),c='k',lw=2)

#guess=(np.amax(new_image)-new_image[0,0],x,y,3,3,0,new_image[0,0])
#popt,pcov=curve_fit(gauss2d,(x_grid,y_grid),new_image.ravel(),p0=guess)
#moffat1d(x,amplitude,xo,gamma,alpha,offset):
#def moffat1d(x,amplitude,xo,alpha,beta,offset):
    #fwhm is 2*gamma
    #return (amplitude*(1.+((x-xo)/alpha)**2.)**-beta)+offset
moffat1d_guess=(np.amax(new_image)-bkg_map_med,0,3,1,bkg_map_med)
popt2,pcov2=curve_fit(moffat1d,radial[0],radial[1],p0=moffat1d_guess)
#moffat2d_guess=(np.amax(new_image)-bkg_map_med,x,y,3.0,3.0,10.**-8.,1.,bkg_map_med)
#popt3,pcov3=curve_fit(moffat2d,(x_grid,y_grid),new_image.ravel(),p0=moffat2d_guess)
# plt.scatter(radial[0],radial[1])
# plt.plot(radial[0],moffat1d(radial[0],popt2[0],popt2[1],popt2[2],popt2[3],popt2[4]))
# plt.show()
ax2.plot(np.linspace(0,new_image_halfwidth,num=50),moffat1d(np.linspace(0,new_image_halfwidth,num=50),popt2[0],popt2[1],popt2[2],popt2[3],popt2[4]),c='r',lw=2)
ax2.text(0.5,0.85,'Moffat FWHM: '+str(2.0*popt2[2]*np.sqrt(2.0**(1./popt2[3])-1.))[:5],transform=ax2.transAxes,fontsize=15*scaling)

ax2.autoscale(False)
ax2.set_xlim([0,new_image_halfwidth])
ax2.set_xlabel('Radius (Pixels)')
ax2.set_ylim([np.amin(radial[1]),np.amax(radial[1])])
ax2.text(0.5,0.93,'Gaussian FWHM: '+str(2.35482*np.mean([popt[3],popt[4]]))[:5],transform=ax2.transAxes,fontsize=15*scaling)
ax2.set_axis_bgcolor('0.85')
ax2.axvline(aperture_radius,linewidth=2,color='hotpink')
ax2.axvline(inner_sky_radius,linewidth=2,color='lime')
ax2.axvline(outer_sky_radius,linewidth=2,color='red')
ax2.axhline(bkg_map_med,linewidth=2,color='yellow')

canvas=FigureCanvasTkAgg(fig,master=root)
plot_widget=canvas.get_tk_widget()
plot_widget.grid(row=0,column=0,columnspan=3)

canvas2=FigureCanvasTkAgg(fig2,master=root)
plot_widget=canvas2.get_tk_widget()
plot_widget.grid(row=0,column=3,columnspan=3)

#initialize label object
label0text = tk.StringVar()
#set label value
label0text.set('Centroid: ')
label0 = tk.Label(root, textvariable=label0text, height=1)
label0.config(font=("Verdana", int(16*scaling),'bold'))
label0.grid(row=1,column=0) #adds the label to the program's grid


#initialize label object
labeltext = tk.StringVar()
#set label value
labeltext.set('X: '+str(x))
label1 = tk.Label(root, textvariable=labeltext, height=1)
label1.config(font=("Verdana", int(16*scaling)))
label1.grid(row=1,column=1) #adds the label to the program's grid

#initialize label object
label2text = tk.StringVar()
#set label value
label2text.set('Y: '+str(y))
label2 = tk.Label(root, textvariable=label2text, height=1)
label2.config(font=("Verdana", int(16*scaling)))
label2.grid(row=1,column=2) #adds the label to the program's grid

#initialize label object
# label3text = tk.StringVar()
# label3text.set('Value: '+str(d.get('data physical '+str(x)+' '+str(y)+' 1 1 yes')))
# label3 = tk.Label(root, textvariable=label3text, height=1)
# label3.config(font=("Verdana", 12))
# label3.grid(row=1,column=2)

label8text = tk.StringVar()
label8text.set('Sky Value: '+str(int(bkg_map_med)))
label8=tk.Label(root, textvariable=label8text, height=1)
label8.config(font=("Verdana", int(16*scaling)))
label8.grid(row=1,column=3)


label9text = tk.StringVar()
label9text.set('Aperture Counts: '+str(int(phot_table['residual_aperture_sum'][0])))
label9=tk.Label(root, textvariable=label9text, height=1)
label9.config(font=("Verdana", int(16*scaling)))
label9.grid(row=1,column=4)

label10text = tk.StringVar()
label10text.set('Mag: '+str(-2.5*np.log10(int(phot_table['residual_aperture_sum'][0]))+25.)[:5])
label10=tk.Label(root, textvariable=label10text, height=1)
label10.config(font=("Verdana", int(16*scaling)))
label10.grid(row=1,column=5)

#canvas.mpl_connect('motion_notify_event',onmove)
canvas.mpl_connect('key_press_event',onpress)
#root.protocol('WM_DELETE_WINDOW',handler)
root.mainloop()
