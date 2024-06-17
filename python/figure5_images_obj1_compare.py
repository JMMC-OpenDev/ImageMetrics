import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import astropy.constants as cst
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.ndimage import rotate
from scipy import ndimage
from scipy import interpolate as intp
from image_utils import *
import cv2
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmcrameri.cm as cmc
from fnmatch import fnmatch

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset




#plot=1 # Drevon Obj1 1st contribution
#plot=2 # Drevon Obj1 2nd contribution
plot=3 # Norris Obj1 contribution
#plot=4 # Model images Obj1 as test
#plot=5 # Drevon Obj2 1st contribution
#plot=6 # Drevon Obj2 2nd contribution
#plot=7 # Norris Obj2 contribution
# plot=8 # Model images Obj2 as test

contribs = ("Drevon_Obj1_1st", "Drevon_Obj1_2nd", "Norris_Obj1", "Model_Obj1", "Drevon_Obj2_1st", "Drevon_Obj2_2nd", "Norris_Obj2", "Model_Obj2")

#########################################################
# Drevon contribution Obj1
if plot ==1:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Drevon/SPIE_contest_final/OBJ1/"
    image_test = ("POLY_OBJ1_PION_cube.fits",\
    "POLY_OBJ1_GRAV_cube.fits",\
    "POLY_OBJ1_MATLM_cube.fits",\
    "POLY_OBJ1_MATN_cube.fits")

#########################################################
# New Drevon contribution Obj1
if plot ==2:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Drevon/SPIE-NEW/OBJ1/"
    image_test = ("POLY_OBJ1_PION_cube.fits",\
    "POLY_OBJ1_GRAV_cube.fits",\
    "POLY_OBJ1_MATLM_cube.fits",\
    "POLY_OBJ1_MATN_cube.fits")

#########################################################
# Norris contribution Obj1
if plot == 3:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Norris/"
    image_test = ("obj1_pionier_norris.fits",\
    "obj1_gravity_norris.fits",\
    "obj1_matisse29_norris.fits",\
    "obj1_matisse8_norris.fits")

#########################################################
# Model images Obj1 as test
if plot == 4:
    dir_test = "/home/fmillour/2024_IBC/Contest_DATA/Obj1/"
    image_test = ("output_Spiral_1710269539282_H.fits",\
    "output_Spiral_1710273340532_K.fits",\
    "output_Spiral_1710273706182_L.fits",\
    "output_Spiral_1710274957235_N.fits")

#########################################################
# Drevon contribution Obj2
if plot == 5:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Drevon/SPIE_contest_final/OBJ2/"
    image_test = ("POLY_OBJ2_PION_cube.fits",\
    "POLY_OBJ2_GRAV_cube.fits",\
    "POLY_OBJ2_MATLM_cube.fits",\
    "POLY_OBJ2_MATN_cube.fits")

#########################################################
# New Drevon contribution Obj2
if plot == 6:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Drevon/SPIE-NEW/OBJ2/"
    image_test = ("POLY_OBJ2_PION_cube.fits",\
    "POLY_OBJ2_GRAV_cube.fits",\
    "POLY_OBJ2_MATLM_cube.fits",\
    "POLY_OBJ2_MATN_cube.fits")
    
#########################################################
# Norris contribution Obj2
if plot == 7:
    dir_test = "/home/fmillour/2024_IBC/RESULTS/Norris/"
    image_test = ("obj2_pionier_norris.fits",\
    "obj2_gravity_norris.fits",\
    "obj2_matisse3_norris.fits",\
    "obj2_matisse8_norris.fits")

#########################################################
# Model images Obj2 as test
if plot == 8:
    dir_test = "/home/fmillour/2024_IBC/Contest_DATA/Obj2/"
    image_test = ("Model_H.fits",\
    "Model_K.fits",\
    "Model_L.fits",\
    "Model_N.fits")
    
    # Data oifits Obj2
    oifits_test = ("Obj2_PIONIER_1.5-1.8.fits",
                    "Obj2_GRAVITY_2.0-2.5.fits",
                    "Obj2_MATISSE_2.9-4.0.fits",
                    "Obj2_MATISSE_8-13.fits")

#########################################################
#########################################################
# Model images Obj1
if plot < 5:
    dir_model = "/home/fmillour/2024_IBC/Contest_DATA/Obj1/"
    image_model = ("output_Spiral_1710269539282_H.fits",\
    "output_Spiral_1710273340532_K.fits",\
    "output_Spiral_1710273706182_L.fits",\
    "output_Spiral_1710274957235_N.fits")

    # Data oifits Obj1
    oifits_model = ("Obj1_PIONIER_1.5-1.8.fits",
                    "Obj1_GRAVITY_2.0-2.5.fits",
                    "Obj1_MATISSE_2.9-4.2.fits",
                    "Obj1_MATISSE_8-13.fits")

#########################################################
# Model images Obj2 
else:
    dir_model = "/home/fmillour/2024_IBC/Contest_DATA/Obj2/"
    image_model = ("Model_H.fits",\
    "Model_K.fits",\
    "Model_L.fits",\
    "Model_N.fits")

    oifits_model = ("Obj2_PIONIER_1.5-1.8.fits",
                    "Obj2_GRAVITY_2.0-2.5.fits",
                    "Obj2_MATISSE_2.9-4.0.fits",
                    "Obj2_MATISSE_8-13.fits")

#########################################################
#########################################################

Bmax = 130;
xunit = 'mas'
lunit='micron'
padfact = 0;

#####################################################
# Set figure
my_dpi=300;

nx = len(image_test)
ny = 4
fig1, ax1b = plt.subplots(ncols=nx, nrows=1, figsize=(10,10/3)) # Notice the equal aspect ratio
ax1 = ax1b.flatten()
fig2, ax2b = plt.subplots(ncols=nx, nrows=1, figsize=(10,10/3)) # Notice the equal aspect ratio
ax2 = ax2b.flatten()
fig3, ax3b = plt.subplots(ncols=nx, nrows=1, figsize=(10,10/3)) # Notice the equal aspect ratio
ax3 = ax3b.flatten()
fig4, ax4b = plt.subplots(ncols=nx, nrows=1, figsize=(10,10/3)) # Notice the equal aspect ratio
ax4 = ax4b.flatten()

for ida, a in enumerate(ax3):
    if ida < 3*nx:
        a.set_aspect('equal')
fig1.subplots_adjust(wspace=0.02, hspace=0.02)

count = 0

for itx,iimg in enumerate(image_test):
    if fnmatch(iimg, 'output_Spiral*') or plot <= 4:
        rng = 45
    else:
        rng = 249
            
    #########################################
    # Load test image
    #print("loading ",iimg)
    # Load test image
    tst_img, tst_px, tst_py, tst_wlen = load_img(dir_test+iimg)

    # Revamp wavelength table for Norris images
    if plot==3 or plot==7:
        # Get wavelengths from the oifits file
        print('Wavelength table screwed up in Norris contribution! Get wavelengths directly from the oifits file')
        print('old wlen',tst_wlen)
        fh = fits.open(dir_model+oifits_model[itx])
        oi_wlen = fh['OI_WAVELENGTH'].data['EFF_WAVE']*1e6
        tst_wlen = np.linspace(oi_wlen[0], oi_wlen[-1], np.shape(tst_img)[0])
        print('new wlen',tst_wlen)
    
    # Rotate and scale image: for the case of the spiral input
    if fnmatch(image_test[itx], 'output_Spiral*'):
        print('Spiral model: rotating image 73 degrees')
        rotangl = -73;
        for wdx,w in enumerate(tst_wlen):
            tst_img[wdx,:,:] = rotate(tst_img[wdx,:,:], rotangl, reshape=False, mode='constant',order=2)
        print('size image ', np.shape(tst_img))

        # Scale image to the one used to produce visibilities
        scaleFact = 10;
        tst_px/=scaleFact
        tst_py/=scaleFact
    elif fnmatch(image_test[itx], 'Model_*'):
        print('Disk model: rotating image 23 degrees')
        rotangl = 23;
        for wdx,w in enumerate(tst_wlen):
            tst_img[wdx,:,:] = rotate(tst_img[wdx,:,:], rotangl, reshape=False, mode='constant',order=2)
    
    
    print('test image size', np.shape(tst_img))

    # Zero pad image to avoid issue when recentering it
    tst_pad_im, tst_padx, tst_pady = pad_img(tst_img, padfact+1, tst_px, tst_py)
    tst_nux = np.shape(tst_pad_im)[1]

    tnbwlen = np.shape(tst_wlen)[0]
    tnpy = tst_nux
    
    # Calculate beam size
    tst_pixsize = np.abs(np.diff(tst_padx)[0])
    tst_beam2, tst_beampix2 = calc_beamsize(tst_wlen, Bmax, tst_pixsize)

    # Blur image to the beam size
    tst_blurred_im3 = blur_imgCube(tst_pad_im, tst_beampix2)

    # recenter image cube
    #tst_cen_im = recenter_imgcube(tst_blurred_im3)
    tst_cen_im = tst_blurred_im3

    #########################################
    # Load model image
    mod_img, mod_px, mod_py, mod_wlen = load_img(dir_model+image_model[itx])
    
    # Renormalize image flux to 1 for each wavelength
    print('Renormalize model image flux to 1 per wavelength...')
    for wdx, w in enumerate(mod_wlen):
        mod_img[wdx,:,:] /= np.sum(mod_img[wdx,:,:])
    
    print('model image size', np.shape(mod_img))
    
    # Rotate image: for the case of the spiral input
    if fnmatch(image_model[itx], 'output_Spiral*'):
        print('Spiral model: rotating image 73 degrees')
        rotangl = -73;
        for wdx,w in enumerate(mod_wlen):
            mod_img[wdx,:,:] = rotate(mod_img[wdx,:,:], rotangl, reshape=False, mode='constant',order=2)
        print('size image ', np.shape(mod_img))

        # Scale image to the one used to produce visibilities
        scaleFact = 10;
        mod_px/=scaleFact
        mod_py/=scaleFact
    if fnmatch(image_model[itx], 'Model_*'):
        print('Disk model: rotating image 23 degrees')
        rotangl = 23;
        for wdx,w in enumerate(mod_wlen):
            mod_img[wdx,:,:] = rotate(mod_img[wdx,:,:], rotangl, reshape=False, mode='constant',order=2)

    # Zero pad image to avoid issue when recentering it
    if plot <= 4:
        mod_pad_im, mod_padx, mod_pady = pad_img(mod_img, padfact, mod_px, mod_py)
    else:
        mod_pad_im, mod_padx, mod_pady = pad_img(mod_img, padfact-1, mod_px, mod_py)
    mod_nux = np.shape(mod_pad_im)[1]

    mnbwlen = np.shape(mod_wlen)[0]
    mnpy = mod_nux
    
    # Calculate beam size
    mod_pixsize = np.abs(np.diff(mod_padx)[0])
    mod_beam2, mod_beampix2 = calc_beamsize(mod_wlen, Bmax, mod_pixsize)

    # Blur image to the beam size
    mod_blurred_im3 = blur_imgCube(mod_pad_im, mod_beampix2)

    # recenter image cube
    #mod_cen_im = recenter_imgcube(mod_blurred_im3)
    mod_cen_im = mod_blurred_im3

    if(tnpy > mnpy):
        newpy = tst_pady
        newpx = tst_padx
    else:
        newpy = mod_pady
        newpx = mod_padx
    if(tnbwlen > mnbwlen):
        newwlen = tst_wlen
    else:
        newwlen = mod_wlen
        
    ############################################################
    # Resample both images to the largest image grid
    print('Resampling both images to the largest image grid...')
    print(mod_pady, mod_padx, mod_wlen, newpy, newpx, newwlen)
    rmod_im = resample_img_xy(mod_cen_im, mod_padx, mod_pady, newpx, newpy)
    rtst_im = resample_img_xy(tst_cen_im, tst_padx, tst_pady, newpx, newpy)
    
    newimg = np.zeros_like(rtst_im)
    newmod = np.copy(newimg)
    distance = np.zeros_like(tst_wlen)
    for iwl, w, in enumerate(tst_wlen):
        print('_______________________________________________________')
        print('Resampling model image to the test image wavelengths...')
        print('wlen', w)
        fi = intp.RegularGridInterpolator((mod_wlen, newpx, newpy), rmod_im, bounds_error=False, fill_value=None)
        nW, nX,nY = np.meshgrid(tst_wlen[iwl], newpx, newpy, indexing='ij')
        rmod_intp = fi((nW, nX, nY))[0,:,:]
        print(np.shape(rmod_intp))
        
        newimg[iwl,:,:], distance[iwl] = compar_images(rmod_intp, rtst_im[iwl,:,:])
        newmod[iwl,:,:] = rmod_intp
        
        # fig, ax = plt.subplots(1,3)
        # ax[0].imshow(rmod_intp)
        # ax[1].imshow(newimg[iwl,:,:])
        # ax[2].imshow(rmod_intp - newimg[iwl,:,:])
        # plt.show()
    
    tst_med = np.median(newimg,axis=0)
    mod_med = np.median(newmod,axis=0)
    
    ####################################################################
    # Calculate the score according to Thiebaut et al.!
    score = 1 - distance / np.sum(newmod, axis=(1,2))
    scoreavg = np.mean(score)
    scorestd = np.std(score)
    
    difmap = newmod - newimg
    difmed = np.median(difmap, axis=0)
    rmsmap = np.sqrt(np.mean((newmod - newimg)**2, axis=0))
    
    mxim = np.max(mod_med)**(1./3.)
    mxlim = np.max(difmed)
    
    ax1[count].imshow(mod_med**(1./3.), origin='lower', extent=[newpx[0],newpx[-1],newpy[0],newpy[-1]], cmap=cmc.batlowW_r, vmin=0, vmax=mxim)
    #hide y-axis
    if(count != 0):
        for k in np.arange(ny):
            ax1[count].get_yaxis().set_visible(False)
    ax1[count].set_xlim(rng, -rng)
    ax1[count].set_ylim(-rng, rng)
    ax1[count].set_xlabel(r'$\alpha$RA [mas]')
    ax1[count].set_ylabel(r'$\delta$DEC [mas]')
            
    ax2[count].imshow(tst_med**(1./3.), origin='lower', extent=[newpx[0],newpx[-1],newpy[0],newpy[-1]], cmap=cmc.batlowW_r, vmin=0, vmax=mxim)
    #hide y-axis
    if(count != 0):
        for k in np.arange(ny):
            ax2[count].get_yaxis().set_visible(False)
    ax2[count].set_xlim(rng, -rng)
    ax2[count].set_ylim(-rng, rng)
    ax2[count].set_xlabel(r'$\alpha$RA [mas]')
    ax2[count].set_ylabel(r'$\delta$DEC [mas]')
            
    ax3[count].imshow((difmed)**(1./3.), origin='lower', extent=[newpx[0],newpx[-1],newpy[0],newpy[-1]], cmap=cmc.vik, vmin=-mxlim, vmax=mxlim)
    #hide y-axis
    if(count != 0):
        for k in np.arange(ny):
            ax3[count].get_yaxis().set_visible(False)
    ax3[count].set_xlim(rng, -rng)
    ax3[count].set_ylim(-rng, rng)
    ax3[count].set_xlabel(r'$\alpha$RA [mas]')
    ax3[count].set_ylabel(r'$\delta$DEC [mas]')
            
    #ax3[count+3*nx].imshow(rmsmap, cmap='hot', origin='lower', extent=[newpx[0],newpx[-1],newpy[0],newpy[-1]])
    
    ax4[count].hist(score, bins=30, edgecolor='black', range=[0,1])
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('final score',score)
    print('average score',scoreavg, 'score std', scorestd)
    print('median score',np.median(score))
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    ax4[count].set_title('Score avg: '+r'{:.2f}'.format(scoreavg)+ ', std: ' + r'{:.2f}'.format(scorestd))
    #hide y-axis
    for k in np.arange(ny):
        ax4[count].get_yaxis().set_visible(False)
    ax4[count].set_xlabel(r'Score')
    
    

    # print(itx)
    # if itx >= 3:
    #plt.show()
        
        
    count+=1

# #hide x-axis
# for k in np.arange(ny):
#     if k != ny-1:
#         for l in np.arange(nx):
#             ax3[k*nx+l].get_xaxis().set_visible(False)
    
    
plt.figure(1)
plt.tight_layout()
plt.savefig('/home/fmillour/2024_IBC/PAPER/'+'figure5_model_'+contribs[plot-1]+'.png', dpi=my_dpi)
plt.figure(2)
plt.tight_layout()
plt.savefig('/home/fmillour/2024_IBC/PAPER/'+'figure5_entry_'+contribs[plot-1]+'.png', dpi=my_dpi)
plt.figure(3)
plt.tight_layout()
plt.savefig('/home/fmillour/2024_IBC/PAPER/'+'figure5_difference_'+contribs[plot-1]+'.png', dpi=my_dpi)
plt.figure(4)
plt.tight_layout()
plt.savefig('/home/fmillour/2024_IBC/PAPER/'+'figure5_score_'+contribs[plot-1]+'.png', dpi=my_dpi)
plt.show()