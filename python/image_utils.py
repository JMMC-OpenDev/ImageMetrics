import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import astropy.constants as cst
from scipy.ndimage import gaussian_filter
from scipy.interpolate import interp1d
from scipy.ndimage import rotate
from scipy import ndimage
from scipy import interpolate as intp
import cv2

rad2mas = 180./np.pi * 3600 * 1000


###############################################################

def resample_img_xy(img, px, py, newpx, newpy):
    print('Resampling image...')
    print(np.shape(img))
    # print('##################################################')
    # print(px)
    # print(py)
    # print('##################################################')
    # print(newpx)
    # print(newpy)
    # print('##################################################')
    
    nwl = np.shape(img)[0]
    wlen = np.arange(nwl)
    # First resample the image in x and y
    newimg = np.zeros((len(wlen), len(newpx), len(newpy)))
    X,Y = np.meshgrid(px, py, indexing='ij')
    nX,nY = np.meshgrid(newpx, newpy, indexing='ij')
    for iwl,w in enumerate(wlen):
        fi = intp.RegularGridInterpolator((px, py), img[iwl,:,:], bounds_error=False, fill_value=0)
        newimg[iwl,:,:] = fi((nX, nY))
    
        
    print(np.shape(newimg))
    print('Done!')
    return newimg

###############################################################

def compar_imageCubes(ref_img, ref_wlen, test_img, test_wlen):
    # Both images MUST have the same size in pixels (see interpolation above)
    print('Comparing images...')
    
    px = py = np.arange(np.shape(ref_img)[1])
    newpx = newpy = px
    # Cross-correlate test image with reference image looping on wavelength
    for wdx,w in enumerate(test_wlen):
        # First interpolate to test image wavelength
        nW, nX,nY = np.meshgrid(w, newpx, newpy, indexing='ij')
        fi = intp.RegularGridInterpolator((px, py), img[iwl,:,:], bounds_error=False, fill_value=0)
        newimg = fi((nW, nX, nY))
        
        plt.figure()
        plt.imshow(newimg)
        plt.show()
        
        
        rwlen_img = imgCube[wdx,:,:]
        twlen_img = imgCube[wdx,:,:]
        corr = cv2.filter2D(twlen_img, ddepth=-1, kernel=rwlen_img)
        max_rows, max_cols = np.where(corr == np.max(corr))
        
        if len(max_rows) != 0:
            twlen_img2 = np.roll(twlen_img, [-max_rows[0]-1-int(nux/2), -max_cols[0]-1-int(nux/2)], axis=[0,1])
            tim_cen[wdx,:,:] = twlen_img2
    
###############################################################

def compar_images(ref_img, test_img):
    """_summary_

    Args:
        ref_img (2D array of floats): Reference image
        test_img (2D array of floats): Test image

    Returns:
        _type_: _description_
    """
    # Both images MUST have the same size in pixels (see interpolation above)
    print('Comparing images...')
    nux = np.shape(test_img)[0]
    nuy = np.shape(test_img)[1]
    
    # Center images
    print('Centering images...')
    corr = cv2.filter2D(test_img, ddepth=-1, kernel=ref_img)
    max_rows, max_cols = np.where(corr == np.max(corr))
    if len(max_rows) != 0:
        test_img2 = np.roll(test_img, [-max_rows[0]-1-int(nux/2), -max_cols[0]-1-int(nuy/2)], axis=[0,1])
    
    l_1   = lambda im1, im2 : np.sum(np.abs(im1 - im2))
    m_s_e = lambda im1, im2 : np.sum(np.abs(im1 - im2)**2)
    
    # Leave the choice to put whatever comparison norm here
    #fnc = m_s_e
    fnc = l_1
    
    
    for i in range(3):
        
        print('#######################################################')
        
        #######################################################
        alpha = 1.0
        delta = 0.1
        grow = 0.5
        oldl1 = 1e99
        count = 0
        newl1 = fnc(ref_img, alpha * test_img2)
        
        print('looping on flux scaling factor...')
        niter = 0;
        while(np.abs(delta) > 1e-6 and niter < 100):
            niter+=1
            if niter > 100:
                print('Failed to converge')
            #print('alpha',alpha, 'delta', delta)
            newl1 = fnc(ref_img, (alpha - delta) * test_img2)
            newl2 = fnc(ref_img, alpha * test_img2)
            newl3 = fnc(ref_img, (alpha + delta) * test_img2)
            #print('oldl1',oldl1)
            #print('newl1',newl1, 'newl2',newl2, 'newl3',newl3)
            
            # test cases
            if(newl1 > newl2 > newl3):
                #print('entering case 1')
                alpha += delta
                delta *= 1.+grow
                oldl1 = newl3
            elif (newl3 > newl2 > newl1):
                #print('entering case 2')
                alpha -= delta
                delta *= -1.-grow
                oldl1 = newl1
            elif newl2 > newl1 > newl3:
                #print('entering case 3a')
                alpha += delta   
                delta *= 1.+grow       
                oldl1 = newl3
            elif newl2 > newl3 > newl1:
                #print('entering case 3b')
                alpha -= delta
                delta *= -1.-grow
                oldl1 = newl1
            elif newl2 <= newl1 and newl2 <= newl3:
                #print('entering case 4')
                delta *= (1. - grow)
            else:
                print('newl1', newl1)
                print('newl2', newl2)
                print('newl3', newl3)
                
                print('error!')
                prout()
        
        # Save image with new normalization
        test_img2 = alpha * test_img2
        print('flux found ',alpha, 'criterion value',fnc(ref_img, test_img2))
        
        #######################################################
        print('looping on x...')
        x=0
        delta = 1;
        switch = 0;
        niter = 0;
        while(delta > 1e-9 and niter < 100):
            niter+=1
            if niter >= 100:
                print('Failed to converge')
            newl1 = fnc(ref_img, np.roll(test_img2,[x - delta , 0], axis=[0,1]))
            newl2 = fnc(ref_img, np.roll(test_img2,[x , 0], axis=[0,1]))
            newl3 = fnc(ref_img, np.roll(test_img2,[x + delta , 0], axis=[0,1]))
            #print('newl1',newl1, 'newl2',newl2, 'newl3',newl3)
            
            # test cases
            if(newl1 > newl2 > newl3):
                if(delta < 0):
                    delta = 1
            elif (newl3 > newl2 > newl1):
                if(delta > 0):
                    delta = -1
            elif newl2 < newl1 and newl2 < newl3:
                break

            x += delta
        test_img2 = np.roll(test_img2,[x , 0], axis=[0,1])
        print('x',x, 'criterion value',fnc(ref_img, test_img2))
        
        #######################################################
        print('looping on y...')
        y=0
        delta = 1;
        switch = 0;
        niter = 0;
        while(delta > 1e-9 and niter < 100):
            niter+=1
            if niter >= 100:
                print('Failed to converge')
            newl1 = fnc(ref_img, np.roll(test_img2,[0,y - delta], axis=[0,1]))
            newl2 = fnc(ref_img, np.roll(test_img2,[0,y], axis=[0,1]))
            newl3 = fnc(ref_img, np.roll(test_img2,[0,y + delta ], axis=[0,1]))
            #print('newl1',newl1, 'newl2',newl2, 'newl3',newl3)
            
            # test cases
            if(newl1 > newl2 > newl3):
                if(delta < 0):
                    delta = 1
            elif (newl3 > newl2 > newl1):
                if(delta > 0):
                    delta = -1
            elif newl2 < newl1 and newl2 < newl3:
                break

            y += delta
        test_img2 = np.roll(test_img2,[0,y], axis=[0,1])
        print('y',y, 'criterion value',fnc(ref_img, test_img2))
        
    return test_img2, fnc(ref_img, test_img2)
    
    # fig, ax = plt.subplots(1,3)
    # ax[0].imshow(ref_img)
    # ax[1].imshow(test_img2)
    # ax[2].imshow(ref_img - test_img2)
    # plt.show()

###############################################################

def recenter_imgcube(imgCube):
    print('Centering image...')
    nux = np.shape(imgCube)[1]
    cen_im = np.zeros(np.shape(imgCube))
    refimg = imgCube[0,:,:]
    # recenter first wavelength to the image max
    max_rows, max_cols = np.where(refimg == np.max(refimg))
    refimg = np.roll(refimg, [-max_rows[0]-int(nux/2), -max_cols[0]-int(nux/2)], axis=[0,1])
    #plt.imshow(refimg)
    #plt.show()
    cen_im[0,:,:] = refimg
    # Cross-correlate and recenter all wavelengths to the first one
    for wdx,w in enumerate(imgCube[:,0,0]):
        if wdx > 0:
            wlen_img = imgCube[wdx,:,:]
            corr = cv2.filter2D(wlen_img, ddepth=-1, kernel=refimg)
            #plt.imshow(corr)
            #plt.show()
            max_rows, max_cols = np.where(corr == np.max(corr))
            if len(max_rows) != 0:
                wlen_img2 = np.roll(wlen_img, [-max_rows[0]-1-int(nux/2), -max_cols[0]-1-int(nux/2)], axis=[0,1])
                cen_im[wdx,:,:] = wlen_img2
            else:
                cen_im[wdx,:,:] = wlen_img
    print('Done!')
    return cen_im

###############################################################

def blur_img(img, beampix):
    # Blur image to the beam size
    blurred_im = np.zeros(np.shape(imgCube))
    blurred_im[wdx,:,:] = gaussian_filter(img, beampix, truncate=20)
    return blurred_im

###############################################################

def blur_imgCube(imgCube, beampix):
    # Blur image to the beam size
    print('Blurring image...')
    blurred_im = np.zeros(np.shape(imgCube))
    for wdx,w in enumerate(imgCube[:,0,0]):
        #print(wdx)
        blurred_im[wdx,:,:] = gaussian_filter(imgCube[wdx,:,:], beampix[wdx], truncate=20)
    print('Done!')
    return blurred_im

###############################################################

def calc_beamsize(wlen, Bmax, pixsize, fac=2):
    # Calculate beam size
    print('Calculating beam size...')
    beam = wlen * 1e-6 / Bmax / fac * rad2mas

    # Calculate beam size in pixels
    beampix = beam / pixsize
    beampix /= 2 * np.sqrt(2* np.log(2))
    #print('beam ',beam, 'beam in pixel', beampix)

    print('Done!')
    return beam, beampix

###############################################################

def pad_img(img, fac, px, py):
    print('Zero padding image...')
    npx = np.shape(img)[1]
    upx = int(np.log10(npx)/np.log10(2)+1) + fac
    nux = 2**upx
    npy = np.shape(img)[2]
    upy = int(np.log10(npy)/np.log10(2)+1) + fac
    nuy = 2**upy
    # Zero pad image to avoid issue when recentering it
    addpadx = int((nux-npx)/2)
    addpady = int((nuy-npy)/2)
    #print(addpad)
    pad_img = np.pad(img, ((0,0),(addpadx,addpadx + npx%2),(addpady,addpady + npy%2)),\
     mode='constant')
    print('pad img size',np.shape(pad_img))

    dpx = np.diff(px)[0]
    #print(dpx)
    #newpx = np.linspace(np.min(px) - addpadx*dpx, np.max(px) + (addpadx + npx%2)*dpx, nux)
    dpy = np.diff(py)[0]
    #newpy = np.linspace(np.min(py) - addpady*dpy, np.max(py) + (addpady + npy%2)*dpy, nux)
    newpx = np.arange(nux)*dpx-nux/2*dpx
    newpy = np.arange(nuy)*dpy-nuy/2*dpy
    
    print('Done!')
    return pad_img, newpx, newpy

###############################################################

def convert_to_mas(inp, unit):
    print('Converting angle unit to mas...')
    if unit=='deg':
        print('Done!') 
        return inp * 3600 * 1000;
    elif unit=='mas':
        print('Done!') 
        return inp;
    elif unit.lower()=='rad':
        print('Done!') 
        return inp * rad2mas;
    else:
        print('failed to find a proper angle unit, dying loudly...')
        error()
    
def convert_to_micron(inp, unit):
    print('Converting wavelength unit to microns...')
    if unit.lower()=='micron':
        print('Done!') 
        return inp;
    elif unit.lower()=='meter':
        print('Done!') 
        return inp * 1e6;
    else:
        print('failed to find a proper length unit, dying loudly...')
        error()

###############################################################

def load_img(filename):
    print('loading',filename)
    fh = fits.open(filename)
    img = fh[0].data
    hdr = fh[0].header
    npx = hdr['NAXIS1']
    cdx = hdr['CDELT1']
    crx = hdr['CRPIX1']
    cvx = hdr['CRVAL1']
    try:
        cux = hdr['CUNIT1']
    except:
        print("NO CUNIT1, setting it to radians")
        cux = 'rad';

    npy = hdr['NAXIS2']
    cdy = hdr['CDELT2']
    cry = hdr['CRPIX2']
    cvy = hdr['CRVAL2']
    try:
        cuy = hdr['CUNIT2']
    except:
        print("NO CUNIT2, setting it to radians")
        cuy = 'rad';

    try:
        nbl = hdr['NAXIS3']
    except:
        print("NO NAXIS3")
        nbl = 1;

    try:
        cdl = hdr['CDELT3']
    except:
        print("NO CDELT3, setting it to 1")
        cdl = 1.;
    try:
        crl = hdr['CRPIX3']
    except:
        print("NO CRPIX3, setting it to 1")
        crl = 1.;
    try:
        cvl = hdr['CRVAL3']
    except:
        print("NO CRVAL3, setting it to 1")
        cvl = 1.;
    try:
        cul = hdr['CUNIT3']
    except:
        print("NO CUNIT3, setting it to microns")
        cul = 'micron';

    '''
    print('nbpx',npx, npy, nbl)
    print('delta',cdx, cdy,cdl)
    print('ref px',crx,cry,crl)
    print('ref val', cvx,cvy,cvl)
    print('units',cux,cuy,cul)
    print('delta in base unit', convert_to_mas(cdx, cux), convert_to_mas(cdy, cuy), convert_to_micron(cdl, cul))
    '''
    px   = convert_to_mas((np.arange(npx) - crx) * cdx + cvx, cux)
    py   = convert_to_mas((np.arange(npy) - cry) * cdy + cvy, cuy)
    if nbl > 1:
        wlen = convert_to_micron((np.arange(nbl) - crl  + 1) * cdl + cvl, cul)
    else:
        wlen = cvl

    if cdl == 1 and crl == 1 and cvl == 1:
        print("Trying Drevon's inputs")
        try:
            wlen = fh[1].data['WAVELENGTH']
            print("...Success!")
            
        except:
            print("...Aaaaaargh!!...")

    #print('nX',npx,'nY',npy,'nlambda',nbl, 'dims',np.shape(img))
    #print("OK")
    if nbl > 1:
        #print("Getting mean image...")
        #img_avg = np.mean(img,axis=0)$
        if any(np.diff(wlen) < 0):
            print("Sorting wavelengths by increasing order")
            # Sort wavelength in increasing order
            idx = np.argsort(wlen)
            wlen = wlen[idx]
            img = img[idx,:,:]

    print('Done loading image!')
    return img, px, py, wlen
