import sys,os, pdb



# Basic system utilities for interacting with files
import glob
import time
import shutil
import warnings
import zipfile
import urllib.request
from astropy import units as u
from astropy.stats import sigma_clip
#from reproject import reproject_exact

# Astropy utilities for opening FITS and ASCII files
from astropy.io import fits
from astropy.io import ascii
from astropy.utils.data import download_file
# Astropy utilities for making plots

# Numpy for doing calculations
import numpy as np

# Matplotlib for making plots
import matplotlib.pyplot as plt
from matplotlib import rc

#import reproject_nirspec,crop_cube_NIRSpec
from skimage import measure
from matplotlib.pyplot import *

def make_mask(sci_cal,spec2_masked_dir,mask_high = 60,mask_negative = -10,run_spec2=False):
    j = []
    for i in sci_cal:
        j.append(fits.open(i)[1].data)

    k = np.median(j,axis=0)  # modified from dither to mean for 2 dithers in the survey program

    mask = np.ones(k.shape)
    mask[np.where(k<mask_negative)] = 0.
    loc_bad = np.where(mask == 0.)
    for i in sci_cal:
        hdu = fits.open(i)
        hdu[1].data[loc_bad] = np.nan
        hdu[2].data[loc_bad] = np.nan
        hdu[3].data[loc_bad] = 1e9

        loc_high = np.where(abs(hdu[1].data)>mask_high)

        hdu[1].data[loc_high] = np.nan
        hdu[2].data[loc_high] = np.nan
        hdu[3].data[loc_high] = 1e9

        hdu.writeto(spec2_masked_dir+i[-52:],overwrite=True)
        hdu.close()
