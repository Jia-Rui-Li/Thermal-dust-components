import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import shutil
import sys
import multiprocessing
import matplotlib.pyplot

print("Inpain_Rebeam_HFI_maps.py, Starting")
root_dir = os.path.abspath('.')

# File name for each channel
def sky_map_file_name(freq_str):
    file_names = {
        "100": "HFI_SkyMap_100_2048_R3.01_full.fits",
        # in K_cmb
        "143": "HFI_SkyMap_143_2048_R3.01_full.fits",
        # in K_cmb
        "217": "HFI_SkyMap_217_2048_R3.01_full.fits",
        # in K_cmb
        "353": "HFI_SkyMap_353_2048_R3.01_full.fits",
        # in K_cmb
        "545": "HFI_SkyMap_545_2048_R3.01_full.fits",
        # in MJy/sr
        "857": "HFI_SkyMap_857_2048_R3.01_full.fits"
        # in MJy/sr        
    }
    # If freq_str exists in file_names, it returns the corresponding value, HFI_SkyMap_XXX_2048_R3.01_full.fits. 
    # If freq_str doesn't exist in file_names, it returns the default value 0.
    return file_names.get(freq_str, 0)

# FWHM from Planck 2018, I. A&A 641, A1(2020), Table 4
# beam_fwhm() gives FWHM for each channel, in radians
def beam_fwhm(freq_str):
    beam_dict = {
        "100": 9.66,
        "143": 7.22,
        "217": 4.90,
        "353": 4.92,
        "545": 4.67,
        "857": 4.22
    }
    beam = beam_dict.get(freq_str, 0)
    return beam / 60 * numpy.pi / 180
    
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# inpainting
# Nside is the nside of HFI map and mask map
Nside = 2048
# If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
neighbour_array = healpy.get_all_neighbours(nside=Nside, theta=numpy.arange(0, 12*Nside**2, 1), phi=None).T

def inpainted_sky_map(freq_str, N_iter):
# Inpaint the points where mask = 0
# N_iter is the iteration count
    N_iter = int(N_iter)
    # sky_map is HFI map
    sky_map = healpy.read_map(root_dir+"/PR3-2018/"+sky_map_file_name(freq_str), h=False, hdu=1, field=0)
    # mask_map
    mask_map_dict = {
        "100": 0,
        "143": 1,
        "217": 2,
        "353": 3,
        "545": 4,
        "857": 5
    }
    mask_field = mask_map_dict.get(freq_str)
    mask_map = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=mask_field)
    # values in mask.fits are Boolean
    # Boolean values to floats
    mask_map = mask_map.astype(numpy.float64)

    SkyMap = sky_map*mask_map
    # Filter out the points with mask = 0.
    ipix_mask = numpy.where(mask_map < 0.9)[0]
    # Get the neighboring pixels of the pixels that need to be inpainted
    neighbours = neighbour_array[ipix_mask]
    # valid_neighbours are the valid neighboring pixels
    # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
    valid_neighbours = neighbours >= 0
    # Count the neighboring pixels
    counts = numpy.sum(valid_neighbours, axis=1)

    for ii in range(0, N_iter):
        # Sum the values of the valid neighboring pixels
        sums = numpy.sum(SkyMap[neighbours] * valid_neighbours, axis=1)
        # Calculate the mean value of the neighboring pixels for each masked pixel
        mean_values = sums / counts
        # Update the SkyMap values for the pixels that need to be inpainted
        SkyMap[ipix_mask] = mean_values
    if freq_str in ["100", "143", "217", "353"]:
    # For 100, 143, 217, 353 GHz, output in muK_CMB
        SkyMap = 1e6*SkyMap
    return SkyMap


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Smooth and rebeam the inpainted HFI maps to FWHM = 9.66'
def rebeamed_inpainted_sky_map(freq_str, N_iter, Nside):
    Nside = int(Nside)
    SkyMap = inpainted_sky_map(freq_str, N_iter)
    # There are bad values in Planck HFI maps, they are set to - 1.63750e+30. This step is taken to remove bad values.
    bad_pixel = numpy.where(SkyMap<-1e+10)[0]
    for pixel in bad_pixel:
        neighbour = healpy.get_all_neighbours(nside=Nside, theta=pixel, phi=None)
        SkyMap[pixel] = numpy.mean(SkyMap[neighbour])

    Lmax = 3*Nside - 1
    beam1 = healpy.gauss_beam(fwhm=beam_fwhm(freq_str), lmax=Lmax)
    beam2 = healpy.gauss_beam(fwhm=9.66/60*numpy.pi/180, lmax=Lmax)
    alm_sky = healpy.map2alm(SkyMap, lmax=Lmax)
    alm_sky_smooth = healpy.smoothalm(alm_sky, beam_window=beam2/beam1)
    SkyMap = healpy.alm2map(alm_sky_smooth, nside=Nside)

    if freq_str in ["100", "143", "217", "353"]:
    # For 100, 143, 217, 353 GHz, in muK_CMB
        healpy.write_map(
            root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits", 
            SkyMap, nest = False, coord = "G", column_units = "muK_CMB", overwrite="True", dtype=numpy.float64)
    elif freq_str in ["545", "857"]:
    # For 545, 857 GHz, in MJy/sr
        healpy.write_map(
            root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits", 
            SkyMap, nest = False, coord = "G", column_units = "MJy_sr", overwrite="True", dtype=numpy.float64)
    return SkyMap

if os.path.exists(root_dir+"/results/"):
   shutil.rmtree(root_dir+"/results/")
os.mkdir(root_dir+"/results/")

print("Inpainting: ")
Nside = 2048
N_iter = int(sys.argv[1])
arguments = [("100", N_iter, Nside), ("143", N_iter, Nside), ("217", N_iter, Nside), ("353", N_iter, Nside), ("545", N_iter, Nside), ("857", N_iter, Nside)]
start_time = time.time()
multiprocessing.Pool(processes=6).starmap(rebeamed_inpainted_sky_map, arguments)
end_time = time.time()
print("Inpaint_Rebeam_HFI_maps.py, Success! \n TIME: ", format( (end_time-start_time)/60, '.2f'),"minutes")

