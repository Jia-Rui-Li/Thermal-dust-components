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

print("Statistics.py, Starting. ")
root_dir = os.path.abspath('.')

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# After removing CMB anisotropies and other foreground components, we have obtained thermal dust component from Planck HFI maps. 
# Before performing statistics, mask the regions contaminated by compact sources. 

def compact_sourece_mask(freq_str, Nside):
    # Compact source mask for single frequency channel
    # FWHM of compact sources from COM_PCCS_030/044/070_R2.04.fits (PCCS2), COM_PCCS_100/143/217/353/545/857_R2.01.fits (PCCS2), COM_PCCS_100/143/217/353/545/857-excluded_R2.01.fits (PCCS2E)
    # PCCS2 covers most of the sky and allows the user to produce subsamples at higher reliabilities than the target 80% integral reliability of the catalogue. 
    # PCCS2E contains sources detected in sky regions where the diffuse emission makes it difficult to quantify the reliability of the detections. 
    # A&A 594, A26 (2016)
    # https://wiki.cosmos.esa.int/planckpla2015/index.php/Catalogues

    if freq_str=="030" or freq_str=="044" or freq_str=="070":
        mask = numpy.zeros(12*Nside**2, dtype=int) + 1
        compact_catalogue = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"_R2.04.fits")
        galactic_longitude=compact_catalogue[1].data.field("GLON") # in unit of degree
        galactic_latitude =compact_catalogue[1].data.field("GLAT") # in unit of degree
        flux         = compact_catalogue[1].data.field("DETFLUX") # "DETFLUX": Flux density of source as determined by detection method, in unit of MJy
        FWHM         = compact_catalogue[1].data.field("GAU_FWHM_EFF") # "GAU_FWHM_EFF": Gaussian fit effective FWHM, in unit arcmin

    elif freq_str=="100" or freq_str=="143" or freq_str=="217" or freq_str=="353" or freq_str=="545" or freq_str=="857":
        mask = numpy.zeros(12*Nside**2, dtype=int) + 1
        compact_catalogue_1 = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"_R2.01.fits")
        # PCCS2, the sources have been detected in regions of the sky where it is possible to estimate the reliability of the detections, 
        # either statistically or by using external catalogues
        galactic_longitude_1=compact_catalogue_1[1].data.field("GLON") # in unit of degree
        galactic_latitude_1 =compact_catalogue_1[1].data.field("GLAT") # in unit of degree
        flux_1         = compact_catalogue_1[1].data.field("DETFLUX") # "DETFLUX": Flux density of source as determined by detection method, in unit of MJy
        FWHM_1         = compact_catalogue_1[1].data.field("GAU_FWHM_EFF") # "GAU_FWHM_EFF": Gaussian fit effective FWHM, in unit arcmin
        compact_catalogue_2 = astropy.io.fits.open(root_dir+"/PR2-2015/COM_PCCS_"+freq_str+"-excluded_R2.01.fits")
        # PCCS2E, the detected sources are located in regions of the sky where it is not possible to make an estimate of their reliability. 
        galactic_longitude_2=compact_catalogue_2[1].data.field("GLON") # in unit of degree
        galactic_latitude_2 =compact_catalogue_2[1].data.field("GLAT") # in unit of degree
        flux_2         = compact_catalogue_2[1].data.field("DETFLUX") # Flux density of source as determined by detection method
        FWHM_2         = compact_catalogue_2[1].data.field("GAU_FWHM_EFF") # GAU_FWHM_EFF, Gaussian fit effective FWHM, in unit arcmin
        # For 100/143/217/353/545/857 GHz, combine PCCS2 and PCCS2E
        galactic_longitude = numpy.concatenate((galactic_longitude_1,galactic_longitude_2), axis=0)
        galactic_latitude  =numpy.concatenate((galactic_latitude_1,galactic_latitude_2 ), axis=0)
        flux           =numpy.concatenate((flux_1, flux_2), axis=0)
        FWHM           =numpy.concatenate((FWHM_1, FWHM_2), axis=0)

    compact_source_data = numpy.array([galactic_longitude, galactic_latitude, flux, FWHM])
    #   Transpose: arrange by a certain row
    # No Transpose: arrange by a certain column
    # numpy.argsort(-compact_source_data[2,:]) returns the indices in reversely sorted order, according to the value of DETFLUX (2nd row of compact_source_data)
    compact_source_data = compact_source_data.T[numpy.argsort(-compact_source_data[2,:])].T
    # "number" is the number of comapact sources in each frequency channel
    number = compact_source_data.shape[1]
    # In each channel, 
    # mask radius = 2 FWHM for sources with the top 10% intensity
    for i in range(0, int(0.1*number), 1):
        # Direction of the i-th compact source
        # lonlat = True, input angles are assumed to be longitude and latitude in degree
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        # Mask radius from arcmin to radians
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=2*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0
    # mask radius = 1.5FWHM for sources with the top 10%-30% intensity
    for i in range(int(0.1*number), int(0.3*number), 1):
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=1.5*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0
    # mask radius = 1  FWHM for others
    for i in range(int(0.3*number), number, 1):
        vector = healpy.ang2vec(compact_source_data[0][i], compact_source_data[1][i], lonlat=True)
        ipix_disc = healpy.query_disc(nside=Nside, vec=vector, radius=1*compact_source_data[3][i]/60*numpy.pi/180)
        mask[ipix_disc] = 0
    return mask

# total_compact_source_mask combines compact source masks for each channel, with nside=2048
total_compact_source_mask = compact_sourece_mask("030", 2048)
for freq_str in ["044", "070", "100", "143", "217", "353", "545", "857"]:
    total_compact_source_mask = total_compact_source_mask*compact_sourece_mask(freq_str, 2048)
sky_fraction = sum(total_compact_source_mask)/(12*2048**2)

healpy.write_map(root_dir+"/results/compact_source_mask.fits", total_compact_source_mask, nest = False, coord = "G", overwrite="True", dtype=int)


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Function Mosaic() is for statistics
def Mosaic(color_correction, name, freq_str_1, freq_str_2, unit, nside1, nside2, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask):
# nside1 for output map
# nside2 for input map
# smooth_degree is the FWHM for smoothing, in degree
# disk_degree is the radius of the disk for statistics, in degree
# 0-th component of the output is "ratio": R(model) or R'(data)
# 1-st component of the output is "correlation": C_xy(model) or C_{x'y'}(data)
# 2-nd component of the output is used to mark points with a correlation of not less than 95% (only for dust data)
# name = "data"/"model"
    smooth_degree=float(smooth_degree)
    disk_degree = float(disk_degree)
    # dust_map_1 and dust_map_2 are (data/model) dust maps in neighbourinng frequency channels
    if color_correction == "yes":
        dust_map_1 = healpy.read_map(root_dir+"/results/dust_"+name+"_"+unit+"_"+freq_str_1+"GHz.fits", nest=False)
        dust_map_2 = healpy.read_map(root_dir+"/results/dust_"+name+"_"+unit+"_"+freq_str_2+"GHz.fits", nest=False)
    elif color_correction == "no":
        dust_map_1 = healpy.read_map(root_dir+"/results/No_color_correction_dust_"+name+"_"+unit+"_"+freq_str_1+"GHz.fits", nest=False)
        dust_map_2 = healpy.read_map(root_dir+"/results/No_color_correction_dust_"+name+"_"+unit+"_"+freq_str_2+"GHz.fits", nest=False)

    # Smooth dust_map_1 with beam = smooth_degree
    dust_map_1 = healpy.smoothing(dust_map_1, fwhm = smooth_degree*numpy.pi/180, pol=False)
    # Smooth dust_map_2 with beam = smooth_degree
    dust_map_2 = healpy.smoothing(dust_map_2, fwhm = smooth_degree*numpy.pi/180, pol=False)
    
    # Mask out galactic plane, with 60% sky coverage
    if galac_mask == "60":
        galaxy_mask = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=2)
    # Mask out galactic plane, with 70% sky coverage
    elif galac_mask == "70":
        galaxy_mask = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=3)
    # Mask out galactic plane, with 80% sky coverage
    elif galac_mask == "80":
        galaxy_mask = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=4)
    # field = 0 -- 20% sky coverage
    # field = 1 -- 40% sky coverage
    # field = 2 -- 60% sky coverage
    # field = 3 -- 70% sky coverage
    # field = 4 -- 80% sky coverage
    # field = 5 -- 90% sky coverage
    mask = galaxy_mask*total_compact_source_mask
    # If exclude the regions with latitude < 30 degree
    if low_mask == "yes": 
        low_galac_mask = numpy.zeros(12*2048**2, dtype=int)
        # disk_list is the region with abs(latitude) > 30 degree
        low_galac_list = healpy.query_strip(2048, 2/3*numpy.pi, 1/3*numpy.pi, inclusive=True)
        # valid region without low latitude region
        low_galac_mask[low_galac_list] = 1
        mask = mask*low_galac_mask
    if zodiacal_mask == "yes": 
        map_Zodiacal = 1e6*healpy.read_map(root_dir+"/PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal.fits")
        map_NoZodiacal = 1e6*healpy.read_map(root_dir+"/PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits")
        ecliptic_mask = numpy.zeros(12*2048**2, dtype=int)
        ecliptic_mask = ecliptic_mask + 1
        ecliptic_pix = numpy.where((map_Zodiacal - map_NoZodiacal) > 5)[0]
        ecliptic_mask[ecliptic_pix] = 0
        mask = mask*ecliptic_mask
    # points in ipix_disc form a sample set
    # ipix_disc with anchor point as the center and disk_degree as the radius
    # Initial values are set to  - 1.63750e+30, bad value for healpix, and will be masked by healpy.mollview.  
    anchor_point = numpy.zeros((3, 12*nside1**2), dtype=numpy.float64)
    anchor_point = anchor_point - 1.63750e+30
    func = lambda theta : 2*numpy.pi*numpy.sin(theta)
    # "area" is the area of ipix_disc  
    area = scipy.integrate.quad(func, 0, numpy.radians(disk_degree))[0]
    # "threshold" is the number of data points contained in 30% of the ipix_disc
    threshold = 0.3*12*nside2**2*area/(4*numpy.pi)
    # for ipix in anchor_point
    for ipix in range(0, 12*nside1**2, 1):
        angular = healpy.pix2ang(nside1, ipix, nest=False, lonlat=False)
        vector = healpy.ang2vec(angular[0], angular[1])
        ipix_disc = healpy.query_disc(nside=nside2, vec=vector, radius=numpy.radians(disk_degree))
        ipix_num = numpy.where(mask[ipix_disc]==1)[0]
        # ipix_area is the un-masked ipix_disc, area of the valid region
        ipix_area = ipix_disc[ipix_num]
        # if the area of valid region larger than 30% of ipix_disc, then the anchor point is valid
        if ipix_area.shape[0] > threshold:
            cov_matrix = numpy.cov(dust_map_1[ipix_area], dust_map_2[ipix_area], ddof=1)
            # the divisor used in the calculation is N - ddof, where N represents the number of elements.
            # numpy.cov(A, B)[0][0]: variance of data set A
            # numpy.cov(A, B)[1][1]: variance of data set B
            # numpy.cov(A, B)[0][1], numpy.cov(A, B)[1][0]: covariance of data set A and B
            # ratio is R(model) or R'(data)
            ratio = cov_matrix[0][1] / cov_matrix[0][0]
            # "correlation" is the correlation between dust_map_1 and dust_map_2
            correlation = cov_matrix[0][1] / numpy.sqrt(cov_matrix[0][0]*cov_matrix[1][1])
            anchor_point[0][ipix] = ratio
            anchor_point[1][ipix] = correlation
            if name == "data":
                if correlation >= 0.95:
                    # mark points with a correlation of not less than 95% (only for dust data)
                    anchor_point[2][ipix] = 1
    numpy.save(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+name+"_"+str(freq_str_1)+"GHz_"+str(freq_str_2)+"GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy",
              anchor_point)
    return 1


print("Statistics: ")
start_time = time.time()

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Mosaic(color_correction, name, freq_str_1, freq_str_2, unit, nside1, nside2, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
# with color correction
# take the galactic plane mask with 80% sky coverage
# smooth_degree=[1, 2, 3] degree, disk_degree=[5, 6, 7] degree
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 1, 5, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 1, 5, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 1, 5, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 1, 5, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 1, 5, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 1, 5, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 1, 6, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 1, 6, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 1, 6, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 1, 6, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 1, 6, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 1, 6, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 1, 7, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 1, 7, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 1, 7, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 1, 7, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 1, 7, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 1, 7, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)

arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 5, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 5, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 5, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 5, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 5, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 5, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 7, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 7, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 7, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 7, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 7, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 7, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)

arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 3, 5, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 3, 5, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 3, 5, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 3, 5, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 3, 5, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 3, 5, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 3, 6, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 3, 6, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 3, 6, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 3, 6, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 3, 6, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 3, 6, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 3, 7, "80", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 3, 7, "80", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 3, 7, "80", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 3, 7, "80", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 3, 7, "80", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 3, 7, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# with color correction
# take the galactic plane mask with 60%, 70% sky coverage
# smooth_degree=2 degree, disk_degree=6 degree
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "60", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "60", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "60", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "60", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "60", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "60", "no", "no"), 
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "70", "no", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "70", "no", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "70", "no", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "70", "no", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "70", "no", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "70", "no", "no")]
multiprocessing.Pool(processes=12).starmap(Mosaic, arguments)

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# with color correction
# take the galactic plane mask with 80% sky coverage
# smooth_degree=2 degree, disk_degree=6 degree
# dust model from Planck 2013 and A&A, 623, A21 (2019)
arguments = [
("yes", "model_2013", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2013", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2013", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2019", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2019", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("yes", "model_2019", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)
    
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# without color correction
# take the galactic plane mask with 80% sky coverage
# smooth_degree=2 degree, disk_degree=6 degree
arguments = [
("no", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2013", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2013", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2013", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2019", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2019", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "no"), 
("no", "model_2019", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "no")] 
multiprocessing.Pool(processes=12).starmap(Mosaic, arguments)

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# with color correction
# take the galactic plane mask with 80% sky coverage
# smooth_degree=2 degree, disk_degree=6 degree
# dust model from Planck 2015
# without low Galactic latitude
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "yes", "no"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "yes", "no"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "yes", "no"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "yes", "no"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "yes", "no"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "yes", "no")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# with color correction
# take the galactic plane mask with 80% sky coverage
# smooth_degree=2 degree, disk_degree=6 degree
# dust model from Planck 2015
# with zodiacal mask
arguments = [
("yes", "model_2015", "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "yes"), 
("yes", "data",       "100", "143", "muK_CMB", 64, 2048, 2, 6, "80", "no", "yes"), 
("yes", "model_2015", "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "yes"), 
("yes", "data",       "217", "353", "muK_CMB", 64, 2048, 2, 6, "80", "no", "yes"), 
("yes", "model_2015", "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "yes"), 
("yes", "data",       "545", "857", "MJy_sr",  64, 2048, 2, 6, "80", "no", "yes")]
multiprocessing.Pool(processes=6).starmap(Mosaic, arguments)

end_time = time.time()
print("Statistics, Succeed ! \n TIME: ", format( (end_time-start_time)/60, '.2f'),"min")
