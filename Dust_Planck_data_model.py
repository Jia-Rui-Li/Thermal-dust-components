import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import sys
import multiprocessing
print("Dust_Planck_data_model.py, Starting")
root_dir = os.path.abspath('.')

# transmission spectra of Planck HFI
# wavenumber in cm^-1
# transmission is normalized to 1 at the maximum for HFI
# https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/The_RIMO
HFI_RIMO = astropy.io.fits.open(root_dir+"/PR3-2018/HFI_RIMO_R3.00.fits")

# converse wavenumber in cm^-1 to frequency in GHz by multiplying by 29979245800/1e9
# 100 GHz
freq_100 = 29979245800/1e9*HFI_RIMO["BANDPASS_F100"].data.field("wavenumber")
trans_100 = HFI_RIMO["BANDPASS_F100"].data.field("transmission")
# 143 GHz
freq_143 = 29979245800/1e9*HFI_RIMO["BANDPASS_F143"].data.field("wavenumber")
trans_143 = HFI_RIMO["BANDPASS_F143"].data.field("transmission")
# 217 GHz
freq_217 = 29979245800/1e9*HFI_RIMO["BANDPASS_F217"].data.field("wavenumber")
trans_217 = HFI_RIMO["BANDPASS_F217"].data.field("transmission")
# 353 GHz
freq_353 = 29979245800/1e9*HFI_RIMO["BANDPASS_F353"].data.field("wavenumber")
trans_353 = HFI_RIMO["BANDPASS_F353"].data.field("transmission")
# 545 GHz
freq_545 = 29979245800/1e9*HFI_RIMO["BANDPASS_F545"].data.field("wavenumber")
trans_545 = HFI_RIMO["BANDPASS_F545"].data.field("transmission")
# 857 GHz
freq_857 = 29979245800/1e9*HFI_RIMO["BANDPASS_F857"].data.field("wavenumber")
trans_857 = HFI_RIMO["BANDPASS_F857"].data.field("transmission")


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*

# When convolving signal with transmission spectra, integration should not start from the minimum value, as it would introduce significant computational errors.
# Configuration related to integration over frequency for HFI channels
## Lower limit of integration is 1/5 of the central frequency
## Upper limit of integration is 5 of the central frequency
# The selection of lower limit and upper limit of integeration ensures that the transmission larger than 1e-6. 
# There are some negative values in trans_143, 
# (num = [229 230 231 232 233 234 235 236 252 253 254 255 256 257 258 275 276 277 278 279])
# we have excluded them. 
# function freq_low_limit() gives the index of the lower integration limit for each channel
#def freq_low_limit(freq_str):
#    if freq_str == "100":
#        return 39
#    elif freq_str == "143":
#        return 56
#    elif freq_str == "217":
#        return 85
#    elif freq_str == "353":
#        return 139
#    elif freq_str == "545":
#        return 215
#    elif freq_str == "857":
#        return 339
#    else:
#        return 1

def freq_low_limit(freq_str):
    if freq_str == "100":
        return 131
    elif freq_str == "143":
        return 280
    elif freq_str == "217":
        return 321
    elif freq_str == "353":
        return 564
    elif freq_str == "545":
        return 262
    elif freq_str == "857":
        return 206
    else:
        return 1


# function freq_high_limit() gives the index of the upper integration limit for each channel
#def freq_high_limit(freq_str):
#    if freq_str == "100":
#        return 630
#    elif freq_str == "143":
#        return 1078
#    elif freq_str == "217":
#        return 1323
#    elif freq_str == "353":
#        return 2099
#    elif freq_str == "545":
#        return 3755
#    elif freq_str == "857":
#        return 7677
#    else:
#        return eval("freq_"+freq_str).shape[0]-1

def freq_high_limit(freq_str):
    if freq_str == "100":
        return 494
    elif freq_str == "143":
        return 789
    elif freq_str == "217":
        return 847
    elif freq_str == "353":
        return 1449
    elif freq_str == "545":
        return 2442
    elif freq_str == "857":
        return 5622
    else:
        return eval("freq_"+freq_str).shape[0]-1

    
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Functions for unit conversion: x_trans(), b_prime_RJ(), b_prime_CMB()

# Planck constant in SI unit
h_Planck = 6.62607015e-34
# Speed of light in SI unit
c_light = 299792458
# Boltzmann constant in SI unit
k_Boltzman = 1.380649e-23
# Average temperature of CMB in SI unit
T_0 = 2.7255
# Reference: A&A 571, A9 (2014) P12

# x_trans = h \nu/kT
def x_trans(freq):
# freq in GHz
    freq = float(freq)
    return h_Planck*1e9*freq/k_Boltzman/T_0

# From Rayleigh-Jeans unit to SI unit
def b_prime_RJ(freq):
# freq in GHz
    freq = float(freq)
    return 2*k_Boltzman*(1e9*freq)**2/c_light**2

# From CMB unit to SI unit
def b_prime_CMB(freq):
# freq in GHz
    freq = float(freq)
    x = x_trans(freq)
    return 2*k_Boltzman*(1e9*freq)**2/c_light**2*x**2*numpy.exp(x)/(numpy.exp(x)-1)**2

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# FWHM from Planck 2018, I. A&A 641, A1(2020), Table 4
# beam_fwhm() gives FWHM for each channel, in radians
def beam_fwhm(freq_str):
    if freq_str == "100":
        beam = 9.66
    elif freq_str == "143":
        beam = 7.22
    elif freq_str == "217":
        beam = 4.90
    elif freq_str == "353":
        beam = 4.92
    elif freq_str == "545":
        beam = 4.67
    elif freq_str == "857":
        beam = 4.22
    else:
        beam = 0
    return beam/60*numpy.pi/180

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of CMB anisotropies
def CMB_map(freq_str, unit, Nside):
    Nside = int(Nside)
    Lmax = 3*Nside - 1
    # Inpainted I intensity map of CMB (SMICA)
    CMB_map = healpy.read_map(root_dir+"/PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits", h=False, field=5) # in unit of K_cmb
    
    if unit == "muK_CMB":
        CMB_map = 1e6*CMB_map
    elif unit == "MJy_sr":
        part1 = 0
        part2 = 0
        # Transmission spectrum of HFI detectors
        freq_array = eval("freq_"+freq_str)
        trans_array= eval("trans_"+freq_str)
        for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
            freq = freq_array[i]
            delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
            trans = trans_array[i]
            # part1 is numerator
            part1 = part1 + b_prime_CMB(freq)*trans*delta_freq
            # part2 is denominator
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        # For CMB anisotropies, T_CMB is independent with frequency
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        CMB_map = 1e-6*1e26*CMB_map*part1/part2
    beam1 = healpy.gauss_beam(fwhm=5/60*numpy.pi/180, lmax=Lmax, pol=False)
    beam2 = healpy.gauss_beam(fwhm=9.66/60*numpy.pi/180, lmax=Lmax, pol=False)
    alm_CMB = healpy.map2alm(CMB_map, lmax=Lmax)
    alm_CMB_smooth = healpy.smoothalm(alm_CMB, beam_window=beam2/beam1)
    CMB_map = healpy.alm2map(alm_CMB_smooth, nside=Nside)
    return CMB_map

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of free-free

# FWHM = 60', Planck 2015, X. A&A 594, A10 (2016), Table 5
# Emission Measure of free-free in pc cm^-6
free_free_EM = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits", h=False, field=0, dtype=numpy.float64)
# Temperature of electrons, in unit of K. 
free_free_Te = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits", h=False, field=3, dtype=numpy.float64)
T4 = free_free_Te/1e4

def gaunt_ff(freq):
# Gaunt factor of free-free, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz, is \nu9 in Planck 2015, X. A&A 594, A10 (2016), Table 4
    freq = float(freq)
    gaunt = numpy.log(numpy.exp(5.960 - numpy.sqrt(3)/numpy.pi*numpy.log(freq/T4**1.5)) + numpy.exp(1))
    return gaunt

def tau_ff(freq):
# Optical depth of free-free, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz, is \nu9 in Planck 2015, X. A&A 594, A10 (2016), Table 4
    freq = float(freq)
    tau = 0.05468/free_free_Te**1.5/freq**2*free_free_EM*gaunt_ff(freq)
    return tau

def K_RJ_ff(freq):
# Brightness temperature of free-free in unit of K_RJ, Planck 2015, X. A&A 594, A10 (2016), Table 4
# freq in GHz
    freq = float(freq)
    s_ff = free_free_Te*(1-numpy.exp(-tau_ff(freq)))
    return s_ff

def intensity_ff(freq):
# Intensity of free-free at freq in W m^-2 Hz^-1 sr^-1
    freq = float(freq)
    return b_prime_RJ(freq)*K_RJ_ff(freq)

def reading_ff(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB or MJy_sr, and the output with Nside
    Nside = int(Nside)
    Lmax = 3*Nside - 1
    # Nside = 256 for COM_CompMap_freefree-commander_0256_R2.00.fits
    part1 = numpy.zeros(12*256**2, dtype=numpy.float64)
    part2 = 0
    # Transmission spectrum of HFI detectors
    freq_array = eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        # part1 is numerator
        part1 = part1 + intensity_ff(freq)*trans*delta_freq
        if unit == "muK_CMB":
            # part2 is denominator
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        map_ff = 1e6*part1/part2
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        map_ff = 1e-6*1e26*part1/part2
    
    # Up_grade from nside=256 to nside=Nside
    alm_ff = healpy.map2alm(map_ff)
    map_ff_new = healpy.alm2map(alm_ff, nside=Nside)
    return map_ff_new


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of synchrotron

# Extension = "SYNC-TEMP" is the template for synchrotron radiation (Intensity VS Frequency)
synchro_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits")

# Reference Frequency of synchrotron Stokes I is 408.0 MHz
synchro_freq_ref = 0.408
synchro_freq = synchro_fits["SYNC-TEMP"].data.field("NU") # GHz
synchro_intensity = synchro_fits["SYNC-TEMP"].data.field("I") # W Hz^-1 m^-2 sr^-1
# synchro_template is the spline interpolation obtained based on synchro_freq and synchro_intensity, representing the relationship between synchrotron intensity and frequency (in log-log).
synchro_template = scipy.interpolate.splrep(numpy.log10(synchro_freq), numpy.log10(synchro_intensity), k=2)
# synchro_intensity_ref is the synchrotron intensity at 408 MHz in SI unit
synchro_intensity_ref = 10**scipy.interpolate.splev(numpy.log10(synchro_freq_ref), synchro_template)


# map of synchrotron radiation in unit of muK_RJ
synchrotron_I = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits", h=False, field=0) # in muK_RJ, at 408 MHz
# From muK_RJ to SI unit
synchrotron_I = 1e-6*b_prime_RJ(synchro_freq_ref)*synchrotron_I

def reading_synch(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB or MJy_sr, and the output with Nside
    Nside = int(Nside)
    part1 = 0
    part2 = 0
    # Transmission spectrum of HFI detectors
    freq_array = eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        # synchro_intensity is synchrotron intensity at freq in SI unit
        synchro_intensity = 10**scipy.interpolate.splev(numpy.log10(freq), synchro_template)
        # Suppose that the spectrum of synchrotron is identical alone all directions
        # part1 is numerator
        part1 = part1 + synchro_intensity/synchro_intensity_ref*trans*delta_freq
        if unit == "muK_CMB":
            # part2 is denominator
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        map_synchro = 1e6*synchrotron_I*part1/part2
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        map_synchro = 1e-6*1e26*synchrotron_I*part1/part2
    
    # Up_grade from nside=256 to nside=Nside
    alm_synchro = healpy.map2alm(map_synchro)
    map_synchro_new = healpy.alm2map(alm_synchro, nside=Nside)
    return map_synchro_new


#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of CO line emission

# Unit conversion of CO line emission, from K_RJ km s^-1 to muK_CMB or MJy_sr
def unit_conversion_CO(freq_str, unit, CO_str):
    # CO_str = "CO_10" / "CO_21" / "CO_32"
    
    # Emission frequencies of CO_10, CO_21, CO_32 are 115.271 GHz, 230.538 GHz, 345.796 GHz, respectively. 
    freq_CO_10 = 115.271
    freq_CO_21 = 230.538
    freq_CO_32 = 345.796
    
    # Emission frequency
    freq_CO = eval("freq_"+CO_str)
    # Transmission spectrum of HFI detectors
    freq_array= eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    # i1 (left), i2 (right) are the indices of the frequencies closest to the CO emission frequency.
    i1 = numpy.where(freq_array < freq_CO)[0][-1]
    i2 = numpy.where(freq_array > freq_CO)[0][0]
    r1 = (freq_array[i2] - freq_CO)/(freq_array[i2] - freq_array[i1])
    r2 = (freq_CO - freq_array[i1])/(freq_array[i2] - freq_array[i1])
    # trans_CO is the transmission at CO emission frequency
    trans_CO = r1*trans_array[i1] + r2*trans_array[i2]

    # freq_CO and delta_freq in GHz
    # speed of light in km/s
    # part1 is numerator
    part1 = b_prime_RJ(freq_CO)*(freq_CO/299792.458)*trans_CO
    # part2 is denominator
    part2 = 0
    for i in range(1, freq_array.shape[0]-2, 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        if unit == "muK_CMB":
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        return 1e6*part1/part2
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        return 1e-6*1e26*part1/part2

def reading_CO(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB or MJy_sr, and the output with Nside, FWHM = 9.66'
    Nside = int(Nside)
    Lmax = 3*Nside - 1
    # High-resolution CO(2-1) line from COMMANDER, in K_RJ km/s, with FWHM = 7.5'
    CO21 = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits", h=False, field=0)
    beam1 = healpy.gauss_beam(fwhm=7.5/60*numpy.pi/180, lmax=Lmax, pol=False)
    beam2 = healpy.gauss_beam(fwhm=9.66/60*numpy.pi/180, lmax=Lmax, pol=False)
    alm_CO21 = healpy.map2alm(CO21, lmax=Lmax)
    alm_CO21_smooth = healpy.smoothalm(alm_CO21, beam_window=beam2/beam1)
    CO21 = healpy.alm2map(alm_CO21_smooth, nside=Nside)
    # Zero out low-signal regions (<1 K_RJ km/s) of CO(2-1) map, ApJ 798, 88 (2015) P3
#    pixel_CO21_low = numpy.where(CO21 < 1)[0]
#    CO21[pixel_CO21_low] = 0
    # CO_2_1 / CO_1_0 = 0.595
    # CO_3_2 / CO_1_0 = 0.297
    # Planck 2013 XIII, P7 and Planck 2018 IV, P4
    CO10_new = unit_conversion_CO(freq_str, unit, "CO_10")*1.000/0.595*CO21
    CO21_new = unit_conversion_CO(freq_str, unit, "CO_21")*0.595/0.595*CO21
    CO32_new = unit_conversion_CO(freq_str, unit, "CO_32")*0.297/0.595*CO21
    CO_map = CO10_new + CO21_new + CO32_new
    return CO_map

##--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Map of 94/100 GHz line
#xline_fits = astropy.io.fits.open(root_dir+"/PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits")
#xline_fits.info()
#xline_fits["comp-map-xline"].header

def reading_xline(freq_str, unit, Nside):
# Reading of HFI detector at each channel, in unit of muK_CMB, and the output with Nside
# Non-zero at 100 GHz
    Nside = int(Nside)
    if unit == "MJy_sr":
        map_xline_new = numpy.zeros(12*Nside**2, dtype=float)
    elif unit == "muK_CMB":
        if freq_str == "100":
            map_xline = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits", h=False, field=0) # muK_CMB, FWHM = 60.0'
            # field = 0: Intensity map
            # field = 1: Mean Intensity
            alm_xline = healpy.map2alm(map_xline)
            map_xline_new = healpy.alm2map(alm_xline, nside=Nside)
        else:
            map_xline_new = numpy.zeros(12*Nside**2, dtype=float)
    return map_xline_new

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# HFI maps removed out CMB anisotropies and other foreground components
def reading_data_dust(freq_str, unit, Nside):
    Nside = int(Nside)
    rebeamed_inpainted_sky_map = healpy.read_map(root_dir+"/results/inpainted_rebeamed_map_"+freq_str+"GHz_Nside_"+str(Nside)+".fits")
    dust_map = rebeamed_inpainted_sky_map - CMB_map(freq_str, unit, Nside) - reading_ff(freq_str, unit, Nside) - reading_synch(freq_str, unit, Nside) - reading_CO(freq_str, unit, Nside) - reading_xline(freq_str, unit, Nside)
    healpy.write_map(root_dir+"/results/dust_data_"+unit+"_"+freq_str+"GHz.fits", 
                     dust_map, nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return dust_map

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Reading of HFI detector at each channel, for thermal dust emission from single component model
def reading_model_dust(freq_str, year, unit, Nside):
# in unit of muK_CMB or MJy_sr, with FWHM = 9.66', and the output with Nside
# If year = "2013", single component model from Planck 2013 release
# If year = "2015", single component model from Planck 2015 release
# If year = "2019", single component model from A&A 623, A21 (2019), https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/ 
    Nside = int(Nside)
    Lmax = 3*Nside - 1
    if year == "2013": 
        # Optical depth of dust at 353 GHz from Planck 2013 release
        dust_tau = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=0)
        # Dust emission spectral index from Planck 2013
        dust_index = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=6)
        # Dust equilibrium temperature of dust from Planck 2013 release, in unit of K
        dust_temperature = healpy.read_map(root_dir+"/PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits", h=False, field=4)        
    elif year == "2015":
        # Optical depth at 353 GHz from Planck 2015, FWHM = 5'
        dust_tau = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits", h=False, field=0)
        # Dust emission spectral index from Planck 2015, FWHM = 5'
        dust_index = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits", h=False, field=0)
        # Dust equilibrium temperature from Planck 2015, in unit of K, FWHM = 5'
        dust_temperature = healpy.read_map(root_dir+"/PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits", h=False, field=0)
    elif year == "2019":
        # Optical depth at 353 GHz from A&A 623, A21 (2019), FWHM = 5', https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/ReadMe
        tau = healpy.read_map(root_dir+"/Melis/tau.fits", h=False, field=0)
        # Dust emission spectral index from A&A 623, A21 (2019), FWHM = 5'
        beta = healpy.read_map(root_dir+"/Melis/beta.fits", h=False, field=0)
        # Dust equilibrium temperature from A&A 623, A21 (2019), in unit of K, FWHM = 5'
        temp = healpy.read_map(root_dir+"/Melis/temp.fits", h=False, field=0)
        # There are NaN values in the model from Melis, we need to inpaint them. 
        dust_model_Melis = numpy.array([beta, tau, temp])
        Nside = 2048
        # If phi is not given or None, theta is interpreted as pixel number, otherwise, theta, phi are angular coordinates in radians.
        neighbour_array = healpy.get_all_neighbours(nside=Nside, theta=numpy.arange(0, 12*Nside**2, 1), phi=None).T
        N_iter = 2000
        # N_iter is the count of iteration for inpainting. 
        pixel_with_nan = numpy.where(numpy.isnan(temp))[0]
        # Get the neighboring pixels of the pixels that need to be inpainted
        neighbours = neighbour_array[pixel_with_nan]
        # valid_neighbours are the valid neighboring pixels
        # Since in the healpix algorithm, a pixel has 8 or 7 neighbours, when a pixel has 7 neighbours, healpy.get_all_neighbours() returns a -1 to mark the invalid neighbour
        valid_neighbours = neighbours >= 0
        # Count the neighboring pixels
        counts = numpy.sum(valid_neighbours, axis=1)
        dust_model_Melis[:, pixel_with_nan] = 0
        for ii in range(0, N_iter):
            # Sum the values of the valid neighboring pixels
            sums = numpy.sum(dust_model_Melis[:, neighbours] * valid_neighbours, axis=2)
            # Calculate the mean value of the neighboring pixels for each masked pixel
            mean_values = sums / counts
            # Update the values with NaN
            dust_model_Melis[:, pixel_with_nan] = mean_values
        beta, tau, temp = dust_model_Melis
        dust_tau = tau
        dust_index = beta
        dust_temperature = temp
    
    # Reference frequency for single component model
    freq_d = 353
    # A&A 596, A109 (2016) Equation (8)
    # https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/ReadMe

    # part1 is numerator
    part1 = numpy.zeros(12*2048**2, dtype=numpy.float64)
    # part2 is denominator
    part2 = 0
    # Transmission spectrum of HFI detectors
    freq_array = eval("freq_"+freq_str)
    trans_array= eval("trans_"+freq_str)
    
    for i in range(freq_low_limit(freq_str), freq_high_limit(freq_str), 1):
        freq = freq_array[i]
        delta_freq = 0.5*(freq_array[i+1] - freq_array[i-1])
        trans = trans_array[i]
        x_trans_freq = h_Planck*(1e9*freq)/k_Boltzman/dust_temperature
        Black_body_spectrum = 2*h_Planck*(1e9*freq)**3/c_light**2/(numpy.exp(x_trans_freq)-1)
        part1 = part1 + dust_tau*(freq/freq_d)**dust_index*Black_body_spectrum*trans*delta_freq
        if unit == "muK_CMB":
            part2 = part2 + b_prime_CMB(freq)*trans*delta_freq
        elif unit == "MJy_sr":
            part2 = part2 + (freq/float(freq_str))**(-1)*trans*delta_freq
        
    if unit == "muK_CMB":
        # From K_CMB to muK_CMB
        dust_map = 1e6*part1/part2
    elif unit == "MJy_sr":
        # From W m^-2 Hz^-1 sr^-1 to MJy/sr
        dust_map = 1e-6*1e26*part1/part2
    # Re-beam to FWHM = 9.66'
    beam1 = healpy.gauss_beam(fwhm=5/60*numpy.pi/180, lmax=Lmax, pol=False)
    beam2 = healpy.gauss_beam(fwhm=9.66/60*numpy.pi/180, lmax=Lmax, pol=False)
    alm_dust = healpy.map2alm(dust_map, lmax=Lmax)
    alm_dust_smooth = healpy.smoothalm(alm_dust, beam_window=beam2/beam1)
    dust_map = healpy.alm2map(alm_dust_smooth, nside=Nside)
    healpy.write_map(root_dir+"/results/dust_model_"+year+"_"+unit+"_"+freq_str+"GHz.fits", 
dust_map, nest = False, coord = "G", column_units = unit, overwrite="True", dtype=numpy.float64)
    return dust_map

#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*

Nside = 2048

print("Dust map from observed data: ")
start_time = time.time()
arguments = [("100", "muK_CMB", Nside), ("143", "muK_CMB", Nside), ("217", "muK_CMB", Nside), ("353", "muK_CMB", Nside), ("545", "MJy_sr", Nside), ("857", "MJy_sr", Nside)]
multiprocessing.Pool(processes=6).starmap(reading_data_dust, arguments)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')

year = "2013"
print("Dust map from model (Planck 2013): ")
start_time = time.time()
arguments = [("100", year, "muK_CMB", Nside), ("143", year, "muK_CMB", Nside), ("217", year, "muK_CMB", Nside), ("353", year, "muK_CMB", Nside), ("353", year, "MJy_sr", Nside), ("545", year, "MJy_sr", Nside), ("857", year, "MJy_sr", Nside)]
multiprocessing.Pool(processes=7).starmap(reading_model_dust, arguments)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')

year = "2015"
print("Dust map from model (Planck 2015): ")
start_time = time.time()
arguments = [("100", year, "muK_CMB", Nside), ("143", year, "muK_CMB", Nside), ("217", year, "muK_CMB", Nside), ("353", year, "muK_CMB", Nside), ("353", year, "MJy_sr", Nside), ("545", year, "MJy_sr", Nside), ("857", year, "MJy_sr", Nside)]
multiprocessing.Pool(processes=7).starmap(reading_model_dust, arguments)
end_time = time.time()
print('TIME: ', format((end_time-start_time)/60, '.2f'), 'min')

year = "2019"
print("Dust map from model (A&A 623, A21, 2019): ")
start_time = time.time()
arguments = [("100", year, "muK_CMB", Nside), ("143", year, "muK_CMB", Nside), ("217", year, "muK_CMB", Nside), ("353", year, "muK_CMB", Nside), ("353", year, "MJy_sr", Nside), ("545", year, "MJy_sr", Nside), ("857", year, "MJy_sr", Nside)]
multiprocessing.Pool(processes=7).starmap(reading_model_dust, arguments)
end_time = time.time()
print("Dust_Planck_data_model.py, Succeed! \n TIME: ", format((end_time-start_time)/60, '.2f'), "min")
