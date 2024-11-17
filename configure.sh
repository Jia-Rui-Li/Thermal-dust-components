#!/bin/bash

sudo apt-get update
sudo apt-get install python3
sudo apt-get install python3-pip
pip3 install numpy
pip3 install matplotlib
pip3 install healpy
pip3 install scipy
pip3 install astropy

# ------------------------------------
# ------------------------------------
# Planck 2013 release
if [[ -d 'PR1-2013' ]]
then
    echo 'You have directory PR1-2013'
else
    echo 'You have not directory PR1-2013, we will create it. '
    mkdir PR1-2013
fi

# ------------------------------------
# Thermal dust model (Planck 2013)
if [[ -f PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits ]]
then
    echo 'You have file HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
else
    echo 'You have not file HFI_CompMap_ThermalDustModel_2048_R1.20.fits, we will download it. '
    wget -cO PR1-2013/HFI_CompMap_ThermalDustModel_2048_R1.20.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/maps/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
fi

if [[ -f PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal.fits ]]
then
    echo 'You have file HFI_SkyMap_217_2048_R1.10_nominal.fits'
else
    echo 'You have not file HFI_SkyMap_217_2048_R1.10_nominal.fits, we will download it. '
    wget -cO PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/maps/HFI_SkyMap_217_2048_R1.10_nominal.fits'
fi

if [[ -f PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits ]]
then
    echo 'You have file HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits'
else
    echo 'You have not file HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits, we will download it. '
    wget -cO PR1-2013/HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_1/all-sky-maps/maps/HFI_SkyMap_217_2048_R1.10_nominal_ZodiCorrected.fits'
fi


# ------------------------------------
# ------------------------------------
# Planck 2015 release

if [[ -d 'PR2-2015' ]]
then
    echo 'You have directory PR2-2015'
else
    echo 'You have not directory PR2-2015, we will create it. '
    mkdir PR2-2015
fi

# ------------------------------------
# Foreground: synchrotron, free-free, CO21
if [[ -f PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_Synchrotron-commander_0256_R2.00.fits'
else
    echo 'You have not file COM_CompMap_Synchrotron-commander_0256_R2.00.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_Synchrotron-commander_0256_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Synchrotron-commander_0256_R2.00.fits'
fi

if [[ -f PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_freefree-commander_0256_R2.00.fits'
else
    echo 'You have not file COM_CompMap_freefree-commander_0256_R2.00.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_freefree-commander_0256_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_freefree-commander_0256_R2.00.fits'
fi

if [[ -f PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_CO21-commander_2048_R2.00.fits'
else
    echo 'You have not file COM_CompMap_CO21-commander_2048_R2.00.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_CO21-commander_2048_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_CO21-commander_2048_R2.00.fits'
fi

if [[ -f PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits ]]
then
    echo 'You have file COM_CompMap_xline-commander_0256_R2.00.fits'
else
    echo 'You have not file COM_CompMap_xline-commander_0256_R2.00.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_xline-commander_0256_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_xline-commander_0256_R2.00.fits'
fi

# ------------------------------------
# Masks
if [[ -f PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits ]]
then
    echo 'You have file HFI_Mask_GalPlane-apo0_2048_R2.00.fits'
else
    echo 'You have not file HFI_Mask_GalPlane-apo0_2048_R2.00.fits, we will download it. '
    wget -cO PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_GalPlane-apo0_2048_R2.00.fits'
fi

if [[ -f PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits ]]
then
    echo 'You have file HFI_Mask_PointSrc_2048_R2.00.fits'
else
    echo 'You have not file HFI_Mask_PointSrc_2048_R2.00.fits, we will download it. '
    wget -cO PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/ancillary-data/masks/HFI_Mask_PointSrc_2048_R2.00.fits'
fi

# ------------------------------------
# Thermal dust model (Planck 2015)
if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
else
    echo 'You have not file COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Opacity_2048_R2.01.fits'
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
else
    echo 'You have not file COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Spectral-Index_2048_R2.01.fits'
fi

if [[ -f PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits ]]
then
    echo 'You have file COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
else
    echo 'You have not file COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/all-sky-maps/maps/component-maps/foregrounds/COM_CompMap_Dust-GNILC-Model-Temperature_2048_R2.01.fits'
fi

# ------------------------------------
# Compact source catalog
if [[ -f PR2-2015/COM_PCCS_030_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_030_R2.04.fits'
else
    echo 'You have not file COM_PCCS_030_R2.04.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_030_R2.04.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_030_R2.04.fits'
fi

if [[ -f PR2-2015/COM_PCCS_044_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_044_R2.04.fits'
else
    echo 'You have not file COM_PCCS_044_R2.04.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_044_R2.04.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_044_R2.04.fits'
fi

if [[ -f PR2-2015/COM_PCCS_070_R2.04.fits ]]
then
    echo 'You have file COM_PCCS_070_R2.04.fits'
else
    echo 'You have not file COM_PCCS_070_R2.04.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_070_R2.04.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_070_R2.04.fits'
fi

if [[ -f PR2-2015/COM_PCCS_100-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_100-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_100-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_100-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_100-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_100_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_100_R2.01.fits'
else
    echo 'You have not file COM_PCCS_100_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_100_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_100_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_143-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_143-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_143-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_143-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_143-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_143_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_143_R2.01.fits'
else
    echo 'You have not file COM_PCCS_143_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_143_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_143_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_217-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_217-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_217-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_217-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_217-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_217_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_217_R2.01.fits'
else
    echo 'You have not file COM_PCCS_217_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_217_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_217_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_353-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_353-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_353-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_353-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_353-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_353_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_353_R2.01.fits'
else
    echo 'You have not file COM_PCCS_353_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_353_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_353_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_545-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_545-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_545-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_545-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_545-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_545_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_545_R2.01.fits'
else
    echo 'You have not file COM_PCCS_545_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_545_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_545_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_857-excluded_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_857-excluded_R2.01.fits'
else
    echo 'You have not file COM_PCCS_857-excluded_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_857-excluded_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_857-excluded_R2.01.fits'
fi

if [[ -f PR2-2015/COM_PCCS_857_R2.01.fits ]]
then
    echo 'You have file COM_PCCS_857_R2.01.fits'
else
    echo 'You have not file COM_PCCS_857_R2.01.fits, we will download it. '
    wget -cO PR2-2015/COM_PCCS_857_R2.01.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/COM_PCCS_857_R2.01.fits'
fi


# ------------------------------------
# ------------------------------------
# Planck 2018 release

if [[ -d 'PR3-2018' ]]
then
    echo 'You have directory PR3-2018'
else
    echo 'You have not directory PR3-2018, we will create it. '
    mkdir PR3-2018
fi

# ------------------------------------
# CMB SMICA
if [[ -f PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits ]]
then
    echo 'You have file COM_CMB_IQU-smica_2048_R3.00_full.fits'
else
    echo 'You have not file COM_CMB_IQU-smica_2048_R3.00_full.fits, we will download it. '
    wget -cO PR3-2018/COM_CMB_IQU-smica_2048_R3.00_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits'
fi

# ------------------------------------
# HFI sky maps
if [[ -f PR3-2018/HFI_SkyMap_100_2048_R3.01_full.fits ]]
then
   echo 'You have file HFI_SkyMap_100_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_100_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_100_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_100_2048_R3.01_full.fits'
fi

if [[ -f PR3-2018/HFI_SkyMap_143_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_143_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_143_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_143_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_143_2048_R3.01_full.fits'
fi

if [[ -f PR3-2018/HFI_SkyMap_217_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_217_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_217_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_217_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_217_2048_R3.01_full.fits'
fi

if [[ -f PR3-2018/HFI_SkyMap_353_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_353_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_353_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_353_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_353_2048_R3.01_full.fits'
fi

if [[ -f PR3-2018/HFI_SkyMap_545_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_545_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_545_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_545_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_545_2048_R3.01_full.fits'
fi

if [[ -f PR3-2018/HFI_SkyMap_857_2048_R3.01_full.fits ]]
then
    echo 'You have file HFI_SkyMap_857_2048_R3.01_full.fits'
else
    echo 'You have not file HFI_SkyMap_857_2048_R3.01_full.fits, we will download it. '
    wget -cO PR3-2018/HFI_SkyMap_857_2048_R3.01_full.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/HFI_SkyMap_857_2048_R3.01_full.fits'
fi

# ------------------------------------
# Transmission spectra
if [[ -f PR3-2018/HFI_RIMO_R3.00.fits ]]
then
    echo 'You have file HFI_RIMO_R3.00.fits'
else
    echo 'You have not file HFI_RIMO_R3.00.fits, we will download it. '
    wget -cO PR3-2018/HFI_RIMO_R3.00.fits 'https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/HFI_RIMO_R3.00.fits'
fi


# ------------------------------------
# ------------------------------------
# Thermal dust model from Melis et al. A&A 623, A21 (2019)
if [[ -d 'Melis' ]]
then
    echo 'You have directory Melis'
else
    echo 'You have not directory Melis, we will create it. '
    mkdir Melis
fi

# ------------------------------------
# Thermal dust model from Melis et al. A&A 623, A21 (2019)
if [[ -f Melis/beta.fits ]]
then
    echo 'You have file beta.fits'
else
    echo 'You have not file beta.fits, we will download it. '
    wget -cO Melis/beta.fits 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/beta.fits'
fi

if [[ -f Melis/tau.fits ]]
then
    echo 'You have file tau.fits'
else
    echo 'You have not file tau.fits, we will download it. '
    wget -cO Melis/tau.fits 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/tau.fits'
fi

if [[ -f Melis/temp.fits ]]
then
    echo 'You have file temp.fits'
else
    echo 'You have not file temp.fits, we will download it. '
    wget -cO Melis/temp.fits 'https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/fits/temp.fits'
fi
