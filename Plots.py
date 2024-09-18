#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Codes for plots in the article

import healpy
import os
import numpy
import scipy
import astropy.io.fits
import time
import shutil
import sys
import matplotlib.pyplot
import multiprocessing

root_dir = os.path.abspath('.')

if os.path.exists(root_dir+"/figure/"):
    shutil.rmtree(root_dir+"/figure/")
os.mkdir(root_dir+"/figure/")

# dust maps at HFI frequencies from Planck HFI maps
def thermal_dust_map(freq_str, unit):
    if unit == "muK_CMB":
        unit_str = "$\mu\mathrm{K_{CMB}}$"
    elif unit == "MJy_sr":
        unit_str = "$\mathrm{MJy\,sr^{-1}}$"

    if freq_str == "100":
        MAX = 200
    elif freq_str == "143":
        MAX = 300
    elif freq_str == "217":
        MAX = 800
    elif freq_str == "353":
        MAX = 5000
    elif freq_str == "545":
        MAX = 7
    elif freq_str == "857":
        MAX = 20
    dust_map = healpy.read_map(root_dir+"/results/dust_data_"+unit+"_"+freq_str+"GHz.fits")
    matplotlib.pyplot.figure(dpi=300, figsize=(3,3))
    healpy.projview(dust_map, cmap="jet", unit=unit_str, norm="none", min=0, max=MAX, extend="neither", hold=True)
    matplotlib.pyplot.title("Dust data map at "+freq_str+"GHz", fontsize=12)
    matplotlib.pyplot.savefig(root_dir+"/figure/thermal_dust_component_in_"+freq_str+"GHz.pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

    dust_map = healpy.smoothing(dust_map, fwhm = 1*numpy.pi/180, pol=False)
    MAX = numpy.percentile(dust_map, 70)
    matplotlib.pyplot.figure(dpi=300, figsize=(3,3))
    healpy.projview(dust_map, cmap="jet", unit=unit_str, norm="none", min=0, max=MAX, format = "%.4g", extend="neither", hold=True)
    matplotlib.pyplot.title("Smoothed dust data map at "+freq_str+"GHz", fontsize=12)
    matplotlib.pyplot.savefig(root_dir+"/figure/smooth_thermal_dust_component_in_"+freq_str+"GHz.pdf", bbox_inches="tight")
    matplotlib.pyplot.close()


# Sky maps for R(model), R'(data), C(model), C'(data), where C is correlation
def map_plot(model_name, freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask): 
    # name = 1: 143 / 100 GHz
    # name = 2: 353 / 217 GHz
    # name = 3: 857 / 545 GHz
    if freq_name == "1":
        mosaic_model = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_100GHz_143GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_100GHz_143GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "100 - 143 GHz pair"
        min1 = 1.5; max1 = 2.5; min2 = 0.95; max2 = 1; min3 = 1.5; max3 = 2.5; min4 = 0.8; max4 = 1;
    elif freq_name == "2":
        mosaic_model = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_217GHz_353GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_217GHz_353GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "217 - 353 GHz pair"
        min1 = 6; max1 = 9; min2 = 0.95; max2 = 1; min3 = 6; max3 = 9; min4 = 0.95; max4 = 1;
    elif freq_name == "3":
        mosaic_model = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_545GHz_857GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_545GHz_857GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "545 - 857 GHz pair"
        min1 = 2.3; max1 = 3.3; min2 = 0.95; max2 = 1; min3 = 2.3; max3 = 3.3; min4 = 0.99; max4 = 1;

    # R (model)
    matplotlib.pyplot.figure(dpi=300, figsize=(5,5))
    healpy.mollview(mosaic_model[0], nlocs=3, nest=False, min=min1, max=max1, cmap="jet", hold=True)
    matplotlib.pyplot.title("$R$ of "+str_name, fontsize=12.5)
    matplotlib.pyplot.savefig(root_dir+"/figure/R_"+model_name+"_color_correction_"+color_correction+"_"+freq_name+"_galactic_mask_"+galac_mask+"_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

    # R' (data)
    matplotlib.pyplot.figure(dpi=300, figsize=(5,5))
    healpy.mollview(mosaic_data[0], nlocs=3, nest=False, min=min3, max=max3, cmap="jet", hold=True)
    matplotlib.pyplot.title("$R'$ of "+str_name, fontsize=12.5)
    matplotlib.pyplot.savefig(root_dir+"/figure/R_data_color_correction_"+color_correction+"_"+freq_name+"_galactic_mask_"+galac_mask+"_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

    # Correlation between two neighbouring frequency channels (data)
    matplotlib.pyplot.figure(dpi=300, figsize=(5,5))
    healpy.mollview(mosaic_data[1], nest=False, min=min4, max=max4, cmap="jet", hold=True)
    matplotlib.pyplot.title("Correlation between "+str_name, fontsize=12.5)
    matplotlib.pyplot.savefig(root_dir+"/figure/C_data_color_correction_"+color_correction+"_"+freq_name+"_galactic_mask_"+galac_mask+"_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()
    return 1

# Valid regions
def region(freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask): 
    if freq_name == "1":
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_100GHz_143GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "100 - 143 GHz pair"
    elif freq_name == "2":
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_217GHz_353GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "217 - 353 GHz pair"
    elif freq_name == "3":
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_545GHz_857GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "545 - 857 GHz pair"

    # total_pixel_list is index of the points without masked (valid region)
    total_pixel_list = numpy.where(mosaic_data[1]>0)[0]
    # pixel_list is index of the points plotted in the scatter plot with C'>0.95 (reliable region, in red)
    pixel_list = numpy.where(mosaic_data[2]==1)[0]
    # fraction_reliable = reliable region / valid region
    fraction_reliable = numpy.array(pixel_list).shape[0]/numpy.array(total_pixel_list).shape[0]
    matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
    healpy.mollview(mosaic_data[2], cmap="jet", min=0, max=2, cbar=False, hold=True)
    matplotlib.pyplot.title("Reliable region of "+str_name+"\n $f_\mathrm{rel} = $"+"{:.1f}".format(100*fraction_reliable)+"%", fontsize=18)
    matplotlib.pyplot.savefig(root_dir+"/figure/Region_color_correction_"+color_correction+"_"+freq_name+"_galactic_mask_"+galac_mask+"_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".pdf", bbox_inches="tight")
    matplotlib.pyplot.close()

# Scatter plots for dust (data) and the sky map for the points in the scatter plots
def scatter_plot(model_name, freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask):
    if freq_name == "1":
        mosaic_model= numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_100GHz_143GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_100GHz_143GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "100 - 143 GHz pair"
    elif freq_name == "2":
        mosaic_model= numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_217GHz_353GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_217GHz_353GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "217 - 353 GHz pair"
    elif freq_name == "3":
        mosaic_model= numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_"+model_name+"_545GHz_857GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        mosaic_data = numpy.load(root_dir+"/results/scatter_points_color_correction_"+color_correction+"_data_545GHz_857GHz_galactic_mask_"+galac_mask+"%_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".npy")
        str_name = "545 - 857 GHz pair"

    pixel_list = numpy.where(mosaic_data[2]==1)[0]
    # R (model) from the points in reliable region
    model= mosaic_model[0][pixel_list]
    # R' (data) from the points in reliable region
    data = mosaic_data[0][pixel_list]
    # Linear fitting of R' vs R
    slope, intercept = numpy.polyfit(model, data, 1)
    x = numpy.arange(0,10,0.1)
    data_fit = slope * x + intercept
    line_name = "y = "+"{:+.3f}".format(slope)+" x "+"{:+.3f}".format(intercept)

    matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
    matplotlib.pyplot.scatter(model, data, s=2, marker="+", color="black")
    matplotlib.pyplot.plot(x, x, label="y = x", linewidth=3)
    matplotlib.pyplot.plot(x, data_fit, label="Linear fitting: ("+line_name+")", linewidth=3)
    matplotlib.pyplot.xlim(xmin=min(model), xmax=max(model))
    matplotlib.pyplot.ylim(ymin=min(data),  ymax=min(data)+1.20*(max(data)-min(data)))
    matplotlib.pyplot.xlabel(r"$R$ (model)", fontsize=17)
    matplotlib.pyplot.ylabel(r"$R'$ (data)", fontsize=17)
    matplotlib.pyplot.legend(fontsize=17, loc="upper left")
    matplotlib.pyplot.title(str_name, fontsize=21)
    matplotlib.pyplot.savefig(root_dir+"/figure/Scatter_"+model_name+"_color_correction_"+color_correction+"_"+freq_name+"_galactic_mask_"+galac_mask+"_smooth_degree_"+"{:01d}".format(int(smooth_degree))+"_disk_degree_"+"{:01d}".format(int(disk_degree))+"_low_galac_mask_"+low_mask+"_zodiacal_mask_"+zodiacal_mask+".png", bbox_inches="tight")
    matplotlib.pyplot.close()


print("Plots: ")
#--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*--*
# Maps of thermal dust emission from Planck data
arguments = [("100", "muK_CMB"), ("143", "muK_CMB"), ("217", "muK_CMB"), ("353", "muK_CMB"), ("545", "MJy_sr"), ("857", "MJy_sr")] 
multiprocessing.Pool(processes=6).starmap(thermal_dust_map, arguments)

source_point_mask = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_PointSrc_2048_R2.00.fits", h=False, hdu=1, field=5)
compact_mask = healpy.read_map(root_dir+"/results/compact_source_mask.fits", h=False, field=0)
galac_plane_mask_60 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=2)
galac_plane_mask_70 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=3)
galac_plane_mask_80 = healpy.read_map(root_dir+"/PR2-2015/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", h=False, field=4)
low_galac_mask = numpy.zeros(12*2048**2, dtype=int)
low_galac_list = healpy.query_strip(2048, 2/3*numpy.pi, 1/3*numpy.pi, inclusive=True)
low_galac_mask[low_galac_list] = 1

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where(source_point_mask > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("Point source mask at 857 GHz", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_point_source_857.pdf", bbox_inches="tight")
matplotlib.pyplot.show()

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where(compact_mask > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("$M_\\mathrm{comp}$", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_compact_source.pdf", bbox_inches="tight")
matplotlib.pyplot.show()

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where((compact_mask*galac_plane_mask_60) > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("$M_\\mathrm{comp} \\times M_{60}$", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_galac_plane_60.pdf", bbox_inches="tight")
matplotlib.pyplot.show()

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where((compact_mask*galac_plane_mask_70) > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("$M_\\mathrm{comp} \\times M_{70}$", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_galac_plane_70.pdf", bbox_inches="tight")
matplotlib.pyplot.show()

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where((compact_mask*galac_plane_mask_80) > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("$M_\\mathrm{comp} \\times M_{80}$", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_galac_plane_80.pdf", bbox_inches="tight")
matplotlib.pyplot.show()

mask = numpy.zeros(12*2048**2, dtype=int)
mask = mask - 1.6375e+30
valid= numpy.where((compact_mask*galac_plane_mask_80*low_galac_mask) > 0.9)[0]
mask[valid] = 1
matplotlib.pyplot.figure(dpi=300, figsize=(7,7))
healpy.mollview(mask, norm="none", cmap="jet", min=0, max=2, cbar=None, hold=True)
matplotlib.pyplot.title("$M_\\mathrm{comp} \\times M_{80} \\times M_\\mathrm{lat}$", fontsize=20)
matplotlib.pyplot.savefig(root_dir+"/figure/Mask_low_galac.pdf", bbox_inches="tight")
matplotlib.pyplot.show()


# map_plot(model_name, freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
# region(freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
# scatter_plot(model_name, freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

nside1 = 64
nside2 = 2048

color_correction = "yes"
smooth_degree = 2
disk_degree = 6
galac_mask = "80"
low_mask = "no"
zodiacal_mask = "no"
for freq_name in ["1", "2", "3"]: 
    map_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    region(freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

color_correction = "yes"
smooth_degree = 2
disk_degree = 6
galac_mask = "80"
low_mask = "yes"
zodiacal_mask = "no"
for freq_name in ["1", "2", "3"]: 
    region(freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

color_correction = "yes"
galac_mask = "80"
low_mask = "no"
zodiacal_mask = "no"
for smooth_degree in [1, 2, 3]:
    for disk_degree in [5, 6, 7]:
        for freq_name in ["1", "2", "3"]: 
            # without excluding the region with latitude < 30 degree, with color correction, with galactic plane mask for 80% coverage
            scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

color_correction = "yes"
smooth_degree = 2
disk_degree = 6
low_mask = "no"
zodiacal_mask = "no"
for galac_mask in ["60", "70", "80"]:
    for freq_name in ["1", "2", "3"]:
        # without excluding the region with latitude < 30 degree, with color correction, with galactic plane mask for 60/70/80% coverage
        scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

color_correction = "yes"
smooth_degree = 2
disk_degree = 6
galac_mask = "80"
low_mask = "no"
zodiacal_mask = "no"
for freq_name in ["1", "2", "3"]:
    # without excluding the region with latitude < 30 degree, without color correction, with galactic plane mask for 80% coverage
    scatter_plot("model_2013", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2019", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

color_correction = "no"
smooth_degree = 2
disk_degree = 6
galac_mask = "80"
low_mask = "no"
zodiacal_mask = "no"
for freq_name in ["1", "2", "3"]:
    # without excluding the region with latitude < 30 degree, without color correction, with galactic plane mask for 80% coverage
    scatter_plot("model_2013", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2019", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)

nside1 = 64
nside2 = 2048
color_correction = "yes"
smooth_degree = 2
disk_degree = 6
galac_mask = "80"
low_mask = "no"
zodiacal_mask = "yes"
for freq_name in ["1", "2", "3"]:
# without excluding the region with latitude < 30 degree, with color correction, with galactic plane mask for 80% coverage
    region(freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
    scatter_plot("model_2015", freq_name, color_correction, smooth_degree, disk_degree, galac_mask, low_mask, zodiacal_mask)
