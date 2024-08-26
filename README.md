# 0. Introduction
This project is used to test the goodness of the single-component modified blackbody thermal dust emission model in CMB experiments. 

More technical details can be found in the following paper: 

# 1. Dependencies and installation
## 1.1 Dependencies
The main dependencies and their versions are listed below: 

`ubuntu 20.04`
`Python 3.8.10`
`astropy 5.2.2`
`healpy 1.16.5`
`matplotlib 3.7.5`
`numpy 1.24.4`
`pip 20.0.2`
`scipy 1.10.1`

The computer configuration should be no less than 12 cores and 64GB of memory. 

## Installation

```
git clone https://github.com/Jia-Rui-Li/Thermal-dust-components.git
cd Thermal-dust-components
bash configure.sh
```

The default behavior of `configure.sh` is a full installation of python packages and promordial .fits file. 

### Installation of python packages
```
sudo apt-get update
sudo apt-get install python3
sudo apt-get install python3-pip
pip3 install numpy
pip3 install matplotlib
pip3 install healpy
pip3 install scipy
pip3 install astropy
```

### Download .fits files 
All of the .fits files about Planck are downloaded from [Caltech database](https://irsa.ipac.caltech.edu/data/Planck/), 
saved in directory /PR1-2013/, /PR2-2015/, /PR3-2018/. 
Files about [M19](https://doi.org/10.1051/0004-6361/201834394) are downloaded from [https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/](https://cdsarc.cds.unistra.fr/ftp/J/A+A/623/A21/), 
saved in directory /Melis/. 

# 2. Running
```
bash run.sh
```
## Step-1
```
python Inpaint_Rebeam_HFI_maps.py
```

1. Masking HFI maps with point source masks. 
2. Iterating 2000 times to in-paint the holes after masking. 

## Step-2
```
python Dust_Planck_data_model.py
```
1. Isolating dust data maps from Planck HFI maps. 
2. Calculating dust model maps from M13, M15, and M19. 

## Step-3
```
python Dust_Planck_data_model_No_color_correction.py
```
Same with Step-2, but without color correction. 

## Step-4
```
python Statistics.py
```
Calculating statistical quantities, including linear ratio R (for model), R' (for data), and correlation coefficient C' (for data). 

## Step-5
```
python Plots.py
```
Plotting. 

# 3. Citation

