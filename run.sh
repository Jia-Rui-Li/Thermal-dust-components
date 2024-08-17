#!/bin/bash
python -u Inpaint_Rebeam_HFI_maps.py 2000
python -u Dust_Planck_data_model.py 
python -u Dust_Planck_data_model_No_color_correction.py
python -u Statistics.py
python -u Plots.py