# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 15:23:46 2020
Test that pandas properly installed in Agisoft python environment
@author: Thomas Van Der Weide
"""

import os
import re
import glob
try:
    import Metashape as PhotoScan
except ImportError:
    import PhotoScan
# import split_in_chunks_python

try:
    # import pandas as pd
    # imgCount_df = pd.DataFrame({"Date": 0,
    #            "Medium" : 0,
    #            "MedPts": 0,
    #            "MedTie": 0,
    #            "MedptsTie": 0,
    #            "MedExtra": 0,}, index=[0])
    
    # print("Pandas has been installed")
    import csv
        mainDict = {"Date": 0,
                "Accuracy" : 0,
                "Key_Limit" : 0,
                "Tie_Limit" : 0,
                "Quality" : 0,
                "Filter" : 0,
                "Aligned" : 0,
                "Total" : 0,
                "Img_Quality_Avg" : 0,
                "RU_Thresh" : 0,
                "PA_Thresh" : 0,
                "RE_Thresh" : 0,
                "Tie_ptsBefore" : 0,
                "Marker_errBefore" : 0,
                "Tie_ptsAfter" : 0,
                "Marker_errAfter" : 0}
    print("CSV can successfully be called")
except:
    print("Failed to use CSV library")

