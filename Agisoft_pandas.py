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
    import pandas as pd
    imgCount_df = pd.DataFrame({"Date": 0,
               "Medium" : 0,
               "MedPts": 0,
               "MedTie": 0,
               "MedptsTie": 0,
               "MedExtra": 0,}, index=[0])
    
    print("Pandas has been installed")
except:
    print("Failed to import Pandas")

