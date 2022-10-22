# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 22:38:16 2021

@author: chaiyut
"""
# Download wget from this url 'https://eternallybored.org/misc/wget/'
import os
import time
import requests
import time
#from furl import furl
link = r'D:\EarthData\link_download\standard\subset_AIRS2RET_7.0_2016.txt'
wget = r'D:\EarthData\software\wget.exe'
outdir = r'D:\EarthData\hdf_standard\2016'

# Download all
command = f' -nc -P {outdir} --load-cookies C:\AD\.urs_cookies --save-cookies C:\AD\.urs_cookies --auth-no-challenge=on --keep-session-cookies --user=skycaptain --password=I=ypp6mT --content-disposition -i '
os.system(wget + command  + link)
