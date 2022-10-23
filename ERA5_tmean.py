'''
#reference: Mateus, P., Catalão, J., Mendes, V. B., & Nico, G. (2020). 
#   An ERA5-Based Hourly Global Pressure and Temperature (HGPT) Model. 
#   Remote Sensing, 12(7). doi:10.3390/rs12071098

#   Processing by Miniconda
#   How to compile on Notepad++
        : %windir%\System32\cmd.exe "/K" C:\Miniconda3\Scripts\activate.bat C:\Miniconda3 && C:\Miniconda3\python.exe -i "$(FULL_CURRENT_PATH)"
'''
#import pandas as pd
#from scipy.stats.stats import pearsonr

# Miniconda 

from math import fabs,floor
from time import perf_counter
import fnmatch
import os,sys
#import matplotlib.pyplot as plt
import pygrib
#from shapely.geometry import shape,Point
import json
import datetime
import pandas as pd
from datetime import timedelta, date
import datetime
dir_app = os.path.dirname(sys.argv[0])

'''
  Grib specific_humidity information
    + r,c : 65,37 start at 0,0
    + A number of layers :
        = 24 hr x 37 pressure level 
        when
            hr: 0, 100, 200, 300,..., 2300
            pressure 37 level: 
                1   2   3   5   7
               10  20  30  50  70
              100 125 150 175 200
              225 250 300 350 400
              450 500 550 600 650
              700 750 775 800 825
              850 875 900 925 950
              975 1000
        = 888 layers start at 1
'''

#grib = pygrib.open(r'C:\Users\Chaiyut\Downloads\temperature_20191230_37levels_24hr.grib')
#for g in grib:
#    v = g.values
#    print(g.longitudeOfFirstGridPoint)
def mkdir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
        
def daterange(start_date='%Y-%m-%d', end_date='%Y-%m-%d'):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
        
def date2doy(date='%Y-%m-%d %H:%M:%S'):
	# input y,m,d,hh,mm,ss : integer
	# output doy decimal : float
	datetimeformat = '%Y-%m-%d %H:%M:%S'
	dd = datetime.datetime.strptime(date, datetimeformat)
	doy = dd.strftime('%j')
	dec = (dd.hour + dd.minute/60.0 + dd.second/3600.0)/24.0
	doydec = int(doy) + dec
	DOY = doydec
	return DOY

def oro_dictionary(gribfile,rcthai_dic):
    g=9.80665
    dic = {};
    grbs = pygrib.open(gribfile)
    grb = grbs[1]
    return grb.values/g

def surtemp_dictionary(gribfile):
    dic = {};
    grbs = pygrib.open(gribfile)
    for grb in grbs:
        dd = grb.date
        tt = grb.time
        key =str(dd) + '-' + str(tt)
        dic[key] = grb.values
    return dic
def geo_dictionary(gribfile):
    g=9.80665;
    dic = {};
    grbs = pygrib.open(gribfile)
    for grb in grbs:
        r,c = grb.values.shape
        pres = grb.levelist
        if pres==1000:
            key =str(pres)
            dic[key] = grb.values/g 
    return dic
def temp_dictionary(gribfile):
    dic = {};
    grbs = pygrib.open(gribfile)
    for grb in grbs:
        pres = grb.levelist
        dd = grb.date
        tt = grb.time
        key =str(pres)+'-'+str(dd) + '-'+ str(tt)
        dic[key] = grb.values
    return dic
    
def rc_thai():
    rc_thai = r'E:\ECMWF\rc_thai.txt'
    dic = {};
    with open(rc_thai,'r') as f:
        for line in f:
            dic[line.rstrip('\n')] =''
        f.close()
    return dic
    
def pointInThai(rc,rcthai_dic):
    inThai = False
    try:
        rcthai_dic[rc]
        inThai = True
    except:
        inThai = False
    return inThai
def testRead():
    grbs_tem = pygrib.open(r'E:\ECMWF\Pressure_Level\temperature\2014\20140101_0023hr.grib')
    grbs_hum = pygrib.open(r'E:\ECMWF\Pressure_Level\specific_humidity\2014\20140101_0023hr.grib')
    for grb in grbs_hum:
        
        pres = grb.level
        dd = grb.date
        tt = grb.time
        hr = str(int(tt/100))                                
        hum_data = f'{grb.level}:{grb.date}:{grb.time/100:.0f}:{grb.values[0][0]}'
        print(hum_data)
        
        grb_temp = grbs_tem.readline()                       
        temp_data = f'{grb_temp.level}:{grb_temp.date}:{grb_temp.time/100:.0f}:{grb_temp.values[0][0]}'  
        print(temp_data+'\n')
def tm(workdir,name,oro_dic,rcthai_dic,outfld):
    mkdir(outfld)
    outfile = os.path.join(outfld,name[0:8] + '.csv')
    if os.path.exists(outfile):
        print('File is Exist: ' + name[0:8] + '.csv')
        return                
    print(name)
    tic = perf_counter()
    with open(outfile,'w') as f:
        surtemp_dic = surtemp_dictionary(os.path.join(r'E:\ECMWF\Singel_Level\2m_temperature',name[0:4]+'\\' + name[0:8] + '.grib'))
        #geo_dic = geo_dictionary(os.path.join(r'E:\ECMWF\Pressure_Level\geopotential',name[0:4] + '\\' + name))
        temp_dic = temp_dictionary(os.path.join(r'E:\ECMWF\Pressure_Level\temperature',name[0:4] + '\\' + name))
        f.write('doy,subdoy,r,c,lat,lon,ee,tm,t_hPa,ts,tm_h,ts_h,date,hr\n')
        grbs = pygrib.open(os.path.join(workdir, name))
        r,c = grbs[1].values.shape
        for i in range(0,r,1):
            for j in range(0,c,1):
                rc = str(i)+'-'+str(j)
                if pointInThai(rc,rcthai_dic):
                    grbs.seek(0)
                    SUM_U, SUM_L = 0.0, 0.0
                    underSurf=False
                    #print(rc);
                    for grb in grbs:
                        # 37 level
                        pres = grb.levelist
                        tt = grb.time
                        dd = str(grb.date)
                        hum = grb.values[i][j]
                        ee = int(pres)*hum/(hum-(hum-1)*0.622)
                        temp_values = temp_dic[str(pres)+'-'+str(dd)+'-'+str(tt)]
                        temp = temp_values[i][j]
                        SUM_U = SUM_U + ee/temp
                        SUM_L = SUM_L + ee/(temp*temp)
                        if pres==1000:
                            lats,lons = grb.latlons()
                            hr = str(int(tt/100))                                
                            ymd_hr = dd[0:4] + '-' + dd[4:6] + '-' + dd[6:8] + ' ' + hr + ':00:00'
                            doyc = date2doy(ymd_hr)  
                            subdoy = doyc-floor(doyc)   
                            
                            ts_h = oro_dic[i][j]
                            tm_h = 0;#(geo_dic[str(pres)])[i][j]
                            tmean = SUM_U/SUM_L
                            ts = (surtemp_dic[str(dd) + '-' + str(tt)])[i][j]
                            lat,lon = lats[i][j],lons[i][j]
                            result  = '{:.7f}'.format(doyc)+','+'{:.7f}'.format(subdoy)
                            result += ','+str(i)+','+str(j)+','+str(lat)+','+str(lon)
                            result += ','+'{:.2f}'.format(ee)
                            result += ','+'{:.2f}'.format(tmean)+','+'{:.2f}'.format(temp)+','+'{:.2f}'.format(ts)
                            result += ','+'{:.0f}'.format(tm_h)+','+'{:.0f}'.format(ts_h)
                            result += ','+str(dd)+','+ hr
                            f.write(result + '\n')
                            #print(result)
                            SUM_U, SUM_L = 0.0, 0.0
                            #print(result);
                            #underSurf=False

        grbs.close()
        f.close()
    toc = perf_counter()
    print(str(int((toc-tic)/60.0)) + ' min')
    
def cal_tm_daily(hm_work,oro_dic,rcthai_dic,outfld):
    for root, dirs, files in os.walk(hm_work, topdown=True):
        for name in files:
            tm(root,name,oro_dic,rcthai_dic,outfld)



rcthai_dic = rc_thai()
oro_dic = oro_dictionary(r'E:\ECMWF\Singel_Level\orography\2019\20190101.grib',rcthai_dic)
hm_work = r'E:\ECMWF\Pressure_Level\specific_humidity\2013\6'
outfld = r'E:\ECMWF\tmean\2013'
tm(hm_work,name,oro_dic,rcthai_dic,outfld) 
