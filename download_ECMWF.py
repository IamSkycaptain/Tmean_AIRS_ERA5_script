import cdsapi
import os
from datetime import timedelta, date
import datetime
import time
from calendar import monthrange

#Ref:
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview
#
import requests
def line_notify(yy):
    url = 'https://notify-api.line.me/api/notify'
    token = 'xxxxxx'
    headers = {'content-type':'application/x-www-form-urlencoded','Authorization':'Bearer '+token}

    msg = 'Load total rain column Complete: ' + yy
    r = requests.post(url, headers=headers, data = {'message':msg})
    return
    
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
def daterange(start_date='%Y-%m-%d', end_date='%Y-%m-%d'):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)
def mkdir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
        
def load_2mterperature_monthly(year):
	c = cdsapi.Client()
	for mm in range(1,13,1):
		outfld = os.path.join(r'E:\ECMWF\Singel_Level\2m_temperature',str(year))
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, str(year) + ('0' + str(mm))[-2:]+'.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve('reanalysis-era5-single-levels', {
				'variable'      : '2m_temperature',
				'product_type'  : 'reanalysis',
                'year': str(year),
                'month': str(mm),
                'day': [
                    '01', '02', '03',
                    '04', '05', '06',
                    '07', '08', '09',
                    '10', '11', '12',
                    '13', '14', '15',
                    '16', '17', '18',
                    '19', '20', '21',
                    '22', '23', '24',
                    '25', '26', '27',
                    '28', '29', '30',
                    '31',
                ],                
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
				'format': 'grib',
			}, outfile)
		time.sleep(1)
        
def load_2mterperature(year):
    c = cdsapi.Client()
    for single_date in daterange(date(year, 1, 1), date(year+1, 1, 1)):
        ymd = single_date.strftime("%Y-%m-%d")
        yy,mm,dd = ymd.split('-')
        outfld = os.path.join(r'D:\ECMWF\Singel_Level\2m_temperature',str(year))
        outfile = os.path.join(outfld, ymd.replace('-','') + '.grib')
        if not os.path.exists(outfld):
            mkdir(outfld)

        if os.path.exists(outfile):
            continue
        print(outfile)
        c.retrieve('reanalysis-era5-single-levels', {
				'variable'      : '2m_temperature',
				'product_type'  : 'reanalysis',
				'year': yy,
				'month':mm,
				'day': dd,
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
				'format': 'grib',
			}, outfile)
            
def load_Orography(year):
    c = cdsapi.Client()
    for single_date in daterange(date(year, 1, 1), date(year+1, 1, 1)):
        ymd = single_date.strftime("%Y-%m-%d")
        yy,mm,dd = ymd.split('-')
        outfld = os.path.join(r'C:\ECMWF\Singel_Level\orography',str(year))
        outfile = os.path.join(outfld, ymd.replace('-','') + '.grib')
        if not os.path.exists(outfld):
            mkdir(outfld)

        if os.path.exists(outfile):
            continue
        print(outfile)
        c.retrieve('reanalysis-era5-single-levels', {
				'variable'      : 'orography',
				'product_type'  : 'reanalysis',
				'year': yy,
				'month':mm,
				'day': dd,
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time':'12:00',
				'format': 'grib'
			}, outfile)
def load_surfacepressure(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1, 1, 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		outfld = r'D:\ECMWF\Singel_Level\surface_pressure' + '\\' + str(year)
		outfile = os.path.join(outfld, ymd.replace('-','') + '.grib')
		if not os.path.exists(outfld):
			mkdir(outfld)
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve('reanalysis-era5-single-levels', {
				'variable'      : 'surface_pressure',
				'product_type'  : 'reanalysis',
				'year': yy,
				'month':mm,
				'day': dd,
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
				'format': 'grib',
			}, outfile)
def load_Surface_Pressure(year):
    c = cdsapi.Client()
    outfile = os.path.join(r'D:\ECMWF\Singel_Level\surface_pressure',str(year) + '.grib')
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': 'surface_pressure',
            'area': [21, 97, 5, 106], # North, West, South, East. Default: global
            'year': '2020',
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
        },
        outfile)
def load_specific_humidity_monthly(year):
	c = cdsapi.Client()
	for mm in range(1,13,1):
		outfld = os.path.join(r'E:\ECMWF\Pressure_Level\specific_humidity',str(year))
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, str(year) + ('0' + str(mm))[-2:]+'.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'specific_humidity',
			'pressure_level': [
				'1', '2', '3',
				'5', '7', '10',
				'20', '30', '50',
				'70', '100', '125',
				'150', '175', '200',
				'225', '250', '300',
				'350', '400', '450',
				'500', '550', '600',
				'650', '700', '750',
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': str(year),
			'month': str(mm),
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
			'time': [
			'00:00', '01:00', '02:00',
			'03:00', '04:00', '05:00',
			'06:00', '07:00', '08:00',
			'09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00',
			'15:00', '16:00', '17:00',
			'18:00', '19:00', '20:00',
			'21:00', '22:00', '23:00',
			],
        },outfile)
		time.sleep(1)    
def load_temperature_monthly(year):
	c = cdsapi.Client()
	for mm in range(1,13,1):
		outfld = os.path.join(r'E:\ECMWF\Pressure_Level\temperature',str(year))
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, str(year) + ('0' + str(mm))[-2:]+'.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'temperature',
			'pressure_level': [
				'1', '2', '3',
				'5', '7', '10',
				'20', '30', '50',
				'70', '100', '125',
				'150', '175', '200',
				'225', '250', '300',
				'350', '400', '450',
				'500', '550', '600',
				'650', '700', '750',
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': str(year),
			'month': str(mm),
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
			'time': [
			'00:00', '01:00', '02:00',
			'03:00', '04:00', '05:00',
			'06:00', '07:00', '08:00',
			'09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00',
			'15:00', '16:00', '17:00',
			'18:00', '19:00', '20:00',
			'21:00', '22:00', '23:00',
			],
        },outfile)
		time.sleep(1)
        
def load_geopotential_monthly(year):
	c = cdsapi.Client()
	for mm in range(1,13,1):
		outfld = os.path.join(r'E:\ECMWF\Pressure_Level\geopotential',str(year))
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, str(year) + ('0' + str(mm))[-2:]+'.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'geopotential',
			'pressure_level':'1000',
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': str(year),
			'month': str(mm),
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
			'time': [
			'00:00', '01:00', '02:00',
			'03:00', '04:00', '05:00',
			'06:00', '07:00', '08:00',
			'09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00',
			'15:00', '16:00', '17:00',
			'18:00', '19:00', '20:00',
			'21:00', '22:00', '23:00',
			],
        },outfile)
		time.sleep(1)
        
def load_specific_humidity(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1,1 , 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		#doy = ('00' + str(int(date2doy(ymd + ' 00:00:00'))))[-3:]
		outfld = os.path.join(r'D:\ECMWF\Pressure_Level\specific_humidity',str(year))
		#hour = ('0'+str(hr))[-2:]
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, ymd.replace('-','') + '_0023hr.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'specific_humidity',
			'pressure_level': [
				'1', '2', '3',
				'5', '7', '10',
				'20', '30', '50',
				'70', '100', '125',
				'150', '175', '200',
				'225', '250', '300',
				'350', '400', '450',
				'500', '550', '600',
				'650', '700', '750',
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': yy,
			'month': mm,
			'day': dd,
			#'time': hour + ':00',
			'time': [
			'00:00', '01:00', '02:00',
			'03:00', '04:00', '05:00',
			'06:00', '07:00', '08:00',
			'09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00',
			'15:00', '16:00', '17:00',
			'18:00', '19:00', '20:00',
			'21:00', '22:00', '23:00',
			],
        },outfile)
		time.sleep(1)

def load_Geopotential(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1,1 , 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		#doy = ('00' + str(int(date2doy(ymd + ' 00:00:00'))))[-3:]
		outfld = os.path.join(r'D:\ECMWF\Pressure_Level\geopotential',str(year))
		#hour = ('0'+str(hr))[-2:]
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, ymd.replace('-','') + '_0023hr.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'geopotential',
			'pressure_level': [
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': yy,
			'month': mm,
			'day': dd,
			#'time': hour + ':00',
			'time': '12:00'
		},outfile)
def load_temperature(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1,1 , 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		#doy = ('00' + str(int(date2doy(ymd + ' 00:00:00'))))[-3:]
		outfld = os.path.join(r'D:\ECMWF\Pressure_Level\temperature',str(year))
		#hour = ('0'+str(hr))[-2:]
		if not os.path.exists(outfld):
			mkdir(outfld)
		outfile = os.path.join(outfld, ymd.replace('-','') + '_0023hr.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
		'reanalysis-era5-pressure-levels',
		{
			'product_type': 'reanalysis',
			'format': 'grib',
			'variable': 'temperature',
			'pressure_level': [
				'1', '2', '3',
				'5', '7', '10',
				'20', '30', '50',
				'70', '100', '125',
				'150', '175', '200',
				'225', '250', '300',
				'350', '400', '450',
				'500', '550', '600',
				'650', '700', '750',
				'775', '800', '825',
				'850', '875', '900',
				'925', '950', '975',
				'1000',
			],
			'area': [21, 97, 5, 106], # North, West, South, East. Default: global
			'year': yy,
			'month': mm,
			'day': dd,
			#'time': hour + ':00',
			'time': [
			'00:00', '01:00', '02:00',
			'03:00', '04:00', '05:00',
			'06:00', '07:00', '08:00',
			'09:00', '10:00', '11:00',
			'12:00', '13:00', '14:00',
			'15:00', '16:00', '17:00',
			'18:00', '19:00', '20:00',
			'21:00', '22:00', '23:00',
			],
		},outfile)

def load_total_pwv_yearly(yy):
    c = cdsapi.Client()
    outfld = os.path.join(r'E:\ECMWF\Singel_Level\total_pwv',str(yy))
    mkdir(outfld)
    outfile = os.path.join(outfld, str(yy)+'0000.grib')
    if os.path.exists(outfile):
        return
    print(outfile)
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': 'total_column_water_vapour',
            'year': yy,
            'area': [21, 97, 5, 106], # North, West, South, East. Default: global
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'format': 'grib',
        },outfile)

def load_total_pwv_yearly(yy):
    c = cdsapi.Client()
    outfld = os.path.join(r'E:\ECMWF\Singel_Level\total_pwv',str(yy))
    mkdir(outfld)
    outfile = os.path.join(outfld, str(yy)+'0000.grib')
    if os.path.exists(outfile):
        return
    print(outfile)
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': 'total_column_water_vapour',
            'year': yy,
            'area': [21, 97, 5, 106], # North, West, South, East. Default: global
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'format': 'grib',
        },outfile)
            
def load_total_pwv(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1, 1, 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		outfld = r'E:\ECMWF\Singel_Level\total_pwv'+'\\' + str(year)
		mkdir(outfld)
		outfile = os.path.join(outfld, ymd.replace('-','') + '.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
			'reanalysis-era5-single-levels',
			{
				'product_type': 'reanalysis',
				'format': 'grib',
				'variable': 'total_column_water_vapour',
				'year': yy,
				'month':mm,
				'day': dd,
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
			},outfile)
def load_rain(year):
	c = cdsapi.Client()
	for single_date in daterange(date(year, 1, 1), date(year+1, 1, 1)):
		ymd = single_date.strftime("%Y-%m-%d")
		yy,mm,dd = ymd.split('-')
		outfld = r'D:\ECMWF\Singel_Level\total_rain_water' + '\\' + str(year)
		mkdir(outfld)
		outfile = os.path.join(outfld, ymd.replace('-','') + '.grib')
		if os.path.exists(outfile):
			continue
		print(outfile)
		c.retrieve(
			'reanalysis-era5-single-levels',
			{
				'product_type': 'reanalysis',
				'format': 'grib',
				'variable': 'total_column_rain_water',
				'year': yy,
				'month':mm,
				'day': dd,
				'area': [21, 97, 5, 106], # North, West, South, East. Default: global
				'time': [
					'00:00', '01:00', '02:00',
					'03:00', '04:00', '05:00',
					'06:00', '07:00', '08:00',
					'09:00', '10:00', '11:00',
					'12:00', '13:00', '14:00',
					'15:00', '16:00', '17:00',
					'18:00', '19:00', '20:00',
					'21:00', '22:00', '23:00',
				],
			},outfile)
def load_rain_yearly(yy):
    c = cdsapi.Client()
    outfld = os.path.join(r'E:\ECMWF\Singel_Level\total_rain_water',str(yy))
    mkdir(outfld)
    outfile = os.path.join(outfld, str(yy)+'0000.grib')
    if os.path.exists(outfile):
        return
    print(outfile)
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'variable': 'total_column_rain_water',
            'year': yy,
            'area': [21, 97, 5, 106], # North, West, South, East. Default: global
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'format': 'grib',
        },outfile)
    line_notify(str(yy))
#load_specific_humidity(2013)
#load_specific_humidity_monthly(2011)
#load_temperature_monthly(2014)
#load_Geopotential(2020)
#load_temperature(2020)
#load_2mterperature(2020)
#load_2mterperature_monthly(2013)
#load_Orography(2019)
#load_rain(2021)
#load_total_pwv(2013)
#load_total_pwv_yearly(1981)
#load_surfacepressure(2020)
load_rain_yearly(2010)