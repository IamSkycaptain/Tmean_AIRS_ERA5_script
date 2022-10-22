import os, sys
import datetime
from datetime import timedelta, date
from numpy import loadtxt,zeros,array
import numpy as np
from math import pi, sqrt, sin, tan, cos, atan2, floor,exp,pow,fabs, fmod
import pyproj
import fnmatch
from jdcal import gcal2jd, jd2gcal
def getNumofYear(year):
    doy = datetime.date(year, 12,31).timetuple().tm_yday
    return doy
def mjdToDoy( dmjd):
	# Output:	Day of Year and second of day   ext 234.32324

    # Julian day
	jd = np.floor (dmjd) + 2400000.5;
    # Integer Julian day
	jdi = np.floor (jd);
    # Fractional part of day
	jdf = jd - jdi + 0.5;
    # Really the next calendar day?
	if (jdf >= 1.0):
		jdf = jdf - 1.0;
		jdi  = jdi + 1;

    # hour = jdf * 24.0;
	i = jdi + 68569;
	n = np.floor (4 * i / 146097);

	i = np.floor (i) - np.floor ((146097 * n + 3) / 4);
	year = np.floor (4000 * (i + 1) / 1461001);
	i = i - (np.floor (1461 * year / 4)) + 31;
	month = np.floor (80 * i / 2447);

	dayOfMonth = i - np.floor (2447 * month / 80);

	i = np.floor (month / 11);

	month = np.floor (month + 2 - 12 * i);
	year = np.floor (100 * (n - 49) + year + i);
	secOfday = dmjd - np.floor(dmjd);
	return date2doy(str(year) + '-' + str(month) + '-' + str(dayOfMonth) + ' 00:00:00') + secOfday

def daterange(start_date='%Y-%m-%d', end_date='%Y-%m-%d'):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + timedelta(n)

def isdaytime(localDOY):
	# input
	#	doy : localDOY
	# output
	# 	true or false
	dec = localDOY-floor(localDOY)
	chk = False
	if dec >=0.25 or dec <=0.75 :
		chk = True
		return chk

def getTm(temp,localDOY):
	#input :
	#   localDOY : float  Date/Time in Thailand
	#   temp : degree Celsius
	#output
	#   Tm : float
	#
	Ts = temp+273.15 # unit Kevin
	Tm = 297.0912 # defualt
	if isdaytime(localDOY):
		Tm =0.6066*Ts+113.2914 # day time
	else:
		Tm =0.7938*Ts+57.4856 # night time
	return Tm

def getZHD(pres,lat,alt):
	# inputgg
	#		press : hPa
	#		alt : meter
	#		lat : radian
	#output
	#		ztd : mm
	#
	zhd = 2.2768*pres/(1.0-(0.00266*cos(2.0*lat))-(alt*0.28*10e-6))
	return zhd


def calPWV(utcDOY,ztd,temp,pres,lat,alt):
	#    input
	#         doy -> fload
	#         ztd ->float (mm)
	#         tempearture -> float(cecious)
	#			pressure
	#         lat -> Latitude (radian unit)
	#         alt -> height above MSL (m)
	#   output
	#        pwv -> unit mm.
	zhd = getZHD(pres,lat,alt) # mm
	Tm = getTm(temp,utcDOY+7.0/24.0)
	cof=ii(Tm)
	pwv = (ztd-zhd)*cof; # unit mm
	return zhd, pwv

def getPWV(utcDOY,ztd,listm,lat,alt):
	#    input
	#         doy -> fload
	#         ztd ->float (mm)
	#         meteorology -> list of meteorology (temp,pres...etc)
	#         lat -> Latitude (radian unit)
	#         alt -> height above MSL (m)
	if len(listm)==0:
		temp = 0.0;
		press = 0.0;
		rain1hr	=0.0;
		zhd = 0.0;
		Tm = 0.0;
		cof = 0.0;
		pwv = 0.0;
	else:
		temp, pres, rain1hr = getmdata(utcDOY,listm)
		zhd = getZHD(pres,lat,alt) # mm
		Tm = getTm(temp,utcDOY+7.0/24.0)
		cof=ii(Tm)
		pwv = 1000.0*(ztd-zhd)*cof; # unit mm
	return zhd, pwv, pres, temp ,rain1hr

def getFldList(workdir):
	paths = []
	for child in os.listdir(workdir):
		path = os.path.join(workdir, child)
		if os.path.isdir(path):
			paths.append(path)
	return paths

def ii(Tm):
	# Tm = unit K
	#
    #II Summary of this function goes here
    #Detailed explanation goes here
    #density = unik kg/m^3
    k1 = 77.6900;   # K/hPa
    k2 = 71.2952;   # K/hPa
    k3 = 375463;    # K^2/hPa
    Mw = 18.01528;# g/mol or 18.016
    Md = 28.9644;   # g/mol
    Rv = 461.45;     # J/(kg?K) % J
    kk = k1*Mw/Md;
    density = 1000.0; #waterdensity(t0); % kg/m^3
    cof = 100.0*(1000000.0)/(density*Rv*(k2-kk+k3/Tm));
    return cof

def readmfile(mfile1):
	listm = []
	with open(mfile1) as f:
		for line in f:
			results = [float(i) for i in line.split()]
			listm.append(results)
	f.close()
	return listm
	'''
	if os.path.isfile(mfile2):
		with open(mfile2) as f:
			for line in f:
				results = [float(i) for i in line.split()]
				listm.append(results)
		f.close()
	'''

def getmdata(cdoy,listm):
	# input doy : float
	#          listm : array loadtxt
	temp = 0.0
	press = 0.0
	rain1hr = 0.0

	for i in range(len(listm)):
		doy = listm[i][0]
		temp = listm[i][1]
		press = listm[i][3]
		rain1hr   = listm[i][4]
		tempC = temp
		pressC = press
		rain1hrC = rain1hr
		if doy >= cdoy:
			if i==0:
				break

			doyB = doy
			doyA = listm[i-1][0]
			tempB = temp
			tempA = listm[i-1][1]
			tempC = tempA

			pressB = press
			pressA = listm[i-1][3]
			pressC = pressA

			#print(pressA, ' ' , pressB)
			rain1hrB = rain1hr
			rain1hrA = listm[i-1][4]
			rain1hrC = rain1hrA

			if doyB-doyA > 0:
				tempC =tempC +(tempB-tempA)*(cdoy-doyA)/(doyB-doyA)
				pressC =pressC +(pressB-pressA)*(cdoy-doyA)/(doyB-doyA)
				rain1hrC =rain1hrC +(rain1hrB-rain1hrA)*(cdoy-doyA)/(doyB-doyA)

			break
	return (tempC, pressC, rain1hrC)

def doy2date(year,doy):
	# input  year : integer, doy : float
	# output yyyy-mm-dd H:M:S : string
	datetimeformat = '%Y-%m-%d %H:%M:%S'
	idoy = floor(doy)

	fhh = (doy-idoy)*24.0
	ihh = floor(fhh)
	fmm = (fhh-ihh)*60.0
	imm = floor(fmm)
	fss = (fmm-imm)*60.0
	iss = floor(fss)

	nd = datetime.datetime.strptime(str(year) + ' ' + str(idoy) + ' ' + str(ihh) + ' ' + str(imm) + ' ' + str(iss), '%Y %j %H %M %S')
	return nd.strftime(datetimeformat)

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


def wsow2date(gpsweek,gpsseconds,leapseconds):
	#'GPS time began on January 6, 1980'
	datetimeformat = '%Y-%m-%d %H:%M:%S'
	leap_options = {
		1998 : 12.0,
		1999 : 13.0,
		2000 : 13.0,
		2001 : 13.0,
		2002 : 13.0,
		2003 : 13.0,
		2004 : 13.0,
		2005 : 14.0,
		2006 : 14.0,
		2007 : 14.0,
		2008 : 15.0,
		2009 : 15.0,
		2010 : 15.0,
		2011 : 15.0,
		2012 : 16.0,
		2013 : 16.0,
		2014 : 16.0,
		2015 : 16.0,
		2016 : 16.0
	}
	epoch = datetime.datetime.strptime('1980-01-06 00:00:00',datetimeformat)
	elapsed = datetime.timedelta(days=(gpsweek*7.0),seconds=(gpsseconds+leapseconds))
	return datetime.datetime.strftime(epoch + elapsed,datetimeformat)

def get_date(date ='%Y-%m-%d %H:%M:%S', addtype='days', adds=0):
	datetimeformat = '%Y-%m-%d %H:%M:%S'
	timeNow = datetime.datetime.strptime(date,datetimeformat)
	options ={
		'days' : timeNow + datetime.timedelta(days=adds),
	    'seconds' : timeNow + datetime.timedelta(seconds=adds),
		'microseconds' : timeNow + datetime.timedelta(microseconds=adds),
		'milliseconds' : timeNow + datetime.timedelta(milliseconds=adds),
		'minutes' : timeNow + datetime.timedelta(minutes=adds),
		'hours' : timeNow + datetime.timedelta(hours=adds)
	}
	if adds!=0:
		anotherTime = options[addtype]
	else:
		anotherTime = timeNow
	return anotherTime.strftime(datetimeformat)

def lowerCaseFilename(workdir):
	for root, dirs, files in os.walk(workdir, topdown=True):
		for name in files:
			srcfile = os.path.join(root, name)
			desfile = os.path.join(root, name.lower())
			os.rename(srcfile, desfile)

def ecef2lla(x, y, z, ellps):
	# Output : radian unit
	# 'GRS80'
	ecef = pyproj.Proj(proj='geocent', ellps=ellps)
	lla = pyproj.Proj(proj='latlong', ellps=ellps)
	lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=True)
	return (lat, lon, alt)


def transformEcefToEnu2(originLonLatAlt, ecef):
    """
    Transform tuple ECEF x,y,z (meters) to tuple E,N,U (meters).

    Based on http://en.wikipedia.org/wiki/Geodetic_system
    """
    x, y, z = ecef
    ox, oy, oz = transformLonLatAltToEcef(originLonLatAlt)
    dx, dy, dz = (x - ox, y - oy, z - oz)
    lonDeg, latDeg, _ = originLonLatAlt
    lon = lonDeg * DEG2RAD
    lat = latDeg * DEG2RAD
    return (-sin(lon) * dx + cos(lon) * dy,
            -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz,
            cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz)
#Referance : http://pydoc.net/Python/geocamUtil/0.1/geocamUtil.geomath/
#                 http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
def transformEcefToEnu(originECEF, originLatLon, ecef):
    """
	Input : originECEF, ecef unit meter
		       originLatLon unit radian
	Modified by Chaiyut
    Transform tuple ECEF x,y,z (meters) to tuple E,N,U (meters).
    Based on http://en.wikipedia.org/wiki/Geodetic_system
    """
    x, y, z = ecef
    ox, oy, oz = originECEF
    dx, dy, dz = (x - ox, y - oy, z - oz)
    lat, lon, _ = originLatLon#ecef2lla(originECEF,'GRS80') # output -> Radian Unit
    return (-sin(lon) * dx + cos(lon) * dy,
            -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz,
            cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz)


#  !!!!     Original Code, might need to modify.
def xyzToENU(originLatLon,xyz):
	#input
	# dx,dy,dz
	# lat, lon unit radius
    dx, dy, dz =xyz
    lat, lon = originLatLon#ecef2lla(originECEF,'GRS80') # output -> Radian Unit
    return (-sin(lon) * dx + cos(lon) * dy,
            -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz,
            cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz)

def ecef2enu(originLonLatAlt, dxyz):
    """
    Transform tuple ECEF dx,dy,dz (meters) to tuple E,N,U (meters).
    Based on http://en.wikipedia.org/wiki/Geodetic_system

	I = originLonLatAlt (lat ,long) unit Radian
	I = dxyz -> dx dy dz in ecef coordinate unit meter
	O = enu local coordinate unit meter
    """
    dx, dy, dz = dxyz
    lat, lon = originLonLatAlt
    return (-sin(lon) * dx + cos(lon) * dy,
            -sin(lat) * cos(lon) * dx - sin(lat) * sin(lon) * dy + cos(lat) * dz,
            cos(lat) * cos(lon) * dx + cos(lat) * sin(lon) * dy + sin(lat) * dz)

def transformEnuToEcef(originLonLatAlt, enu):
    """
    Transform tuple E,N,U (meters) to tuple ECEF x,y,z (meters).

    Based on http://en.wikipedia.org/wiki/Geodetic_system
    """
    e, n, u = enu
    lon, lat, _ = originLonLatAlt
    lon = lonDeg * DEG2RAD
    lat = latDeg * DEG2RAD
    ox, oy, oz = transformLonLatAltToEcef(originLonLatAlt)
    return (ox - sin(lon) * e - cos(lon) * sin(lat) * n + cos(lon) * cos(lat) * u,
            oy + cos(lon) * e - sin(lon) * sin(lat) * n + cos(lat) * sin(lon) * u,
            oz + cos(lat) * n + sin(lat) * u)



def chkTime(hms ='H:M:S',chk = 'M:S'):
	hms = hms.split(':')
	ms   = chk.split(':')
	css = ms[1]
	ss = hms[2]
	result = False
	if css == ss:
		if ms[0]=='00':
			result =True
		elif float(hms[1])%float(ms[0])==00:
			result = True
	return result

def 	vmf1_ht(ah,aw,dmjd,dlat,hell,zd):
	doy = dmjd  - 44239.0 + 1.0 - 28.0;
	bh = 0.0029;
	c0h = 0.062;
	if (dlat<0):
		phh  = pi;
		c11h = 0.007;
		c10h = 0.002;
	else:
		phh  = 0.0;
		c11h = 0.005;
		c10h = 0.001;

	ch = c0h + ((cos(doy/365.250*2.0*pi + phh)+1.0)*c11h/2.0 + c10h)*(1.0-cos(dlat));
	sine   = sin(pi/2.0 - zd);
	beta   = bh/( sine + ch  );
	gamma  = ah/( sine + beta);
	topcon = (1.0 + ah/(1.0 + bh/(1.0 + ch)));
	vmf1h   = topcon/(sine+gamma);

	#% C  height correction for hydrotatic part [Niell, 1996]
	a_ht = 2.53e-5;
	b_ht = 5.49e-3;
	c_ht = 1.14e-3;
	hs_km  = hell/1000.0;
	beta         = b_ht/( sine + c_ht);
	gamma        = a_ht/( sine + beta);
	topcon       = (1.0 + a_ht/(1.0 + b_ht/(1.0 + c_ht)));
	ht_corr_coef = 1.0/sine - topcon/(sine + gamma);
	ht_corr      = ht_corr_coef * hs_km;
	mfh = vmf1h + ht_corr; # hydrostatic mappping function

	bw = 0.00146;
	cw = 0.04391;
	beta   = bw/( sine + cw );
	gamma  = aw/( sine + beta);
	topcon = (1.0 + aw/(1.0 + bw/(1.0 + cw)));
	mfw   = topcon/(sine+gamma); # wet mapping function
	return (mfh,mfw)
def gpsFromUTC(year, month, day, hour, min, sec, leapSecs=18):
	 # http://leapsecond.com/java/gpsclock.htm
	'''converts UTC to: gpsWeek, secsOfWeek, gpsDay, secsOfDay
	 96
	 97      a good reference is:  http://www.oc.nps.navy.mil/~jclynch/timsys.html
	 98
	 99      This is based on the following facts (see reference above):
	100
	101      GPS time is basically measured in (atomic) seconds since
	102      January 6, 1980, 00:00:00.0  (the GPS Epoch)
	103
	104      The GPS week starts on Saturday midnight (Sunday morning), and runs
	105      for 604800 seconds.
	106
	107      Currently, GPS time is 13 seconds ahead of UTC (see above reference).
	108      While GPS SVs transmit this difference and the date when another leap
	109      second takes effect, the use of leap seconds cannot be predicted.  This
	110      routine is precise until the next leap second is introduced and has to be
	111      updated after that.
	112
	113      SOW = Seconds of Week
	114      SOD = Seconds of Day
	115
	116      Note:  Python represents time in integer seconds, fractions are lost!!!
	'''
	gpsEpoch = (1980, 1, 6, 0, 0, 0)  # (year, month, day, hh, mm, ss)
	secsInWeek = 604800
	secsInDay = 86400

	secFract = sec % 1
	epochTuple = gpsEpoch + (-1, -1, 0)
	t0 = time.mktime(epochTuple)
	t = time.mktime((year, month, day, hour, min, sec, -1, -1, 0))
	# Note: time.mktime strictly works in localtime and to yield UTC, it should be
	#       corrected with time.timezone
	#       However, since we use the difference, this correction is unnecessary.
	# Warning:  trouble if daylight savings flag is set to -1 or 1 !!!
	t = t + leapSecs
	tdiff = t - t0
	gpsSOW = (tdiff % secsInWeek)  + secFract
	gpsWeek = int(floor(tdiff/secsInWeek))
	gpsDay = int(floor(gpsSOW/secsInDay))
	gpsSOD = (gpsSOW % secsInDay)
	return (gpsWeek, gpsSOW, gpsDay, gpsSOD)

def UTCFromGps(gpsWeek, SOW, leapSecs=14):
	'''converts gps week and seconds to UTC
	137
	138      see comments of inverse function!
	139
	140      SOW = seconds of week
	141      gpsWeek is the full number (not modulo 1024)
	'''
	secFract = SOW % 1
	epochTuple = gpsEpoch + (-1, -1, 0)
	t0 = time.mktime(epochTuple) - time.timezone  #mktime is localtime, correct for UTC
	tdiff = (gpsWeek * secsInWeek) + SOW - leapSecs
	t = t0 + tdiff
	(year, month, day, hh, mm, ss, dayOfWeek, julianDay, daylightsaving) = time.gmtime(t)
	#use gmtime since localtime does not allow to switch off daylighsavings correction!!!
	return (year, month, day, hh, mm, ss + secFract)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def mkdir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def createDic(data):
	dic = {}
	for d in data:
		key = '{:.6f}'.format(d[0])
		dic[key] = d[1]
	return dic

def getZTD(ref_ztd,doy):
	ztd = 0.0
	for dataztd in ref_ztd:
		doyztd = dataztd[0]
		tropztd = dataztd[1]
		if doy==doyztd:
			ztd = tropztd
			break
	return ztd;
def sphericalHM (dlat,dlon):
	'''
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%   Create by  : Chaiyut Charoenphon
	%   Create Date: 05/01/2016
	%
	%       Input: ellipsoidal Coordination
	%       -------------------------------
	%             dlat: latitude in radians
	%             dlon: longitude in radians
	%
	%       Output: longitude and latitude-related functions
	%       -------------------------------
	%             aP : Pnm(sin?)cos(m?)
	%             bP : Pnm(sin?)sin(m?)
	%   Reference:
	%   - J. Böhm, R. Heinkelmann, H. Schuh, Short Note: A Global Model of Pressure
	% and Temperature for Geodetic Applications, Journal of Geodesy,
	% doi:10.1007/s00190-007-0135-3, 2007.
	%   - Heiskanen and Moritz, Physical Geodesy, 1967
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% parameter t
	'''
	t = sin(dlat);

	#% degree n and order m
	n = 9;
	m = 9;

	#% determine n!  (faktorielle)  moved by 1
	dfac = [x for x in range(1,22,1)]
	dfac[1] = 1;
	P = [
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0],
		[0,0,0,0,0,0,0,0,0,0,0]
	]
	for i in range(1,(2*n + 2),1):
		dfac[i+1] = dfac[i]*i;
	#% determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
	for i in range(0,n+1,1):
		for j in range(0,min(i,m)+1,1):
			ir = floor((i - j)/2);
			sum = 0;
			#%(Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-67)
			for k in range(0,ir+1,1):
				sum = sum + pow(-1,k)*dfac[2*i - 2*k + 1]/dfac[k + 1]/dfac[i - k + 1]/dfac[i - j - 2*k + 1]*pow(t,(i - j - 2*k));
			#print(sum);
			#% Legendre functions moved by 1
			P[i + 1][j + 1] = 1.0/pow(2,i)*sqrt(pow((1 - pow(t,2)),(j)))*sum

	# #% spherical harmonics
	i = 0;
	aP = [x for x in range(0,56,1)]
	bP = [x for x in range(0,56,1)]
	for n in range(0,9+1,1):
		for m in range(0,n+1,1):
			i = i + 1;
			aP[i] = (P[n+1][m+1])*cos(m*dlon);
			bP[i] = (P[n+1][m+1])*sin(m*dlon);
	return aP,bP

def gpt( dmjd,dlat,dlon,dhgt,is_msl=0):
	'''
	% Modified JD = JD ? 2400000.5
	% This subroutine determines Global Pressure and Temperature
	% based on Spherical Harmonics up to degree and order 9
	%
	% input data
	% ----------
	% dmjd: modified julian date
	% dlat: ellipsoidal latitude in radians
	% dlon: longitude in radians
	% dhgt: ellipsoidal height in m
	% is_msl: 0 is ellipsoidal height, 1 is orthometric height (default=0)
	% output data
	% -----------
	% pres: pressure in hPa
	% temp: temperature in Celsius
	% undu: Geoid undulation in m (from a 9x9 EGM based model)
	%
	% Johannes Boehm, 2006 June 12
	% rev 2006 June 16: geoid undulation is accounted for
	% rev 2011 July 21: latitude -> ellipsoidal latitude (J. Boehm)
	%
	% Reference:
	% J. Böhm, R. Heinkelmann, H. Schuh, Short Note: A Global Model of Pressure
	% and Temperature for Geodetic Applications, Journal of Geodesy,
	% doi:10.1007/s00190-007-0135-3, 2007.

	% reference day is 28 January
	% this is taken from Niell (1996) to be consistent
	'''
	doy = dmjd  - 44239.0 + 1 - 28;
	#%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	a_geoid = np.array([
		-5.6195e-001,-6.0794e-002,-2.0125e-001,-6.4180e-002,-3.6997e-002,
		+1.0098e+001,+1.6436e+001,+1.4065e+001,+1.9881e+000,+6.4414e-001,
		-4.7482e+000,-3.2290e+000,+5.0652e-001,+3.8279e-001,-2.6646e-002,
		+1.7224e+000,-2.7970e-001,+6.8177e-001,-9.6658e-002,-1.5113e-002,
		+2.9206e-003,-3.4621e+000,-3.8198e-001,+3.2306e-002,+6.9915e-003,
		-2.3068e-003,-1.3548e-003,+4.7324e-006,+2.3527e+000,+1.2985e+000,
		+2.1232e-001,+2.2571e-002,-3.7855e-003,+2.9449e-005,-1.6265e-004,
		+1.1711e-007,+1.6732e+000,+1.9858e-001,+2.3975e-002,-9.0013e-004,
		-2.2475e-003,-3.3095e-005,-1.2040e-005,+2.2010e-006,-1.0083e-006,
		+8.6297e-001,+5.8231e-001,+2.0545e-002,-7.8110e-003,-1.4085e-004,
		-8.8459e-006,+5.7256e-006,-1.5068e-006,+4.0095e-007,-2.4185e-008
	]);

	b_geoid = np.array([
		+0.0000e+000,+0.0000e+000,-6.5993e-002,+0.0000e+000,+6.5364e-002,
		-5.8320e+000,+0.0000e+000,+1.6961e+000,-1.3557e+000,+1.2694e+000,
		+0.0000e+000,-2.9310e+000,+9.4805e-001,-7.6243e-002,+4.1076e-002,
		+0.0000e+000,-5.1808e-001,-3.4583e-001,-4.3632e-002,+2.2101e-003,
		-1.0663e-002,+0.0000e+000,+1.0927e-001,-2.9463e-001,+1.4371e-003,
		-1.1452e-002,-2.8156e-003,-3.5330e-004,+0.0000e+000,+4.4049e-001,
		+5.5653e-002,-2.0396e-002,-1.7312e-003,+3.5805e-005,+7.2682e-005,
		+2.2535e-006,+0.0000e+000,+1.9502e-002,+2.7919e-002,-8.1812e-003,
		+4.4540e-004,+8.8663e-005,+5.5596e-005,+2.4826e-006,+1.0279e-006,
		+0.0000e+000,+6.0529e-002,-3.5824e-002,-5.1367e-003,+3.0119e-005,
		-2.9911e-005,+1.9844e-005,-1.2349e-006,-7.6756e-009,+5.0100e-008
	]);
	'''
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	% Pressure
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	'''
	ap_mean = np.array([
	+1.0108e+003,+8.4886e+000,+1.4799e+000,-1.3897e+001,+3.7516e-003,
	-1.4936e-001,+1.2232e+001,-7.6615e-001,-6.7699e-002,+8.1002e-003,
	-1.5874e+001,+3.6614e-001,-6.7807e-002,-3.6309e-003,+5.9966e-004,
	+4.8163e+000,-3.7363e-001,-7.2071e-002,+1.9998e-003,-6.2385e-004,
	-3.7916e-004,+4.7609e+000,-3.9534e-001,+8.6667e-003,+1.1569e-002,
	+1.1441e-003,-1.4193e-004,-8.5723e-005,+6.5008e-001,-5.0889e-001,
	-1.5754e-002,-2.8305e-003,+5.7458e-004,+3.2577e-005,-9.6052e-006,
	-2.7974e-006,+1.3530e+000,-2.7271e-001,-3.0276e-004,+3.6286e-003,
	-2.0398e-004,+1.5846e-005,-7.7787e-006,+1.1210e-006,+9.9020e-008,
	+5.5046e-001,-2.7312e-001,+3.2532e-003,-2.4277e-003,+1.1596e-004,
	+2.6421e-007,-1.3263e-006,+2.7322e-007,+1.4058e-007,+4.9414e-009
	])

	bp_mean = np.array([
	+0.0000e+000,+0.0000e+000,-1.2878e+000,+0.0000e+000,+7.0444e-001,
	+3.3222e-001,+0.0000e+000,-2.9636e-001,+7.2248e-003,+7.9655e-003,
	+0.0000e+000,+1.0854e+000,+1.1145e-002,-3.6513e-002,+3.1527e-003,
	+0.0000e+000,-4.8434e-001,+5.2023e-002,-1.3091e-002,+1.8515e-003,
	+1.5422e-004,+0.0000e+000,+6.8298e-001,+2.5261e-003,-9.9703e-004,
	-1.0829e-003,+1.7688e-004,-3.1418e-005,+0.0000e+000,-3.7018e-001,
	+4.3234e-002,+7.2559e-003,+3.1516e-004,+2.0024e-005,-8.0581e-006,
	-2.3653e-006,+0.0000e+000,+1.0298e-001,-1.5086e-002,+5.6186e-003,
	+3.2613e-005,+4.0567e-005,-1.3925e-006,-3.6219e-007,-2.0176e-008,
	+0.0000e+000,-1.8364e-001,+1.8508e-002,+7.5016e-004,-9.6139e-005,
	-3.1995e-006,+1.3868e-007,-1.9486e-007,+3.0165e-010,-6.4376e-010
	]);

	ap_amp = np.array([
	-1.0444e-001,+1.6618e-001,-6.3974e-002,+1.0922e+000,+5.7472e-001,
	-3.0277e-001,-3.5087e+000,+7.1264e-003,-1.4030e-001,+3.7050e-002,
	+4.0208e-001,-3.0431e-001,-1.3292e-001,+4.6746e-003,-1.5902e-004,
	+2.8624e+000,-3.9315e-001,-6.4371e-002,+1.6444e-002,-2.3403e-003,
	+4.2127e-005,+1.9945e+000,-6.0907e-001,-3.5386e-002,-1.0910e-003,
	-1.2799e-004,+4.0970e-005,+2.2131e-005,-5.3292e-001,-2.9765e-001,
	-3.2877e-002,+1.7691e-003,+5.9692e-005,+3.1725e-005,+2.0741e-005,
	-3.7622e-007,+2.6372e+000,-3.1165e-001,+1.6439e-002,+2.1633e-004,
	+1.7485e-004,+2.1587e-005,+6.1064e-006,-1.3755e-008,-7.8748e-008,
	-5.9152e-001,-1.7676e-001,+8.1807e-003,+1.0445e-003,+2.3432e-004,
	+9.3421e-006,+2.8104e-006,-1.5788e-007,-3.0648e-008,+2.6421e-010
	]);

	bp_amp = np.array([
	+0.0000e+000,+0.0000e+000,+9.3340e-001,+0.0000e+000,+8.2346e-001,
	+2.2082e-001,+0.0000e+000,+9.6177e-001,-1.5650e-002,+1.2708e-003,
	+0.0000e+000,-3.9913e-001,+2.8020e-002,+2.8334e-002,+8.5980e-004,
	+0.0000e+000,+3.0545e-001,-2.1691e-002,+6.4067e-004,-3.6528e-005,
	-1.1166e-004,+0.0000e+000,-7.6974e-002,-1.8986e-002,+5.6896e-003,
	-2.4159e-004,-2.3033e-004,-9.6783e-006,+0.0000e+000,-1.0218e-001,
	-1.3916e-002,-4.1025e-003,-5.1340e-005,-7.0114e-005,-3.3152e-007,
	+1.6901e-006,+0.0000e+000,-1.2422e-002,+2.5072e-003,+1.1205e-003,
	-1.3034e-004,-2.3971e-005,-2.6622e-006,+5.7852e-007,+4.5847e-008,
	+0.0000e+000,+4.4777e-002,-3.0421e-003,+2.6062e-005,-7.2421e-005,
	+1.9119e-006,+3.9236e-007,+2.2390e-007,+2.9765e-009,-4.6452e-009
	]);

	'''
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	% Temperature
	%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	'''
	at_mean = np.array([
	+1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000,
	+9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002,
	-4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003,
	+7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003,
	+8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003,
	+2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001,
	+5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006,
	-7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003,
	-2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008,
	-5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004,
	-1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009
	]);

	bt_mean = np.array([
	+0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001,
	+1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002,
	+0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003,
	+0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003,
	+5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003,
	-1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002,
	-9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006,
	-1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003,
	+5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008,
	+0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004,
	+1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009
	]);

	at_amp = np.array([
	-1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002,
	+2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003,
	+1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003,
	-6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004,
	+3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004,
	-9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002,
	-8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006,
	-1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004,
	+2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008,
	+7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006,
	-4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009
	]);

	bt_amp = np.array([
	+0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000,
	-1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002,
	+0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003,
	+0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003,
	-1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003,
	-3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001,
	+6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006,
	-1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004,
	-8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007,
	+0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005,
	+1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009
	]);

	#%+++++++++++++++++++++++++++++++++++++++++++++++++++++
	[aP,bP]= sphericalHM(dlat, dlon);
	aP = np.array(aP)
	bP = np.array(bP)
	#% Geoid height
	n = len(a_geoid)+1;
	undu = aP[1:n].dot(a_geoid) + bP[1:n].dot(b_geoid);

	n = len(ap_mean)+1;
	apm =aP[1:n].dot(ap_mean) + bP[1:n].dot(bp_mean);
	apa = aP[1:n].dot(ap_amp) + bP[1:n].dot(bp_amp);
	pres0  = apm + apa*cos(doy/365.25*2.0*pi);

	n = len(at_mean)+1;
	atm= aP[1:n].dot(at_mean) + bP[1:n].dot(bt_mean);
	ata = aP[1:n].dot(at_amp) + bP[1:n].dot(bt_amp);
	temp0 =  atm + ata*cos(doy/365.25*2.0*pi);
	#% height correction for pressure
	#% orthometric height
	hort = dhgt if is_msl==1 else dhgt - undu
	p = pres0*pow(1.0-0.00002260*hort,5.2250);
	t = temp0 - 0.00650*hort;
	return p,t,undu,hort

def gpt2w_THA (dmjd,dlat,dlon,hell,nstat,it,is_msl:0):
	#function [p,T,dT,Tm,e,ah,aw,la,undu] = gpt2_1w_THA (dmjd,dlat,dlon,hell,nstat,it)

	# (c) Department of Geodesy and Geoinformation, Vienna University of
	# Technology, 2013
	#
	# The copyright in this document is vested in the Department of Geodesy and
	# Geoinformation (GEO), Vienna University of Technology, Austria. This document
	# may only be reproduced in whole or in part, or stored in a retrieval
	# system, or transmitted in any form, or by any means electronic,
	# mechanical, photocopying or otherwise, either with the prior permission
	# of GEO or in accordance with the terms of ESTEC Contract No.
	# 4000107329/12/NL/LvH.
	# ---
	#
	# This subroutine determines pressure, temperature, temperature lapse rate,
	# mean temperature of the water vapor, water vapor pressure, hydrostatic
	# and wet mapping function coefficients ah and aw, water vapour decrease
	# factor and geoid undulation for specific sites near the Earth surface.
	# It is based on a 1 x 1 degree external grid file ('gpt2_1wA.grd') with mean
	# values as well as sine and cosine amplitudes for the annual and
	# semiannual variation of the coefficients.
	#
	# c Reference:
	# J. B?hm, G. M?ller, M. Schindelegger, G. Pain, R. Weber, Development of an
	# improved blind model for slant delays in the troposphere (GPT2w),
	# GPS Solutions, 2014, doi:10.1007/s10291-014-0403-7
	#
	# input parameters:
	#
	# dmjd:  modified Julian date (scalar, only one epoch per call is possible)
	# dlat:  ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
	# dlon:  longitude in radians [-pi:pi] or [0:2pi] (vector)
	# hell:  ellipsoidal height in m (vector)
	# nstat: number of stations in dlat, dlon, and hell
	#        maximum possible: not relevant for Matlab version
	# it:    case 1: no time variation but static quantities
	#        case 0: with time variation (annual and semiannual terms)
	#
	# output parameters:
	#
	# p:    pressure in hPa (vector of length nstat)
	# T:    temperature in degrees Celsius (vector of length nstat)
	# dT:   temperature lapse rate in degrees per km (vector of length nstat)
	# Tm:   mean temperature of the water vapor in degrees Kelvin (vector of length nstat)
	# e:    water vapor pressure in hPa (vector of length nstat)
	# ah:   hydrostatic mapping function coefficient at zero height (VMF1)
	#       (vector of length nstat)
	# aw:   wet mapping function coefficient (VMF1) (vector of length nstat)
	# la:   water vapor decrease factor (vector of length nstat)
	# undu: geoid undulation in m (vector of length nstat)
	# is_msl : 0: ellipsoidal height, 1 : orthometric height (default:0)

	# The hydrostatic mapping function coefficients have to be used with the
	# height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
	# coefficients refer to zero height.
	#
	# Example 1 (Vienna, 2 August 2012, with time variation,grid file 'gpt2_1wA.grd):
	#
	#dmjd = 56141.0;
	#dlat = 10.200*pi/180.0;
	#dlon = 100.400*pi/180.0;
	#hell = 156.0;
	#nstat = 1;
	#it = 0;
	#
	# output:
	# p = 1002.788 hPa
	# T = 22.060 deg Celsius
	# dT = -6.230 deg / km
	# Tm = 281.304 K
	# e = 16.742 hPa
	# ah = 0.0012646
	# aw = 0.0005752
	# la = 2.6530
	# undu = 45.76 m
	#
	# Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):
	#
	# dmjd = 56141.d0
	# dlat(1) = 48.20d0*pi/180.d0
	# dlon(1) = 16.37d0*pi/180.d0
	# hell(1) = 156.d0
	# nstat = 1
	# it = 1
	#
	# output:
	# p = 1003.709 hPa
	# T = 11.79 deg Celsius
	# dT = -5.49 deg / km
	# Tm = 273.22 K
	# e = 10.26 hPa
	# ah = 0.0012396
	# aw = 0.0005753
	# la = 2.6358
	# undu = 45.76 m
	#
	#
	# Klemens Lagler, 2 August 2012
	# Johannes Boehm, 6 August 2012, revision
	# Klemens Lagler, 21 August 2012, epoch change to January 1 2000
	# Johannes Boehm, 23 August 2012, adding possibility to determine constant field
	# Johannes Boehm, 27 December 2012, reference added
	# Johannes Boehm, 10 January 2013, correction for dlat = -90 degrees
	#                                  (problem found by Changyong He)
	# Johannes Boehm, 21 May 2013, bug with dmjd removed (input parameter dmjd was replaced
	#                 unintentionally; problem found by Dennis Ferguson)
	# Gregory Pain,   17 June 2013, adding water vapor decrease factor la
	# Gregory Pain,   21 June 2013, using the 1 degree grid : better for calculating zenith wet delays (la)
	# Gregory Pain,   01 July 2013, adding mean temperature of the water vapor Tm
	# Gregory Pain,   30 July 2013, changing the method to calculate the water vapor partial pressure (e)
	# Gregory Pain,   31 July 2013, correction for (dlat = -90 degrees, dlon = 360 degrees)
	# Johannes Boehm, 27 December 2013, copyright notice added
	# Johannes Boehm, 25 August 2014, default input file changed to
	#                 gpt2_1wA.grd (slightly different humidity values)
	# Johannes Boehm, 25 August 2014, reference changed to Boehm et al. in GPS
	#                 Solutions
	# ---
	#  lat    lon   p:a0    A1   B1   A2   B2  T:a0    A1   B1   A2   B2  Q:a0    A1    B1    A2    B2 dT:a0    A1   B1   A2   B2    undu       Hs   h:a0      A1      B1      A2      B2    w:a0      A1      B1      A2      B2  lam:a0      A1      B1      A2      B2    Tm:a0    A1   B1   A2   B2
	data=[
	[21.5,	96.5,	92024,	400,	-7,	-13,	-45,	296.3,	-1.5,	1.5,	-2,	-0.3,	13.28,	-4.34,	-3.01,	0.23,	-0.7,	-5.8,	-0.6,	-0.9,	0.3,	-0.1,	-45.35,	801.33,	1.2817,	-0.0053,	-0.0012,	-0.0008,	-0.0006,	0.5674,	-0.0752,	-0.0342,	-0.0041,	0.0021,	3.0122,	1.0544,	0.2471,	0.3619,	-0.0198,	284.7,	0.3,	1,	-0.8,	-0.4],
	[21.5,	97.5,	88076,	354,	-14,	-29,	-54,	294.2,	-1.7,	1.4,	-2,	-0.2,	12.55,	-4.14,	-2.93,	0.41,	-0.66,	-6.4,	-1.6,	-1,	0.2,	-0.1,	-42.37,	1181.94,	1.2814,	-0.0058,	-0.0017,	-0.0006,	-0.0005,	0.5496,	-0.0735,	-0.0302,	-0.0049,	0.0012,	3.1104,	1.0644,	0.1526,	0.4327,	-0.0338,	282.9,	0.2,	0.7,	-0.7,	-0.3],
	[21.5,	98.5,	89055,	373,	-32,	-22,	-62,	294.6,	-1.6,	1.4,	-1.9,	0,	12.6,	-4.07,	-2.81,	0.55,	-0.55,	-5.9,	-1.7,	-0.8,	0.1,	-0.1,	-40.37,	1086.33,	1.2815,	-0.0058,	-0.0014,	-0.0008,	-0.0003,	0.556,	-0.0722,	-0.0269,	-0.0063,	-0.0003,	2.8987,	0.89,	0.0937,	0.3752,	-0.0106,	283.3,	0.2,	0.7,	-0.6,	-0.2],
	[21.5,	99.5,	85734,	340,	-32,	-43,	-67,	292.7,	-1.9,	1.2,	-1.8,	0.2,	11.99,	-3.75,	-2.14,	0.33,	-0.47,	-6,	-1.6,	-0.8,	0.2,	-0.1,	-38.26,	1416.63,	1.2811,	-0.0063,	-0.0019,	-0.0004,	-0.0003,	0.5383,	-0.0715,	-0.0267,	-0.0049,	-0.0019,	3.1182,	1.0291,	0.1683,	0.3873,	-0.004,	281.7,	-0.1,	0.6,	-0.6,	0],
	[21.5,	100.5,	89119,	394,	-47,	-34,	-77,	293.8,	-2,	0.8,	-1.5,	0.2,	12.72,	-3.75,	-1.55,	0.1,	-0.28,	-4.7,	-1.1,	0.1,	-0.1,	0.1,	-36.51,	1083.01,	1.2811,	-0.006,	-0.0017,	-0.0006,	-0.0002,	0.5559,	-0.0692,	-0.0259,	-0.0045,	-0.0031,	2.836,	0.7846,	0.2264,	0.2133,	0.0494,	283.1,	-0.1,	0.8,	-0.7,	0.1],
	[21.5,	101.5,	89464,	413,	-45,	-38,	-81,	293.5,	-2.5,	0.3,	-1.4,	0.2,	12.97,	-3.77,	-1.15,	-0.09,	-0.17,	-4.8,	-1.1,	0.4,	0.1,	0.2,	-34.78,	1051.32,	1.2809,	-0.0062,	-0.0019,	-0.0005,	-0.0002,	0.5568,	-0.0677,	-0.0251,	-0.0035,	-0.0039,	2.8511,	0.7368,	0.2313,	0.1612,	0.051,	283.1,	-0.4,	0.7,	-0.7,	0.2],
	[21.5,	102.5,	89645,	433,	-40,	-40,	-84,	293.1,	-3,	-0.1,	-1.4,	0.1,	13.05,	-3.91,	-0.97,	-0.24,	-0.06,	-4.9,	-0.6,	0.4,	0.3,	0.3,	-32.9,	1036.26,	1.2806,	-0.0063,	-0.0021,	-0.0004,	-0.0002,	0.5573,	-0.0658,	-0.0233,	-0.0028,	-0.0045,	2.8503,	0.66,	0.1983,	0.1124,	0.0446,	283,	-0.7,	0.5,	-0.8,	0.2],
	[21.5,	103.5,	90829,	483,	-35,	-42,	-86,	293.2,	-4,	-0.4,	-1.3,	0,	13.22,	-4.18,	-0.95,	-0.34,	0.01,	-4.2,	-0.8,	0.9,	-0.3,	0.6,	-31.16,	925.87,	1.2803,	-0.0066,	-0.0021,	-0.0004,	-0.0002,	0.5634,	-0.0629,	-0.0209,	-0.0025,	-0.005,	2.7757,	0.5392,	0.1497,	0.0672,	0.0319,	283.1,	-1.2,	0.3,	-0.8,	0.2],
	[21.5,	104.5,	94263,	600,	-26,	-41,	-86,	294.5,	-5.1,	-0.9,	-1.3,	-0.3,	13.54,	-4.52,	-0.93,	-0.35,	0.1,	-5.8,	-0.1,	0.4,	0,	0.4,	-29.74,	608.36,	1.2799,	-0.007,	-0.0021,	-0.0004,	-0.0003,	0.5817,	-0.0589,	-0.0183,	-0.0026,	-0.0056,	2.5629,	0.3807,	0.112,	0.038,	0.037,	284.1,	-1.9,	0.1,	-0.8,	0.2],
	[21.5,	105.5,	99391,	768,	12,	-24,	-77,	296.7,	-5.8,	-1.3,	-1.2,	-0.6,	14.31,	-5.12,	-1.07,	-0.37,	0.19,	-6.2,	-0.1,	0.1,	-0.2,	0.1,	-28.45,	148.95,	1.2795,	-0.0072,	-0.0022,	-0.0006,	-0.0005,	0.6086,	-0.0543,	-0.0164,	-0.0027,	-0.0056,	2.3658,	0.2184,	0.0559,	0.0097,	0.0292,	285.8,	-2.6,	-0.1,	-0.8,	0.1],
	[21.5,	106.5,	98030,	745,	28,	-28,	-77,	295.4,	-5.9,	-1.7,	-1.1,	-0.7,	13.96,	-5.24,	-1.09,	-0.36,	0.22,	-5.7,	0,	0.3,	0.1,	0.3,	-26.23,	269.42,	1.2796,	-0.0073,	-0.0024,	-0.0006,	-0.0005,	0.6013,	-0.0529,	-0.0159,	-0.0029,	-0.0062,	2.4573,	0.1899,	0.0371,	0.0053,	0.0279,	285.1,	-2.6,	-0.2,	-0.7,	0.1],
	[20.5,	96.5,	90470,	359,	0,	-16,	-41,	295.5,	-1.3,	1.5,	-1.9,	-0.3,	13.05,	-4.33,	-3.17,	0.32,	-0.84,	-6,	-0.5,	-0.7,	0.1,	-0.1,	-45.2,	950.84,	1.2822,	-0.0044,	-0.0009,	-0.0007,	-0.0006,	0.5609,	-0.0747,	-0.0354,	-0.0027,	0.0016,	3.1985,	1.1815,	0.2994,	0.3984,	0.0161,	284.1,	0.7,	1.1,	-0.8,	-0.4],
	[20.5,	97.5,	89160,	351,	-14,	-18,	-49,	295.1,	-1.2,	1.7,	-1.9,	-0.1,	12.85,	-4.28,	-3.23,	0.59,	-0.75,	-6.5,	-1.5,	-0.9,	0.1,	-0.1,	-42.11,	1076.81,	1.2822,	-0.0047,	-0.0008,	-0.0007,	-0.0004,	0.5562,	-0.0737,	-0.0304,	-0.0051,	0.0005,	3.0484,	0.9883,	0.1091,	0.4323,	-0.0131,	283.6,	0.6,	1,	-0.7,	-0.3],
	[20.5,	98.5,	90350,	372,	-32,	-10,	-58,	295.7,	-1.3,	1.7,	-1.9,	0.1,	13.14,	-4.03,	-2.89,	0.6,	-0.65,	-4.4,	-0.3,	0.3,	-0.2,	0.1,	-40.04,	961.59,	1.2822,	-0.0046,	-0.0006,	-0.0009,	-0.0003,	0.5632,	-0.0726,	-0.0272,	-0.0061,	-0.0011,	2.8931,	0.8865,	0.0894,	0.3737,	0.0107,	284.2,	0.6,	1,	-0.7,	-0.1],
	[20.5,	99.5,	90348,	385,	-35,	-15,	-67,	295.3,	-1.3,	1.3,	-1.6,	0.2,	13.3,	-3.79,	-1.85,	0.19,	-0.43,	-4.9,	-1.1,	0.1,	-0.3,	0.1,	-37.89,	962.84,	1.282,	-0.0048,	-0.0009,	-0.0007,	-0.0002,	0.5622,	-0.0715,	-0.0275,	-0.0044,	-0.0029,	2.88,	0.8449,	0.2241,	0.2528,	0.0578,	284.1,	0.5,	1.1,	-0.7,	0],
	[20.5,	100.5,	92013,	413,	-35,	-16,	-74,	295.5,	-1.5,	0.8,	-1.4,	0.2,	13.77,	-3.68,	-1.12,	-0.14,	-0.21,	-5,	-0.8,	0.2,	-0.2,	0.1,	-36.01,	806.06,	1.2817,	-0.0047,	-0.0011,	-0.0006,	-0.0002,	0.5704,	-0.0694,	-0.0269,	-0.0031,	-0.0041,	2.77,	0.7394,	0.2674,	0.1509,	0.0765,	284.7,	0.3,	1.1,	-0.8,	0.2],
	[20.5,	101.5,	90648,	406,	-27,	-24,	-75,	294.3,	-2,	0.5,	-1.3,	0.2,	13.6,	-3.7,	-0.83,	-0.23,	-0.08,	-4.5,	-0.9,	0.4,	-0.4,	0.2,	-34.28,	937.84,	1.2815,	-0.005,	-0.0015,	-0.0005,	-0.0002,	0.5625,	-0.0682,	-0.0254,	-0.0024,	-0.0046,	2.8687,	0.7338,	0.2671,	0.1347,	0.0715,	283.9,	0.1,	0.9,	-0.8,	0.2],
	[20.5,	102.5,	90406,	416,	-23,	-31,	-78,	293.8,	-2.4,	0.2,	-1.3,	0.2,	13.4,	-3.81,	-0.74,	-0.3,	0.01,	-4.1,	-0.6,	0.7,	-0.5,	0.3,	-32.45,	964.1,	1.2813,	-0.0052,	-0.0016,	-0.0004,	-0.0003,	0.561,	-0.0657,	-0.0229,	-0.0022,	-0.0048,	2.839,	0.6411,	0.2169,	0.102,	0.0512,	283.6,	-0.2,	0.7,	-0.8,	0.2],
	[20.5,	103.5,	89625,	434,	-21,	-41,	-81,	292.9,	-3.5,	0,	-1.3,	0.1,	13.12,	-3.93,	-0.77,	-0.36,	0.05,	-4.8,	-0.1,	0.5,	0.3,	0.4,	-30.7,	1042.08,	1.2811,	-0.0057,	-0.0017,	-0.0003,	-0.0003,	0.5568,	-0.0638,	-0.0206,	-0.0024,	-0.0051,	2.8775,	0.5945,	0.1719,	0.0885,	0.029,	283,	-0.6,	0.6,	-0.7,	0.2],
	[20.5,	104.5,	93540,	559,	-15,	-37,	-84,	294.4,	-4.8,	-0.5,	-1.1,	-0.3,	13.55,	-4.15,	-0.81,	-0.46,	0.1,	-6.2,	0.3,	0.2,	0.1,	0.2,	-29.12,	675.78,	1.2807,	-0.0062,	-0.0016,	-0.0003,	-0.0003,	0.578,	-0.0595,	-0.0176,	-0.0031,	-0.0051,	2.6167,	0.434,	0.109,	0.0462,	0.024,	284.2,	-1.4,	0.4,	-0.7,	0.2],
	[20.5,	105.5,	98607,	730,	19,	-27,	-77,	296.5,	-5.7,	-1.2,	-0.9,	-0.6,	14.46,	-4.65,	-0.92,	-0.54,	0.14,	-6.1,	0.2,	0.5,	-0.2,	0.2,	-27.16,	217.44,	1.2804,	-0.0066,	-0.0017,	-0.0004,	-0.0005,	0.6045,	-0.0552,	-0.0158,	-0.0035,	-0.0051,	2.4526,	0.282,	0.0699,	0.0148,	0.0218,	285.8,	-2.2,	0.1,	-0.7,	0],
	[20.5,	106.5,	100517,	797,	58,	-19,	-70,	297.1,	-5.7,	-2,	-0.6,	-0.8,	15.2,	-5.1,	-1.03,	-0.44,	0.06,	-6.2,	-0.4,	0.9,	-0.4,	0.2,	-24.72,	48.55,	1.2802,	-0.0066,	-0.0019,	-0.0004,	-0.0006,	0.6134,	-0.0528,	-0.0148,	-0.004,	-0.005,	2.5273,	0.2205,	0.0496,	0.0305,	0.003,	286.4,	-2.3,	-0.2,	-0.6,	0],
	[19.5,	96.5,	93985,	365,	-5,	6,	-30,	297.4,	-0.6,	1.7,	-1.7,	-0.3,	14.01,	-4.46,	-3.29,	0.35,	-0.99,	-5.9,	-0.5,	-0.6,	0,	-0.2,	-45.43,	619.89,	1.2827,	-0.003,	0.0001,	-0.001,	-0.0005,	0.5781,	-0.0756,	-0.0367,	-0.0023,	0.0003,	3.0356,	1.0005,	0.3268,	0.2593,	0.0336,	285.8,	1.2,	1.5,	-0.9,	-0.4],
	[19.5,	97.5,	92903,	360,	-20,	5,	-41,	297.4,	-0.7,	1.8,	-1.8,	-0.1,	13.88,	-4.38,	-3.42,	0.71,	-0.86,	-5.9,	0.6,	-0.2,	0.3,	-0.1,	-42.16,	721.04,	1.2828,	-0.0032,	0.0002,	-0.001,	-0.0004,	0.5749,	-0.0741,	-0.0313,	-0.0047,	-0.0005,	2.9015,	0.8702,	0.1294,	0.3547,	0.0182,	285.4,	1.1,	1.4,	-0.7,	-0.3],
	[19.5,	98.5,	91040,	364,	-26,	-2,	-52,	296.1,	-1,	1.8,	-1.8,	0.1,	13.63,	-4.02,	-2.8,	0.58,	-0.69,	-5.9,	-1.4,	-0.2,	-0.1,	-0.1,	-39.92,	896.24,	1.2827,	-0.0037,	-0.0001,	-0.0009,	-0.0003,	0.5663,	-0.0735,	-0.0284,	-0.0049,	-0.0016,	2.9506,	0.8879,	0.1172,	0.3639,	0.0236,	284.6,	1,	1.3,	-0.7,	-0.1],
	[19.5,	99.5,	93323,	398,	-33,	2,	-62,	297,	-1,	1.4,	-1.5,	0.2,	14.31,	-3.7,	-1.62,	0.06,	-0.38,	-5.1,	-0.5,	0.3,	-0.3,	0.1,	-37.87,	681.21,	1.2825,	-0.0037,	-0.0002,	-0.0008,	-0.0002,	0.5773,	-0.0715,	-0.0279,	-0.0035,	-0.0034,	2.8019,	0.7788,	0.253,	0.1919,	0.0824,	285.6,	0.9,	1.4,	-0.8,	0.1],
	[19.5,	100.5,	93111,	403,	-23,	-4,	-65,	296.3,	-1.2,	0.9,	-1.3,	0.1,	14.38,	-3.54,	-0.84,	-0.3,	-0.16,	-4.6,	-0.9,	0.3,	-0.5,	0.1,	-35.72,	702.61,	1.2822,	-0.0038,	-0.0007,	-0.0006,	-0.0002,	0.5755,	-0.0701,	-0.0272,	-0.0022,	-0.0041,	2.8053,	0.7428,	0.3057,	0.1212,	0.0861,	285.4,	0.7,	1.3,	-0.8,	0.1],
	[19.5,	101.5,	91478,	393,	-10,	-15,	-65,	295.1,	-1.3,	0.6,	-1.2,	0.1,	13.96,	-3.54,	-0.63,	-0.33,	-0.02,	-5.5,	-1.7,	0.2,	-0.3,	0.2,	-33.79,	859.61,	1.282,	-0.0041,	-0.0011,	-0.0005,	-0.0003,	0.5666,	-0.0687,	-0.0251,	-0.0018,	-0.0043,	2.8608,	0.7191,	0.289,	0.1157,	0.0736,	284.5,	0.5,	1.1,	-0.8,	0.1],
	[19.5,	102.5,	89194,	375,	-2,	-29,	-66,	293.5,	-1.9,	0.4,	-1.3,	0,	13.39,	-3.58,	-0.52,	-0.39,	0.06,	-5.4,	-1.2,	0.2,	0.1,	0.2,	-31.83,	1081.76,	1.2818,	-0.0045,	-0.0014,	-0.0003,	-0.0004,	0.5544,	-0.0669,	-0.0229,	-0.0019,	-0.0046,	2.9833,	0.7083,	0.2778,	0.1073,	0.0576,	283.3,	0.3,	0.9,	-0.7,	0.1],
	[19.5,	103.5,	88469,	388,	-1,	-37,	-71,	292.5,	-2.9,	0.3,	-1.3,	0,	13.06,	-3.62,	-0.6,	-0.39,	0.1,	-5,	-0.2,	0.4,	0.2,	0.4,	-29.78,	1154.42,	1.2817,	-0.0048,	-0.0014,	-0.0002,	-0.0004,	0.5506,	-0.0645,	-0.0207,	-0.0027,	-0.0049,	3.0197,	0.6868,	0.2109,	0.1217,	0.0322,	282.7,	-0.1,	0.7,	-0.7,	0.2],
	[19.5,	104.5,	93712,	532,	2,	-25,	-78,	294.6,	-4.5,	-0.1,	-1.1,	-0.2,	13.82,	-3.81,	-0.76,	-0.54,	0.09,	-5.9,	0.4,	0.2,	0,	0.3,	-28.05,	659.43,	1.2813,	-0.0054,	-0.0012,	-0.0003,	-0.0004,	0.5784,	-0.0596,	-0.0177,	-0.004,	-0.0046,	2.6791,	0.4871,	0.1141,	0.0633,	0.0165,	284.5,	-1,	0.6,	-0.7,	0.1],
	[19.5,	105.5,	99832,	727,	35,	-14,	-72,	297.4,	-5.4,	-1.1,	-0.6,	-0.6,	15.01,	-4.17,	-0.89,	-0.63,	0.02,	-6.4,	-0.3,	0.8,	-0.4,	0.3,	-25.6,	109.01,	1.2809,	-0.0061,	-0.0012,	-0.0003,	-0.0005,	0.6101,	-0.0549,	-0.0159,	-0.0045,	-0.0044,	2.4741,	0.3244,	0.0839,	0.0321,	0.0098,	286.5,	-1.9,	0.3,	-0.6,	0],
	[19.5,	106.5,	101035,	773,	69,	-12,	-66,	297.9,	-5.3,	-2.1,	-0.2,	-0.7,	15.95,	-4.51,	-0.88,	-0.48,	-0.08,	-6.4,	-0.4,	2.4,	-1,	0.2,	-23.07,	1.38,	1.2809,	-0.0061,	-0.0014,	-0.0003,	-0.0005,	0.6149,	-0.0533,	-0.0154,	-0.0052,	-0.0046,	2.6379,	0.2881,	0.1115,	0.059,	0.0078,	287,	-2,	0,	-0.5,	-0.1],
	[18.5,	96.5,	96843,	358,	-9,	21,	-19,	298.6,	-0.2,	1.6,	-1.5,	-0.3,	14.94,	-4.35,	-3.05,	0.16,	-1.09,	-5.8,	-1.2,	-0.3,	-0.5,	0,	-45.59,	357.98,	1.2831,	-0.0019,	0.0008,	-0.001,	-0.0004,	0.5922,	-0.0756,	-0.038,	-0.0011,	-0.0007,	2.969,	0.9243,	0.3851,	0.1486,	0.0559,	287,	1.6,	1.7,	-0.9,	-0.4],
	[18.5,	97.5,	94165,	343,	-15,	14,	-32,	297.6,	-0.3,	1.8,	-1.6,	-0.1,	14.49,	-4.3,	-3.09,	0.48,	-0.96,	-5.3,	0.1,	-0.2,	-0.1,	0,	-42.33,	604.14,	1.283,	-0.0023,	0.0006,	-0.001,	-0.0003,	0.5807,	-0.0744,	-0.0333,	-0.0027,	-0.0009,	2.9404,	0.8571,	0.2005,	0.2771,	0.0177,	285.9,	1.5,	1.6,	-0.8,	-0.3],
	[18.5,	98.5,	93698,	362,	-25,	12,	-46,	297.4,	-0.6,	1.7,	-1.6,	0,	14.6,	-3.98,	-2.52,	0.43,	-0.66,	-4.8,	0.2,	0.4,	-0.3,	0.1,	-39.91,	646.86,	1.283,	-0.0027,	0.0004,	-0.0009,	-0.0002,	0.5794,	-0.0731,	-0.0294,	-0.0036,	-0.0021,	2.8967,	0.8067,	0.166,	0.2841,	0.0465,	285.8,	1.3,	1.5,	-0.8,	-0.1],
	[18.5,	99.5,	95392,	394,	-28,	14,	-55,	298.1,	-0.7,	1.4,	-1.4,	0.1,	15.04,	-3.57,	-1.39,	-0.1,	-0.29,	-5.2,	-0.1,	0,	0,	0.1,	-37.83,	490.38,	1.2829,	-0.0029,	0.0002,	-0.0008,	-0.0002,	0.5877,	-0.0712,	-0.0281,	-0.0027,	-0.0034,	2.7741,	0.7349,	0.2657,	0.1536,	0.0861,	286.6,	1.2,	1.6,	-0.8,	0.1],
	[18.5,	100.5,	94894,	398,	-14,	8,	-56,	297.4,	-0.8,	0.9,	-1.2,	0,	14.92,	-3.41,	-0.7,	-0.44,	-0.09,	-4.7,	-0.1,	-0.1,	0.1,	0,	-35.42,	537.64,	1.2826,	-0.0031,	-0.0003,	-0.0006,	-0.0003,	0.5848,	-0.0696,	-0.0268,	-0.0018,	-0.0037,	2.7687,	0.7014,	0.3067,	0.093,	0.0817,	286.2,	1,	1.4,	-0.9,	0.1],
	[18.5,	101.5,	94829,	409,	-3,	3,	-58,	297,	-1.2,	0.6,	-1.1,	-0.1,	14.63,	-3.45,	-0.58,	-0.48,	-0.01,	-5.1,	-0.1,	0,	-0.1,	-0.1,	-33.24,	546.59,	1.2823,	-0.0032,	-0.0006,	-0.0005,	-0.0004,	0.5849,	-0.0666,	-0.0237,	-0.002,	-0.0033,	2.6825,	0.6063,	0.2332,	0.0762,	0.0409,	285.9,	0.8,	1.2,	-0.8,	0.1],
	[18.5,	102.5,	95299,	427,	6,	2,	-61,	297.1,	-1.7,	0.4,	-1.2,	-0.1,	14.32,	-3.63,	-0.63,	-0.48,	0.07,	-5.6,	0.2,	0.1,	-0.1,	0.1,	-31.6,	506.14,	1.2822,	-0.0034,	-0.0008,	-0.0005,	-0.0004,	0.5878,	-0.0631,	-0.0211,	-0.0026,	-0.0029,	2.5761,	0.4998,	0.1619,	0.0621,	0.0144,	285.9,	0.4,	1,	-0.8,	0.1],
	[18.5,	103.5,	94919,	444,	11,	0,	-63,	296.4,	-2.4,	0.5,	-1.4,	-0.2,	14.06,	-3.75,	-0.73,	-0.44,	0.12,	-5.8,	0.1,	-0.1,	0.2,	0.1,	-29.85,	543.47,	1.2821,	-0.0036,	-0.0008,	-0.0005,	-0.0005,	0.5852,	-0.0608,	-0.0196,	-0.0035,	-0.0032,	2.5986,	0.4744,	0.1232,	0.075,	0.0084,	285.5,	0.1,	0.9,	-0.8,	0],
	[18.5,	104.5,	93414,	469,	17,	-9,	-65,	294.7,	-3.6,	0.2,	-1.3,	-0.3,	13.96,	-3.56,	-0.75,	-0.48,	0.08,	-5.4,	0.4,	0.3,	0.1,	0.5,	-27.25,	684.34,	1.2819,	-0.0043,	-0.0009,	-0.0004,	-0.0004,	0.5762,	-0.0598,	-0.019,	-0.0047,	-0.0042,	2.7876,	0.5636,	0.1328,	0.1003,	0.0192,	284.6,	-0.4,	0.8,	-0.7,	0.1],
	[18.5,	105.5,	98221,	617,	34,	-3,	-66,	296.9,	-4.7,	-0.4,	-0.9,	-0.4,	14.89,	-3.44,	-0.85,	-0.7,	-0.03,	-6.4,	0.1,	0.5,	-0.1,	0.1,	-23.81,	250.64,	1.2815,	-0.0051,	-0.0008,	-0.0003,	-0.0004,	0.6009,	-0.0557,	-0.0172,	-0.0056,	-0.0043,	2.5683,	0.4474,	0.1104,	0.047,	0.0207,	286.2,	-1.2,	0.5,	-0.6,	0],
	[18.5,	106.5,	100967,	708,	66,	-1,	-60,	298.5,	-4.8,	-1.4,	-0.3,	-0.6,	15.99,	-3.62,	-0.77,	-0.65,	-0.14,	-6.7,	-0.4,	1.9,	-0.8,	0.1,	-21.14,	6.47,	1.2814,	-0.0054,	-0.001,	-0.0002,	-0.0004,	0.6144,	-0.0526,	-0.0166,	-0.0061,	-0.0045,	2.6076,	0.3678,	0.1449,	0.051,	0.0268,	287.3,	-1.6,	0.3,	-0.5,	-0.1],
	[17.5,	96.5,	99787,	346,	-11,	31,	-7,	299.7,	0,	1.4,	-1.3,	-0.5,	15.93,	-4.02,	-2.38,	-0.25,	-0.99,	-4.4,	1.1,	0.3,	0.2,	0.5,	-45.58,	94.3,	1.2832,	-0.0009,	0.0012,	-0.0009,	-0.0003,	0.6075,	-0.0742,	-0.0384,	-0.0001,	-0.002,	2.9188,	0.8974,	0.491,	0.0327,	0.1255,	288.2,	1.9,	1.9,	-0.9,	-0.3],
	[17.5,	97.5,	96771,	324,	-11,	22,	-20,	298.7,	0.4,	1.6,	-1.3,	-0.1,	15.15,	-4.14,	-2.35,	-0.07,	-0.92,	-5.6,	-0.8,	-0.3,	-0.3,	0,	-42.3,	366.33,	1.2832,	-0.0013,	0.0009,	-0.0009,	-0.0003,	0.594,	-0.0733,	-0.0346,	-0.0011,	-0.0018,	2.8643,	0.7869,	0.3359,	0.1077,	0.0492,	287,	1.8,	1.8,	-0.8,	-0.2],
	[17.5,	98.5,	93918,	335,	-10,	12,	-35,	297.2,	-0.1,	1.6,	-1.3,	0,	14.65,	-3.71,	-1.82,	-0.04,	-0.6,	-5.8,	-1.3,	-0.2,	-0.4,	0.1,	-39.56,	627.99,	1.2832,	-0.0021,	0.0004,	-0.0007,	-0.0002,	0.5806,	-0.0727,	-0.0309,	-0.0022,	-0.0024,	2.9172,	0.808,	0.272,	0.1837,	0.0486,	285.8,	1.6,	1.6,	-0.8,	-0.1],
	[17.5,	99.5,	97494,	382,	-22,	23,	-45,	299.2,	-0.4,	1.3,	-1.2,	0,	15.59,	-3.28,	-1.08,	-0.29,	-0.2,	-5.2,	-0.9,	0,	0,	0,	-37.35,	300.53,	1.2831,	-0.0022,	0.0004,	-0.0007,	-0.0002,	0.5989,	-0.0695,	-0.0279,	-0.0025,	-0.0033,	2.707,	0.6952,	0.2648,	0.1232,	0.0778,	287.4,	1.4,	1.7,	-0.8,	0.1],
	[17.5,	100.5,	95371,	375,	-5,	13,	-46,	297.8,	-0.7,	0.9,	-1.1,	-0.1,	15.1,	-3.22,	-0.48,	-0.58,	0.01,	-5.3,	-0.9,	0.3,	-0.6,	0.2,	-34.86,	494.78,	1.2829,	-0.0026,	-0.0001,	-0.0005,	-0.0004,	0.5878,	-0.0682,	-0.0267,	-0.0017,	-0.0033,	2.7863,	0.6743,	0.3336,	0.0682,	0.0879,	286.4,	1.2,	1.5,	-0.8,	0.1],
	[17.5,	101.5,	95447,	397,	3,	9,	-51,	297.5,	-1.3,	0.6,	-1,	-0.1,	14.85,	-3.26,	-0.37,	-0.63,	0.06,	-5.2,	-0.2,	0.4,	-0.5,	0.2,	-32.45,	490.45,	1.2827,	-0.0029,	-0.0004,	-0.0004,	-0.0004,	0.589,	-0.0644,	-0.023,	-0.0023,	-0.0027,	2.6854,	0.5662,	0.2491,	0.0525,	0.0393,	286.2,	0.8,	1.3,	-0.8,	0],
	[17.5,	102.5,	98747,	453,	5,	19,	-56,	299.3,	-1.9,	0.6,	-1.1,	-0.1,	14.93,	-3.53,	-0.54,	-0.66,	0.08,	-6.3,	-0.5,	0.1,	-0.2,	0.2,	-30.96,	194.34,	1.2827,	-0.0029,	-0.0003,	-0.0005,	-0.0004,	0.607,	-0.0595,	-0.0197,	-0.0031,	-0.0018,	2.4335,	0.403,	0.1392,	0.0286,	0.0026,	287.3,	0.4,	1.1,	-0.8,	0],
	[17.5,	103.5,	98962,	464,	10,	24,	-54,	299.1,	-2.2,	0.7,	-1.3,	-0.2,	14.83,	-3.71,	-0.79,	-0.58,	0.1,	-5.8,	-0.1,	0.2,	-0.3,	0.4,	-29.22,	175.46,	1.2826,	-0.0029,	-0.0003,	-0.0006,	-0.0005,	0.6076,	-0.0577,	-0.0187,	-0.0042,	-0.0021,	2.4383,	0.3802,	0.094,	0.0438,	0.0019,	287.2,	0.2,	1,	-0.8,	0],
	[17.5,	104.5,	97962,	471,	21,	21,	-52,	297.8,	-2.8,	0.6,	-1.3,	-0.3,	14.78,	-3.56,	-0.83,	-0.55,	0.06,	-6.2,	0,	0,	-0.1,	0.4,	-26.22,	267.43,	1.2823,	-0.0031,	-0.0003,	-0.0006,	-0.0005,	0.6012,	-0.0567,	-0.0184,	-0.0055,	-0.0028,	2.527,	0.4325,	0.0967,	0.0709,	0.0117,	286.7,	0,	0.9,	-0.7,	0],
	[17.5,	105.5,	96102,	488,	34,	10,	-54,	296.1,	-3.5,	0.2,	-1.1,	-0.4,	14.85,	-3.17,	-0.8,	-0.6,	-0.06,	-5.5,	0.3,	0.5,	-0.4,	0.4,	-22.25,	436.92,	1.2821,	-0.0038,	-0.0005,	-0.0004,	-0.0004,	0.5904,	-0.0558,	-0.0186,	-0.0067,	-0.004,	2.7269,	0.5366,	0.1341,	0.094,	0.0326,	285.8,	-0.4,	0.8,	-0.6,	0],
	[17.5,	106.5,	99148,	582,	54,	9,	-53,	297.8,	-4.2,	-0.4,	-0.6,	-0.5,	15.46,	-2.87,	-0.73,	-0.74,	-0.18,	-6.3,	0.2,	0.8,	-0.6,	0.2,	-19.16,	166.2,	1.2818,	-0.0045,	-0.0005,	-0.0003,	-0.0004,	0.6059,	-0.052,	-0.0178,	-0.0071,	-0.0044,	2.5923,	0.4676,	0.1497,	0.0538,	0.0351,	286.8,	-1,	0.6,	-0.5,	-0.1],
	[16.5,	96.5,	100783,	330,	3,	26,	2,	300.3,	-0.2,	0.7,	-0.9,	-0.7,	16.94,	-2.85,	-0.98,	-0.73,	-0.52,	-5.2,	2.1,	0.5,	0.4,	1.1,	-45.22,	7.43,	1.2833,	-0.0005,	0.0011,	-0.0007,	-0.0003,	0.6141,	-0.0724,	-0.038,	-0.0001,	-0.0033,	3.139,	1.1744,	0.816,	-0.0402,	0.3351,	288.5,	2,	1.9,	-0.8,	-0.3],
	[16.5,	97.5,	99893,	305,	-4,	27,	-5,	299.8,	0.2,	1.1,	-1,	-0.3,	16.46,	-3.28,	-1.24,	-0.59,	-0.6,	-4.9,	1,	0,	0.5,	0.2,	-41.8,	87.86,	1.2831,	-0.0004,	0.001,	-0.0007,	-0.0003,	0.6108,	-0.071,	-0.035,	-0.0008,	-0.0033,	2.9101,	0.916,	0.5447,	0.0176,	0.1769,	288.1,	2.1,	1.8,	-0.7,	-0.2],
	[16.5,	98.5,	95436,	306,	-1,	15,	-23,	297.7,	0.3,	1.2,	-1,	-0.1,	15.14,	-3.38,	-1.2,	-0.4,	-0.48,	-4.7,	-1,	0.2,	-0.5,	0.1,	-38.59,	490.04,	1.2832,	-0.0014,	0.0005,	-0.0006,	-0.0003,	0.5899,	-0.0703,	-0.0305,	-0.0024,	-0.003,	2.8593,	0.7606,	0.3448,	0.1,	0.0659,	286.3,	1.8,	1.7,	-0.7,	-0.1],
	[16.5,	99.5,	97942,	354,	-12,	24,	-36,	299.4,	-0.2,	1.2,	-1,	0,	15.61,	-2.95,	-0.82,	-0.42,	-0.11,	-6.3,	-0.7,	-0.1,	0,	0.1,	-36.26,	262.11,	1.2832,	-0.0018,	0.0004,	-0.0005,	-0.0003,	0.6033,	-0.0668,	-0.0266,	-0.0031,	-0.003,	2.6468,	0.6531,	0.2439,	0.1153,	0.0632,	287.4,	1.5,	1.6,	-0.7,	0],
	[16.5,	100.5,	95808,	353,	1,	17,	-37,	298.1,	-0.6,	0.9,	-0.9,	-0.1,	15.18,	-2.95,	-0.37,	-0.68,	0.06,	-6,	-1.1,	0.2,	-0.4,	0.3,	-33.82,	455.54,	1.2832,	-0.0023,	0,	-0.0004,	-0.0004,	0.592,	-0.0654,	-0.0257,	-0.0024,	-0.0028,	2.7561,	0.6258,	0.309,	0.0571,	0.0753,	286.5,	1.3,	1.5,	-0.7,	0],
	[16.5,	101.5,	97868,	385,	5,	19,	-42,	298.9,	-1.2,	0.7,	-0.9,	-0.1,	15.41,	-3.04,	-0.26,	-0.79,	0.07,	-6.1,	-0.1,	0,	-0.1,	0.2,	-31.42,	271.77,	1.2829,	-0.0023,	-0.0001,	-0.0004,	-0.0004,	0.6032,	-0.0605,	-0.0221,	-0.003,	-0.0019,	2.5762,	0.4844,	0.2278,	0.025,	0.029,	287.1,	0.9,	1.3,	-0.7,	0],
	[16.5,	102.5,	98587,	417,	8,	21,	-48,	299.3,	-1.8,	0.7,	-0.9,	-0.1,	15.13,	-3.16,	-0.41,	-0.78,	0.05,	-5.9,	-0.1,	0.3,	-0.3,	0.3,	-29.93,	207.82,	1.283,	-0.0026,	-0.0001,	-0.0004,	-0.0004,	0.6073,	-0.0576,	-0.0196,	-0.0038,	-0.0014,	2.4839,	0.4113,	0.1597,	0.0231,	0.0049,	287.3,	0.5,	1.2,	-0.7,	0],
	[16.5,	103.5,	98749,	423,	12,	24,	-46,	299.2,	-2.1,	0.8,	-1,	-0.1,	14.91,	-3.21,	-0.68,	-0.67,	0.01,	-5.9,	0.1,	0.3,	-0.3,	0.5,	-27.92,	192.95,	1.283,	-0.0026,	0,	-0.0005,	-0.0004,	0.6081,	-0.0562,	-0.0188,	-0.005,	-0.0016,	2.4574,	0.4045,	0.1076,	0.052,	0.0007,	287.3,	0.4,	1.1,	-0.6,	0],
	[16.5,	104.5,	98715,	423,	21,	27,	-43,	298.8,	-2.2,	0.8,	-1.1,	-0.2,	15,	-3.27,	-0.8,	-0.59,	-0.04,	-5.8,	0,	0.3,	-0.3,	0.5,	-25.1,	197.36,	1.2828,	-0.0025,	0,	-0.0005,	-0.0004,	0.6072,	-0.0547,	-0.0184,	-0.0062,	-0.0023,	2.4863,	0.4093,	0.0903,	0.077,	0.0054,	287.1,	0.3,	1.1,	-0.6,	0],
	[16.5,	105.5,	98489,	431,	33,	28,	-40,	298.2,	-2.4,	0.6,	-1.1,	-0.3,	15.21,	-3.18,	-0.8,	-0.59,	-0.11,	-6,	-0.1,	0.2,	-0.3,	0.3,	-21.37,	220.51,	1.2825,	-0.0026,	0,	-0.0005,	-0.0004,	0.6049,	-0.0531,	-0.0185,	-0.0071,	-0.0033,	2.5407,	0.427,	0.1138,	0.0858,	0.0214,	286.9,	0.1,	1,	-0.6,	0],
	[16.5,	106.5,	95336,	428,	49,	11,	-44,	296.1,	-3,	0.3,	-1,	-0.3,	14.86,	-2.71,	-0.72,	-0.55,	-0.22,	-5.8,	0,	0.3,	-0.2,	0.4,	-17.37,	507.83,	1.2824,	-0.0033,	-0.0003,	-0.0004,	-0.0003,	0.588,	-0.0516,	-0.0191,	-0.0078,	-0.0043,	2.7548,	0.5394,	0.1612,	0.1241,	0.036,	285.6,	-0.2,	0.9,	-0.4,	0],
	[15.5,	96.5,	100888,	310,	18,	18,	11,	300.5,	-0.4,	0.1,	-0.5,	-0.8,	17.34,	-2.01,	-0.15,	-0.78,	-0.23,	-7.3,	1.6,	0.6,	0.4,	0.9,	-44.67,	0,	1.2833,	-0.0005,	0.0008,	-0.0005,	-0.0003,	0.617,	-0.0696,	-0.0364,	-0.0015,	-0.0046,	3.1986,	1.2626,	0.9218,	-0.0226,	0.4054,	288.3,	2,	1.7,	-0.5,	-0.2],
	[15.5,	97.5,	100096,	280,	11,	20,	6,	299.8,	0.2,	0.4,	-0.6,	-0.5,	16.93,	-2.28,	-0.39,	-0.71,	-0.34,	-5.8,	0.6,	-0.2,	0.7,	0.3,	-41,	71.56,	1.2831,	-0.0003,	0.0007,	-0.0005,	-0.0003,	0.6144,	-0.0683,	-0.0336,	-0.0021,	-0.0046,	2.9904,	1.0521,	0.6537,	0.0363,	0.2413,	288,	2.1,	1.7,	-0.5,	-0.1],
	[15.5,	98.5,	95582,	277,	8,	13,	-12,	297.4,	0.4,	0.8,	-0.7,	-0.2,	15.32,	-2.99,	-0.61,	-0.65,	-0.31,	-4.7,	-1,	0.5,	-0.6,	0.2,	-37.49,	478.4,	1.2832,	-0.0012,	0.0004,	-0.0004,	-0.0003,	0.5938,	-0.0668,	-0.0291,	-0.0036,	-0.0038,	2.8503,	0.7326,	0.4102,	0.0583,	0.1005,	286.1,	1.8,	1.6,	-0.5,	0],
	[15.5,	99.5,	98585,	328,	-5,	24,	-27,	299.7,	-0.2,	1.1,	-0.9,	0,	15.53,	-2.7,	-0.58,	-0.55,	-0.07,	-7,	-0.8,	-0.2,	0.1,	0.2,	-35.11,	206.63,	1.2832,	-0.0015,	0.0004,	-0.0004,	-0.0003,	0.6097,	-0.0626,	-0.0249,	-0.0043,	-0.0032,	2.5365,	0.5681,	0.2173,	0.0973,	0.0564,	287.5,	1.5,	1.6,	-0.6,	0.1],
	[15.5,	100.5,	99407,	349,	-3,	28,	-29,	300.1,	-0.8,	0.9,	-0.8,	-0.1,	15.75,	-2.67,	-0.37,	-0.79,	0.05,	-6,	-0.7,	0.3,	-0.4,	0.4,	-32.75,	130.99,	1.2834,	-0.0018,	0.0003,	-0.0003,	-0.0004,	0.6136,	-0.0597,	-0.0236,	-0.0038,	-0.0028,	2.5272,	0.506,	0.2179,	0.0444,	0.053,	287.8,	1.2,	1.5,	-0.6,	0],
	[15.5,	101.5,	98458,	351,	8,	23,	-30,	299.3,	-1,	0.7,	-0.7,	-0.1,	15.56,	-2.78,	-0.12,	-0.94,	0.05,	-5.4,	0.1,	0.5,	-0.6,	0.4,	-30.32,	217.87,	1.2832,	-0.002,	0.0001,	-0.0003,	-0.0004,	0.6083,	-0.0572,	-0.0218,	-0.0039,	-0.0021,	2.5582,	0.4477,	0.2366,	0.0071,	0.0357,	287.3,	1,	1.4,	-0.6,	0],
	[15.5,	102.5,	99131,	372,	11,	23,	-36,	299.7,	-1.5,	0.8,	-0.8,	0,	15.39,	-2.86,	-0.3,	-0.92,	-0.03,	-5.4,	0.2,	0.6,	-0.6,	0.4,	-28.64,	158.8,	1.2832,	-0.0021,	0.0001,	-0.0003,	-0.0003,	0.612,	-0.0546,	-0.0198,	-0.0044,	-0.0014,	2.4776,	0.388,	0.1762,	0.009,	0.0076,	287.5,	0.7,	1.3,	-0.6,	0],
	[15.5,	103.5,	99439,	376,	14,	26,	-36,	299.8,	-1.6,	0.9,	-0.8,	0,	15.26,	-2.92,	-0.58,	-0.81,	-0.13,	-5.9,	0.3,	0.4,	-0.3,	0.4,	-26.58,	130.95,	1.2833,	-0.0021,	0.0002,	-0.0003,	-0.0003,	0.6139,	-0.0534,	-0.019,	-0.0053,	-0.0015,	2.443,	0.3742,	0.124,	0.0348,	-0.0066,	287.6,	0.6,	1.2,	-0.5,	0],
	[15.5,	104.5,	99441,	369,	21,	28,	-34,	299.7,	-1.5,	0.9,	-0.9,	0,	15.33,	-3.07,	-0.74,	-0.71,	-0.22,	-5.9,	0.2,	0.2,	-0.2,	0.3,	-23.62,	131.27,	1.2832,	-0.0019,	0.0003,	-0.0004,	-0.0003,	0.6136,	-0.052,	-0.0187,	-0.0061,	-0.002,	2.4527,	0.36,	0.1013,	0.06,	-0.0114,	287.5,	0.6,	1.2,	-0.5,	0],
	[15.5,	105.5,	98236,	351,	33,	26,	-30,	298.7,	-1.3,	0.9,	-1,	-0.1,	15.21,	-3.1,	-0.79,	-0.6,	-0.28,	-6,	-0.1,	0.1,	-0.2,	0.1,	-19.83,	241.12,	1.2829,	-0.0018,	0.0002,	-0.0004,	-0.0004,	0.6067,	-0.0504,	-0.0189,	-0.0068,	-0.003,	2.5055,	0.3659,	0.1093,	0.0882,	-0.0027,	287,	0.6,	1.1,	-0.5,	-0.1],
	[15.5,	106.5,	93010,	320,	52,	8,	-31,	295.4,	-1.7,	0.8,	-1,	-0.2,	14.47,	-2.75,	-0.73,	-0.47,	-0.3,	-5.4,	-0.3,	0.1,	-0.2,	0.1,	-15.32,	721.58,	1.2829,	-0.0023,	-0.0001,	-0.0003,	-0.0004,	0.5789,	-0.0497,	-0.0201,	-0.0073,	-0.0044,	2.8258,	0.5074,	0.1623,	0.1604,	0.0205,	284.8,	0.4,	1,	-0.4,	0],
	[14.5,	96.5,	100902,	284,	24,	14,	18,	300.7,	-0.3,	0,	-0.3,	-0.7,	17.24,	-1.88,	-0.02,	-0.7,	-0.18,	-8.5,	0.8,	0.2,	0.3,	0.5,	-44.11,	0,	1.2833,	-0.0005,	0.0007,	-0.0003,	-0.0003,	0.6202,	-0.0655,	-0.0341,	-0.0031,	-0.0062,	3.0399,	1.0733,	0.7751,	0.011,	0.344,	288.1,	1.9,	1.6,	-0.4,	-0.1],
	[14.5,	97.5,	100561,	258,	18,	17,	16,	300.2,	0.1,	0.2,	-0.4,	-0.6,	17.08,	-1.84,	-0.03,	-0.69,	-0.15,	-7.8,	0.4,	0,	0.3,	0.2,	-40.18,	31.51,	1.2831,	-0.0003,	0.0006,	-0.0003,	-0.0004,	0.6201,	-0.0641,	-0.0317,	-0.0037,	-0.0061,	2.8915,	0.9716,	0.5943,	0.0666,	0.2405,	287.9,	2,	1.6,	-0.4,	-0.1],
	[14.5,	98.5,	97178,	255,	12,	15,	-1,	298.1,	0.4,	0.6,	-0.6,	-0.3,	15.81,	-2.62,	-0.15,	-0.8,	-0.13,	-5.2,	-0.9,	0.6,	-0.5,	0.5,	-36.58,	334.03,	1.2831,	-0.0008,	0.0004,	-0.0003,	-0.0004,	0.6052,	-0.0618,	-0.028,	-0.0045,	-0.0051,	2.7624,	0.6722,	0.4256,	0.041,	0.14,	286.5,	1.8,	1.5,	-0.4,	0],
	[14.5,	99.5,	98841,	296,	2,	23,	-16,	299.7,	-0.2,	0.9,	-0.7,	-0.1,	15.56,	-2.54,	-0.32,	-0.71,	-0.05,	-6.5,	-0.1,	0,	0.1,	0.1,	-34.15,	184.07,	1.2833,	-0.0013,	0.0004,	-0.0003,	-0.0003,	0.6141,	-0.0584,	-0.0246,	-0.0051,	-0.0043,	2.5021,	0.5053,	0.2421,	0.0674,	0.0746,	287.3,	1.4,	1.5,	-0.5,	0.1],
	[14.5,	100.5,	100088,	316,	4,	27,	-17,	300.4,	-0.8,	0.8,	-0.6,	-0.1,	15.91,	-2.42,	-0.17,	-0.89,	0.03,	-5.5,	0.9,	-0.2,	0.1,	-0.1,	-31.72,	70.9,	1.2835,	-0.0015,	0.0003,	-0.0002,	-0.0004,	0.6202,	-0.0557,	-0.0233,	-0.0049,	-0.0037,	2.4888,	0.4631,	0.2373,	0.0295,	0.0697,	287.9,	1.2,	1.5,	-0.5,	0.1],
	[14.5,	101.5,	97747,	301,	17,	19,	-17,	298.8,	-0.8,	0.6,	-0.5,	-0.1,	15.54,	-2.48,	0.03,	-0.96,	0.01,	-5.5,	0.1,	0.6,	-0.7,	0.4,	-28.85,	281.53,	1.2833,	-0.0016,	0.0001,	-0.0002,	-0.0004,	0.6076,	-0.0544,	-0.0222,	-0.0049,	-0.003,	2.5951,	0.4442,	0.2739,	0.0119,	0.0575,	286.9,	1.1,	1.4,	-0.4,	0],
	[14.5,	102.5,	98477,	317,	18,	21,	-23,	299.2,	-1,	0.7,	-0.5,	0,	15.56,	-2.61,	-0.2,	-0.95,	-0.13,	-5.3,	0.2,	0.6,	-0.6,	0.4,	-26.39,	216.41,	1.2833,	-0.0017,	0.0002,	-0.0002,	-0.0003,	0.611,	-0.0521,	-0.0207,	-0.0051,	-0.0023,	2.5351,	0.3818,	0.2062,	0.0145,	0.0186,	287.2,	1,	1.3,	-0.4,	0],
	[14.5,	103.5,	99484,	322,	19,	24,	-25,	299.9,	-1,	0.9,	-0.6,	0.1,	15.61,	-2.74,	-0.5,	-0.93,	-0.29,	-6,	0.2,	0.4,	-0.4,	0.2,	-23.99,	126.4,	1.2833,	-0.0016,	0.0004,	-0.0002,	-0.0003,	0.6163,	-0.0504,	-0.02,	-0.0053,	-0.0022,	2.4655,	0.3366,	0.1501,	0.0205,	-0.0072,	287.6,	0.9,	1.3,	-0.4,	0],
	[14.5,	104.5,	99471,	313,	24,	25,	-23,	299.9,	-0.7,	1,	-0.6,	0.1,	15.65,	-2.92,	-0.7,	-0.85,	-0.43,	-6.1,	-0.1,	0.1,	-0.3,	0,	-21.01,	127.66,	1.2833,	-0.0014,	0.0004,	-0.0002,	-0.0002,	0.6161,	-0.0488,	-0.0198,	-0.0055,	-0.0024,	2.4636,	0.3072,	0.1206,	0.036,	-0.0233,	287.6,	0.9,	1.3,	-0.4,	0],
	[14.5,	105.5,	99112,	302,	32,	26,	-19,	299.6,	-0.4,	1.1,	-0.7,	0,	15.66,	-3.1,	-0.85,	-0.75,	-0.49,	-5.9,	-0.3,	-0.1,	-0.3,	-0.1,	-17.27,	161,	1.2831,	-0.0011,	0.0004,	-0.0003,	-0.0003,	0.6139,	-0.0466,	-0.0198,	-0.0056,	-0.0031,	2.4769,	0.2852,	0.1076,	0.0525,	-0.0197,	287.4,	0.8,	1.2,	-0.4,	-0.1],
	[14.5,	106.5,	97020,	296,	46,	23,	-15,	298.1,	-0.5,	1,	-0.9,	-0.2,	15.22,	-3.03,	-0.92,	-0.55,	-0.39,	-6,	-0.7,	-0.2,	-0.2,	0.1,	-12.68,	351.14,	1.283,	-0.0013,	0.0003,	-0.0003,	-0.0004,	0.6022,	-0.0448,	-0.0205,	-0.0057,	-0.0043,	2.5563,	0.3086,	0.1249,	0.0891,	0.0172,	286.5,	0.6,	1.2,	-0.4,	-0.1],
	[13.5,	96.5,	100912,	257,	25,	12,	24,	300.8,	-0.2,	0.1,	-0.3,	-0.6,	17.14,	-1.81,	-0.02,	-0.64,	-0.19,	-9.1,	0.2,	-0.1,	0.2,	0.2,	-43.54,	0,	1.2832,	-0.0005,	0.0006,	-0.0002,	-0.0003,	0.6236,	-0.0603,	-0.032,	-0.0041,	-0.0077,	2.8702,	0.8667,	0.6097,	0.0376,	0.2775,	288,	1.8,	1.6,	-0.3,	0],
	[13.5,	97.5,	100921,	238,	20,	15,	23,	300.6,	0.1,	0.2,	-0.3,	-0.6,	17.19,	-1.62,	0.09,	-0.63,	-0.08,	-8.5,	0.3,	-0.1,	0.2,	0.2,	-39.42,	0,	1.283,	-0.0003,	0.0005,	-0.0002,	-0.0004,	0.6253,	-0.0589,	-0.0299,	-0.0046,	-0.0074,	2.7962,	0.8507,	0.5118,	0.086,	0.2208,	287.9,	1.9,	1.5,	-0.3,	0],
	[13.5,	98.5,	98084,	228,	15,	15,	10,	298.5,	0.3,	0.5,	-0.5,	-0.4,	16.21,	-2.23,	0.03,	-0.8,	-0.06,	-5.8,	0.4,	0.1,	0.1,	0.3,	-35.9,	252.98,	1.2831,	-0.0006,	0.0004,	-0.0002,	-0.0004,	0.6129,	-0.0568,	-0.0272,	-0.0051,	-0.0066,	2.7367,	0.6518,	0.4201,	0.0541,	0.161,	286.7,	1.7,	1.5,	-0.3,	0],
	[13.5,	99.5,	99039,	258,	10,	21,	-3,	299.5,	-0.1,	0.7,	-0.6,	-0.2,	15.94,	-2.26,	-0.08,	-0.78,	-0.04,	-5.9,	0.3,	0.1,	0.1,	0.2,	-33.05,	167.27,	1.2832,	-0.001,	0.0004,	-0.0002,	-0.0004,	0.6181,	-0.0541,	-0.0248,	-0.0055,	-0.0057,	2.5463,	0.5007,	0.2959,	0.0567,	0.107,	287.2,	1.4,	1.5,	-0.3,	0.1],
	[13.5,	100.5,	100813,	282,	14,	25,	-4,	300.7,	-0.7,	0.5,	-0.4,	-0.2,	16.56,	-2.15,	0.08,	-0.93,	0.03,	-6,	0.9,	-0.1,	0.1,	0.1,	-30.6,	7.85,	1.2834,	-0.0012,	0.0003,	-0.0001,	-0.0004,	0.6267,	-0.0516,	-0.0235,	-0.0055,	-0.005,	2.5368,	0.4585,	0.2959,	0.0248,	0.1056,	288,	1.3,	1.4,	-0.3,	0.1],
	[13.5,	101.5,	99852,	270,	23,	21,	-3,	300,	-0.6,	0.4,	-0.3,	-0.2,	16.51,	-2.28,	0.15,	-1,	0.01,	-6.1,	0.6,	0.3,	-0.1,	0.1,	-27.38,	93.67,	1.2832,	-0.0011,	0.0002,	-0.0001,	-0.0004,	0.6213,	-0.05,	-0.0223,	-0.0055,	-0.0043,	2.5741,	0.4119,	0.2931,	0.0095,	0.0874,	287.6,	1.2,	1.3,	-0.3,	0],
	[13.5,	102.5,	99791,	275,	23,	21,	-9,	299.9,	-0.5,	0.6,	-0.4,	-0.1,	16.38,	-2.41,	-0.09,	-1.04,	-0.19,	-5.9,	0.4,	0.4,	-0.2,	0.2,	-24.23,	99.47,	1.2832,	-0.0012,	0.0003,	0,	-0.0003,	0.6204,	-0.0483,	-0.0211,	-0.0053,	-0.0037,	2.5442,	0.3585,	0.2346,	0.0051,	0.0434,	287.6,	1.1,	1.3,	-0.3,	0],
	[13.5,	103.5,	100553,	282,	21,	24,	-13,	300.3,	-0.5,	0.8,	-0.3,	0.1,	16.36,	-2.61,	-0.36,	-1.1,	-0.44,	-5.4,	0.7,	0.2,	-0.2,	0.1,	-21.55,	31.62,	1.2833,	-0.0011,	0.0005,	0,	-0.0002,	0.6239,	-0.0466,	-0.0208,	-0.0048,	-0.0033,	2.49,	0.2983,	0.1895,	-0.0092,	0.0059,	287.9,	1,	1.3,	-0.4,	0],
	[13.5,	104.5,	100136,	276,	26,	23,	-12,	300.1,	-0.2,	0.9,	-0.3,	0.1,	16.2,	-2.84,	-0.63,	-1.04,	-0.6,	-5.2,	0.6,	0.1,	-0.3,	-0.1,	-18.5,	68.12,	1.2833,	-0.001,	0.0005,	0,	-0.0002,	0.6214,	-0.0451,	-0.0211,	-0.0042,	-0.0034,	2.485,	0.2553,	0.1542,	-0.0062,	-0.0172,	287.8,	1,	1.3,	-0.4,	0],
	[13.5,	105.5,	100001,	272,	32,	26,	-8,	299.9,	-0.1,	1,	-0.5,	0,	16.24,	-3.04,	-0.85,	-0.92,	-0.62,	-5.4,	-0.1,	-0.1,	0,	-0.1,	-14.76,	80.93,	1.2832,	-0.0008,	0.0006,	-0.0001,	-0.0002,	0.6201,	-0.0429,	-0.0214,	-0.0037,	-0.004,	2.5,	0.2408,	0.1321,	0.0089,	-0.0091,	287.7,	0.9,	1.3,	-0.4,	-0.1],
	[13.5,	106.5,	99519,	276,	40,	30,	-2,	299.4,	-0.2,	1,	-0.8,	-0.3,	15.96,	-3.07,	-0.99,	-0.68,	-0.48,	-5.8,	-0.8,	-0.1,	-0.3,	0.1,	-10.24,	125.7,	1.2831,	-0.0008,	0.0006,	-0.0003,	-0.0004,	0.6168,	-0.0406,	-0.0219,	-0.0034,	-0.0051,	2.495,	0.2473,	0.1366,	0.0364,	0.0331,	287.5,	0.7,	1.3,	-0.5,	-0.1],
	[12.5,	96.5,	100919,	229,	22,	11,	29,	300.9,	-0.1,	0.2,	-0.2,	-0.5,	17.17,	-1.61,	-0.03,	-0.55,	-0.21,	-9.2,	0,	-0.2,	0,	0.1,	-42.75,	0,	1.2831,	-0.0005,	0.0006,	-0.0002,	-0.0003,	0.6269,	-0.0543,	-0.0299,	-0.0046,	-0.0089,	2.7611,	0.7159,	0.4975,	0.058,	0.2373,	287.8,	1.6,	1.5,	-0.2,	0],
	[12.5,	97.5,	100922,	217,	19,	14,	28,	300.8,	0,	0.2,	-0.3,	-0.5,	17.28,	-1.43,	0.1,	-0.56,	-0.12,	-8.8,	0.3,	-0.1,	0.2,	0.1,	-38.13,	0,	1.283,	-0.0004,	0.0005,	-0.0001,	-0.0004,	0.6284,	-0.0532,	-0.0283,	-0.0049,	-0.0084,	2.7418,	0.7374,	0.4527,	0.0911,	0.2047,	287.8,	1.7,	1.5,	-0.2,	0],
	[12.5,	98.5,	99523,	206,	16,	16,	21,	299.4,	0.2,	0.4,	-0.4,	-0.4,	16.76,	-1.79,	0.04,	-0.68,	-0.11,	-6.2,	0,	0,	-0.1,	0.3,	-34.72,	124.9,	1.283,	-0.0004,	0.0005,	-0.0001,	-0.0004,	0.6229,	-0.0512,	-0.0264,	-0.0052,	-0.0078,	2.6841,	0.6161,	0.3837,	0.0745,	0.1601,	287.1,	1.6,	1.5,	-0.2,	0.1],
	[12.5,	99.5,	97987,	222,	19,	17,	9,	298.6,	0,	0.6,	-0.5,	-0.3,	16.15,	-1.97,	0,	-0.75,	-0.07,	-6.2,	0,	0.1,	-0.1,	0.2,	-31.5,	261.26,	1.2832,	-0.0009,	0.0005,	-0.0001,	-0.0004,	0.6158,	-0.0499,	-0.025,	-0.0054,	-0.0072,	2.6701,	0.52,	0.3529,	0.0553,	0.1516,	286.5,	1.4,	1.4,	-0.3,	0.1],
	[12.5,	100.5,	100877,	250,	23,	21,	9,	300.8,	-0.6,	0.3,	-0.4,	-0.3,	17.12,	-1.63,	0.2,	-0.77,	0,	-7.4,	0.6,	0.1,	0.1,	0.1,	-29.05,	3.2,	1.2832,	-0.001,	0.0004,	-0.0001,	-0.0004,	0.6295,	-0.0477,	-0.0237,	-0.0053,	-0.0065,	2.6273,	0.5164,	0.3439,	0.0552,	0.144,	287.8,	1.2,	1.4,	-0.2,	0.1],
	[12.5,	101.5,	100420,	238,	30,	19,	12,	300.5,	-0.5,	0.2,	-0.2,	-0.3,	17.14,	-1.67,	0.17,	-0.79,	-0.03,	-7.4,	0.8,	0,	0.2,	0.1,	-26,	43.45,	1.2831,	-0.0009,	0.0003,	0,	-0.0005,	0.627,	-0.0465,	-0.0229,	-0.0051,	-0.0059,	2.6371,	0.4854,	0.3202,	0.0485,	0.1225,	287.7,	1.2,	1.3,	-0.2,	0],
	[12.5,	102.5,	98576,	226,	31,	16,	6,	299.2,	-0.1,	0.4,	-0.2,	-0.2,	16.47,	-1.98,	0.01,	-0.9,	-0.17,	-6.4,	0.1,	0.4,	-0.5,	0.1,	-22.43,	208.82,	1.2831,	-0.0009,	0.0003,	0,	-0.0004,	0.6175,	-0.0456,	-0.0222,	-0.0047,	-0.0054,	2.631,	0.4044,	0.2855,	0.0206,	0.0908,	286.9,	1.2,	1.3,	-0.3,	0],
	[12.5,	103.5,	98902,	234,	30,	17,	-1,	299.4,	0,	0.7,	-0.2,	0,	16.25,	-2.23,	-0.25,	-0.95,	-0.42,	-5.5,	-0.4,	0.2,	-0.5,	-0.1,	-19.42,	180.36,	1.2831,	-0.001,	0.0004,	0,	-0.0003,	0.6183,	-0.0438,	-0.0221,	-0.0038,	-0.0049,	2.553,	0.3208,	0.2291,	0.0045,	0.0453,	287.1,	1.1,	1.3,	-0.3,	0],
	[12.5,	104.5,	100229,	245,	30,	20,	-3,	300.1,	0,	0.8,	-0.2,	0.1,	16.49,	-2.57,	-0.53,	-1,	-0.65,	-5,	0.2,	0,	-0.2,	0,	-16.42,	60.86,	1.2832,	-0.0008,	0.0005,	0.0001,	-0.0002,	0.6244,	-0.0414,	-0.0223,	-0.0028,	-0.0049,	2.497,	0.2461,	0.1843,	-0.0127,	0.01,	287.7,	0.9,	1.3,	-0.4,	0],
	[12.5,	105.5,	100358,	242,	34,	22,	2,	299.9,	0,	0.9,	-0.4,	-0.1,	16.59,	-2.84,	-0.77,	-0.93,	-0.67,	-4.7,	0.1,	-0.2,	-0.3,	-0.1,	-12.67,	49.03,	1.2831,	-0.0006,	0.0006,	0,	-0.0002,	0.6246,	-0.0391,	-0.0227,	-0.0019,	-0.0054,	2.5021,	0.2083,	0.1579,	-0.0138,	0.0122,	287.7,	0.8,	1.3,	-0.4,	0],
	[12.5,	106.5,	99161,	233,	40,	22,	6,	299.2,	-0.1,	1,	-0.7,	-0.3,	16.07,	-2.82,	-0.9,	-0.68,	-0.53,	-5.7,	-0.7,	-0.1,	-0.2,	0.1,	-8.22,	157.46,	1.2831,	-0.0005,	0.0007,	-0.0002,	-0.0003,	0.6182,	-0.0368,	-0.0231,	-0.0012,	-0.0063,	2.501,	0.2236,	0.1617,	0.0135,	0.0517,	287.2,	0.7,	1.3,	-0.5,	-0.1],
	[11.5,	96.5,	100923,	201,	17,	10,	32,	300.9,	-0.1,	0.3,	-0.2,	-0.4,	17.24,	-1.38,	-0.02,	-0.46,	-0.25,	-9.2,	0,	-0.2,	0,	0,	-41.47,	0,	1.283,	-0.0006,	0.0006,	-0.0001,	-0.0003,	0.6301,	-0.0479,	-0.0278,	-0.0048,	-0.0099,	2.6919,	0.5934,	0.4292,	0.0624,	0.2193,	287.7,	1.5,	1.5,	-0.1,	0.1],
	[11.5,	97.5,	100924,	194,	15,	14,	31,	300.8,	-0.1,	0.3,	-0.2,	-0.4,	17.33,	-1.27,	0.06,	-0.5,	-0.2,	-9,	0.2,	-0.1,	0.1,	0,	-36.3,	0,	1.283,	-0.0005,	0.0006,	-0.0001,	-0.0003,	0.6314,	-0.047,	-0.0268,	-0.0047,	-0.0095,	2.6827,	0.6101,	0.4053,	0.0785,	0.1958,	287.6,	1.5,	1.4,	-0.2,	0.1],
	[11.5,	98.5,	100056,	186,	13,	16,	27,	299.8,	0,	0.5,	-0.3,	-0.3,	17.08,	-1.47,	-0.03,	-0.58,	-0.23,	-7.1,	-0.5,	-0.2,	0.4,	0,	-32.7,	77.21,	1.2829,	-0.0004,	0.0006,	-0.0001,	-0.0004,	0.6282,	-0.0454,	-0.0255,	-0.0046,	-0.009,	2.6589,	0.543,	0.3565,	0.0657,	0.1626,	287.2,	1.4,	1.4,	-0.2,	0.1],
	[11.5,	99.5,	99565,	200,	18,	18,	21,	299.4,	-0.1,	0.5,	-0.4,	-0.3,	16.98,	-1.48,	0,	-0.61,	-0.18,	-6.4,	-0.2,	0,	-0.2,	0.2,	-29.56,	120.55,	1.283,	-0.0007,	0.0006,	-0.0001,	-0.0003,	0.6259,	-0.0441,	-0.0246,	-0.0045,	-0.0084,	2.6724,	0.5014,	0.3517,	0.0535,	0.1642,	287,	1.2,	1.4,	-0.2,	0.1],
	[11.5,	100.5,	100921,	218,	26,	19,	21,	300.8,	-0.4,	0.3,	-0.3,	-0.4,	17.45,	-1.22,	0.15,	-0.57,	-0.1,	-8.5,	0.4,	0,	0.1,	0.1,	-27.05,	0,	1.2831,	-0.0009,	0.0004,	-0.0001,	-0.0004,	0.6321,	-0.043,	-0.024,	-0.0042,	-0.0078,	2.6518,	0.4991,	0.3583,	0.0535,	0.1656,	287.7,	1.1,	1.3,	-0.2,	0.1],
	[11.5,	101.5,	100920,	213,	32,	17,	23,	301,	-0.4,	0.2,	-0.2,	-0.4,	17.37,	-1.26,	0.1,	-0.57,	-0.12,	-8.9,	0.4,	0,	0.1,	0,	-24.19,	0,	1.283,	-0.0009,	0.0004,	-0.0001,	-0.0004,	0.632,	-0.0422,	-0.0235,	-0.0037,	-0.0073,	2.609,	0.464,	0.3262,	0.0474,	0.1414,	287.7,	1.1,	1.3,	-0.2,	0],
	[11.5,	102.5,	100521,	199,	34,	15,	22,	300.4,	-0.1,	0.2,	-0.1,	-0.3,	17.36,	-1.35,	0.04,	-0.64,	-0.19,	-8.1,	0.3,	-0.1,	0.2,	0,	-20.55,	36.32,	1.2828,	-0.0007,	0.0003,	0,	-0.0004,	0.6298,	-0.041,	-0.0231,	-0.0031,	-0.0069,	2.6068,	0.4418,	0.2956,	0.0413,	0.1166,	287.5,	1.2,	1.3,	-0.3,	0],
	[11.5,	103.5,	97761,	194,	36,	13,	13,	298.4,	0.1,	0.5,	-0.1,	-0.2,	16.42,	-1.95,	-0.18,	-0.8,	-0.35,	-5.7,	-0.1,	0,	-0.4,	-0.1,	-17.14,	282.16,	1.283,	-0.0008,	0.0004,	0,	-0.0004,	0.6159,	-0.0401,	-0.0235,	-0.002,	-0.0067,	2.6403,	0.3248,	0.2773,	-0.0047,	0.0972,	286.4,	1,	1.3,	-0.3,	0],
	[11.5,	104.5,	99999,	212,	34,	17,	10,	299.7,	0,	0.8,	-0.3,	0,	16.75,	-2.21,	-0.51,	-0.79,	-0.59,	-5.4,	0,	-0.1,	0.1,	0,	-14.16,	82.06,	1.283,	-0.0007,	0.0005,	0,	-0.0002,	0.6264,	-0.0372,	-0.0233,	-0.0011,	-0.0066,	2.5125,	0.2406,	0.202,	-0.0066,	0.049,	287.4,	0.9,	1.3,	-0.4,	0],
	[11.5,	105.5,	100766,	214,	35,	18,	13,	300,	0,	0.9,	-0.4,	-0.1,	16.8,	-2.42,	-0.76,	-0.74,	-0.66,	-4.8,	0.1,	-0.4,	0,	0,	-10.51,	13.3,	1.283,	-0.0006,	0.0006,	0,	-0.0002,	0.63,	-0.0346,	-0.0234,	-0.0001,	-0.0068,	2.4617,	0.1842,	0.1689,	-0.0157,	0.0403,	287.7,	0.7,	1.3,	-0.4,	0],
	[11.5,	106.5,	100020,	198,	39,	17,	16,	299.7,	0.1,	0.9,	-0.5,	-0.3,	16.49,	-2.47,	-0.85,	-0.62,	-0.58,	-5.7,	-0.2,	-0.4,	0.2,	0,	-6.05,	80.37,	1.283,	-0.0004,	0.0007,	-0.0001,	-0.0003,	0.6262,	-0.032,	-0.0236,	0.0009,	-0.0072,	2.4588,	0.1798,	0.1792,	-0.0182,	0.061,	287.4,	0.7,	1.3,	-0.5,	0],
	[10.5,	96.5,	100925,	172,	10,	9,	34,	301,	-0.2,	0.4,	-0.2,	-0.3,	17.3,	-1.14,	-0.01,	-0.38,	-0.3,	-9.3,	-0.1,	-0.1,	0,	0,	-39.66,	0,	1.283,	-0.0006,	0.0006,	-0.0001,	-0.0002,	0.6332,	-0.0413,	-0.0258,	-0.0045,	-0.0108,	2.6409,	0.486,	0.3799,	0.0558,	0.213,	287.5,	1.3,	1.4,	-0.1,	0.2],
	[10.5,	97.5,	100926,	168,	9,	12,	33,	300.9,	-0.2,	0.4,	-0.2,	-0.3,	17.4,	-1.05,	0.04,	-0.4,	-0.28,	-9.2,	0.1,	-0.1,	0,	0,	-34.36,	0,	1.283,	-0.0005,	0.0006,	-0.0001,	-0.0003,	0.6341,	-0.0405,	-0.0252,	-0.0043,	-0.0104,	2.6409,	0.4989,	0.3694,	0.0602,	0.1968,	287.5,	1.2,	1.4,	-0.1,	0.2],
	[10.5,	98.5,	100004,	165,	9,	15,	29,	299.8,	-0.2,	0.5,	-0.3,	-0.2,	17.22,	-1.13,	-0.08,	-0.46,	-0.35,	-7.2,	-0.4,	-0.2,	0,	-0.1,	-30.54,	81.63,	1.2829,	-0.0006,	0.0007,	-0.0002,	-0.0003,	0.6304,	-0.0394,	-0.0245,	-0.0038,	-0.01,	2.6525,	0.4653,	0.3397,	0.0469,	0.171,	287,	1.2,	1.4,	-0.2,	0.1],
	[10.5,	99.5,	100657,	180,	14,	18,	28,	300.2,	-0.3,	0.5,	-0.3,	-0.2,	17.51,	-1.06,	-0.02,	-0.49,	-0.31,	-7.6,	0.4,	-0.1,	0,	0,	-27.44,	24.01,	1.2829,	-0.0008,	0.0006,	-0.0002,	-0.0003,	0.6333,	-0.038,	-0.024,	-0.0032,	-0.0095,	2.6411,	0.4299,	0.3418,	0.0255,	0.1671,	287.3,	1.1,	1.3,	-0.2,	0.1],
	[10.5,	100.5,	100926,	191,	23,	18,	29,	300.8,	-0.4,	0.3,	-0.2,	-0.3,	17.52,	-1,	0.1,	-0.47,	-0.22,	-9,	0.2,	0,	0,	0,	-24.85,	0,	1.2829,	-0.001,	0.0005,	-0.0002,	-0.0003,	0.6344,	-0.0373,	-0.024,	-0.0026,	-0.009,	2.6253,	0.4104,	0.364,	0.0108,	0.1759,	287.5,	1,	1.3,	-0.2,	0.1],
	[10.5,	101.5,	100925,	189,	29,	16,	30,	301.1,	-0.4,	0.2,	-0.1,	-0.3,	17.29,	-1.11,	0.02,	-0.48,	-0.24,	-9.3,	0.1,	-0.1,	0,	0,	-21.95,	0,	1.2829,	-0.0009,	0.0005,	-0.0002,	-0.0004,	0.6343,	-0.0367,	-0.0239,	-0.0019,	-0.0086,	2.5621,	0.3723,	0.3327,	0.0048,	0.154,	287.5,	1,	1.3,	-0.2,	0],
	[10.5,	102.5,	100928,	184,	32,	16,	31,	301,	-0.3,	0.2,	-0.1,	-0.4,	17.33,	-1.13,	-0.02,	-0.5,	-0.24,	-8.9,	0.2,	-0.1,	0,	0,	-18.41,	0,	1.2829,	-0.0008,	0.0005,	-0.0001,	-0.0004,	0.634,	-0.0358,	-0.0239,	-0.001,	-0.0082,	2.5556,	0.3763,	0.3169,	0.0045,	0.1416,	287.6,	0.9,	1.3,	-0.3,	0],
	[10.5,	103.5,	100518,	179,	34,	16,	28,	300.5,	-0.2,	0.4,	-0.1,	-0.3,	17.32,	-1.31,	-0.15,	-0.55,	-0.33,	-8,	0.3,	-0.1,	0.1,	-0.1,	-14.81,	36.54,	1.2828,	-0.0007,	0.0005,	-0.0001,	-0.0004,	0.6317,	-0.0344,	-0.0239,	0,	-0.008,	2.5735,	0.3544,	0.2929,	0.0005,	0.1253,	287.4,	0.9,	1.3,	-0.3,	0],
	[10.5,	104.5,	100446,	185,	33,	17,	23,	300.1,	-0.2,	0.7,	-0.3,	-0.2,	17.08,	-1.66,	-0.46,	-0.55,	-0.51,	-6.4,	0.5,	-0.3,	0.2,	-0.1,	-11.64,	42.5,	1.2829,	-0.0007,	0.0006,	-0.0001,	-0.0003,	0.6313,	-0.0323,	-0.024,	0.001,	-0.0078,	2.5228,	0.2657,	0.2386,	-0.0122,	0.0937,	287.4,	0.7,	1.3,	-0.4,	0],
	[10.5,	105.5,	100874,	189,	35,	17,	24,	300.1,	-0.2,	0.8,	-0.4,	-0.2,	16.96,	-1.85,	-0.69,	-0.5,	-0.61,	-5.5,	0.2,	-0.6,	0.3,	-0.1,	-8.04,	4.52,	1.2829,	-0.0007,	0.0007,	-0.0002,	-0.0002,	0.6333,	-0.0299,	-0.024,	0.0022,	-0.0078,	2.4673,	0.2053,	0.2138,	-0.0285,	0.0794,	287.5,	0.6,	1.3,	-0.5,	0],
	[10.5,	106.5,	100807,	186,	43,	15,	28,	300.1,	-0.3,	0.6,	-0.4,	-0.5,	17.1,	-1.77,	-0.58,	-0.49,	-0.54,	-6.4,	0.2,	-0.5,	0.3,	0.1,	-3.91,	10.77,	1.283,	-0.0005,	0.0008,	-0.0002,	-0.0003,	0.6325,	-0.0274,	-0.0246,	0.0037,	-0.0079,	2.5327,	0.2411,	0.2848,	-0.058,	0.1076,	287.5,	0.5,	1.2,	-0.5,	-0.1],
	[9.5,	96.5,	100924,	144,	3,	8,	35,	301,	-0.2,	0.4,	-0.2,	-0.2,	17.33,	-0.92,	-0.01,	-0.29,	-0.34,	-9.3,	0,	-0.1,	0,	0,	-37.38,	0,	1.2829,	-0.0007,	0.0006,	-0.0001,	-0.0002,	0.6363,	-0.0344,	-0.0235,	-0.0039,	-0.0114,	2.5898,	0.3758,	0.3288,	0.0462,	0.2066,	287.4,	1,	1.3,	-0.1,	0.2],
	[9.5,	97.5,	100925,	139,	0,	11,	34,	300.9,	-0.2,	0.5,	-0.2,	-0.2,	17.42,	-0.87,	-0.04,	-0.3,	-0.37,	-9,	0,	-0.1,	0,	0,	-32.29,	0,	1.2829,	-0.0006,	0.0007,	-0.0001,	-0.0002,	0.6369,	-0.0337,	-0.023,	-0.0036,	-0.011,	2.5915,	0.3805,	0.3109,	0.0493,	0.1856,	287.4,	1,	1.3,	-0.1,	0.2],
	[9.5,	98.5,	99593,	136,	1,	14,	30,	299.5,	-0.2,	0.7,	-0.3,	-0.1,	17.12,	-0.91,	-0.19,	-0.32,	-0.47,	-6.9,	0.2,	-0.3,	-0.2,	-0.2,	-28.14,	118.1,	1.2829,	-0.0007,	0.0007,	-0.0002,	-0.0002,	0.631,	-0.0328,	-0.0227,	-0.003,	-0.0106,	2.6257,	0.3645,	0.2906,	0.0395,	0.1616,	286.8,	1,	1.3,	-0.2,	0.2],
	[9.5,	99.5,	100496,	154,	8,	17,	31,	299.9,	-0.4,	0.5,	-0.3,	-0.1,	17.44,	-0.82,	-0.11,	-0.39,	-0.48,	-7.1,	0.2,	-0.1,	-0.1,	-0.1,	-24.73,	38.34,	1.2828,	-0.0009,	0.0007,	-0.0002,	-0.0002,	0.6347,	-0.0316,	-0.0227,	-0.002,	-0.0102,	2.6044,	0.3374,	0.3081,	0.0014,	0.1525,	287.2,	0.8,	1.3,	-0.2,	0.2],
	[9.5,	100.5,	100929,	168,	19,	17,	34,	300.8,	-0.6,	0.3,	-0.1,	-0.3,	17.49,	-0.76,	0.04,	-0.42,	-0.39,	-8.8,	0.1,	-0.1,	-0.1,	0,	-22.11,	0,	1.2829,	-0.0011,	0.0006,	-0.0002,	-0.0002,	0.6363,	-0.031,	-0.0233,	-0.001,	-0.0098,	2.5942,	0.327,	0.3504,	-0.0251,	0.1671,	287.4,	0.7,	1.2,	-0.3,	0.1],
	[9.5,	101.5,	100928,	168,	25,	16,	36,	301.1,	-0.6,	0.2,	-0.1,	-0.3,	17.26,	-0.86,	0,	-0.4,	-0.35,	-9.3,	0,	-0.1,	0,	-0.1,	-19.36,	0,	1.283,	-0.001,	0.0006,	-0.0002,	-0.0003,	0.6361,	-0.0306,	-0.0236,	-0.0001,	-0.0095,	2.5514,	0.3131,	0.3473,	-0.0318,	0.1637,	287.4,	0.7,	1.2,	-0.3,	0.1],
	[9.5,	102.5,	100931,	167,	29,	16,	37,	301.1,	-0.6,	0.2,	-0.1,	-0.4,	17.21,	-0.93,	-0.04,	-0.41,	-0.33,	-9.3,	0,	-0.1,	0,	-0.1,	-15.92,	0,	1.2829,	-0.0009,	0.0006,	-0.0002,	-0.0003,	0.6358,	-0.0299,	-0.0238,	0.001,	-0.0091,	2.5395,	0.3155,	0.3455,	-0.0381,	0.1592,	287.4,	0.7,	1.2,	-0.3,	0],
	[9.5,	103.5,	100933,	167,	32,	17,	37,	301,	-0.6,	0.2,	-0.1,	-0.4,	17.29,	-0.97,	-0.1,	-0.4,	-0.34,	-9.1,	0.2,	-0.2,	0,	-0.1,	-12.25,	0,	1.2829,	-0.0008,	0.0006,	-0.0002,	-0.0004,	0.6357,	-0.0287,	-0.024,	0.0021,	-0.0088,	2.5475,	0.3224,	0.3395,	-0.0377,	0.1542,	287.4,	0.6,	1.2,	-0.4,	0],
	[9.5,	104.5,	100880,	167,	32,	17,	36,	300.6,	-0.6,	0.3,	-0.3,	-0.3,	17.26,	-1.1,	-0.27,	-0.37,	-0.42,	-8.3,	0.5,	-0.2,	0.2,	-0.1,	-8.84,	4.73,	1.2829,	-0.0008,	0.0007,	-0.0002,	-0.0003,	0.6355,	-0.027,	-0.024,	0.0032,	-0.0084,	2.5395,	0.295,	0.3176,	-0.0428,	0.1384,	287.4,	0.5,	1.2,	-0.5,	0],
	[9.5,	105.5,	100698,	168,	35,	17,	35,	300.2,	-0.7,	0.4,	-0.4,	-0.4,	17.16,	-1.25,	-0.4,	-0.37,	-0.5,	-7.3,	0.3,	-0.4,	0.3,	-0.1,	-5.42,	20.88,	1.2829,	-0.0007,	0.0007,	-0.0003,	-0.0003,	0.6346,	-0.025,	-0.0242,	0.0046,	-0.0081,	2.5442,	0.2654,	0.321,	-0.0646,	0.1308,	287.3,	0.4,	1.2,	-0.5,	0],
	[9.5,	106.5,	100897,	174,	44,	15,	37,	300.1,	-0.9,	0.1,	-0.4,	-0.6,	17.37,	-1.15,	-0.25,	-0.38,	-0.45,	-7.7,	0.2,	-0.2,	0.2,	0.1,	-1.6,	3.63,	1.2829,	-0.0007,	0.0007,	-0.0003,	-0.0003,	0.6353,	-0.0228,	-0.0247,	0.0061,	-0.008,	2.611,	0.3025,	0.3985,	-0.0983,	0.1544,	287.3,	0.2,	1.2,	-0.6,	-0.1],
	[8.5,	96.5,	100923,	119,	-4,	5,	35,	301.1,	-0.3,	0.5,	-0.1,	-0.1,	17.35,	-0.76,	0,	-0.22,	-0.37,	-9.3,	-0.1,	-0.1,	0,	0,	-34.84,	0,	1.2829,	-0.0008,	0.0007,	-0.0001,	-0.0001,	0.6395,	-0.0275,	-0.0211,	-0.0032,	-0.0116,	2.5427,	0.2651,	0.2813,	0.0376,	0.1986,	287.2,	0.8,	1.2,	-0.1,	0.2],
	[8.5,	97.5,	100921,	115,	-6,	9,	35,	301.1,	-0.3,	0.5,	-0.2,	-0.1,	17.4,	-0.73,	-0.05,	-0.23,	-0.42,	-9,	0.1,	-0.1,	0,	0,	-30.19,	0,	1.283,	-0.0007,	0.0007,	-0.0001,	-0.0001,	0.6397,	-0.0269,	-0.0207,	-0.0029,	-0.0112,	2.548,	0.2713,	0.2648,	0.0395,	0.1779,	287.3,	0.8,	1.2,	-0.1,	0.2],
	[8.5,	98.5,	99670,	111,	-6,	12,	32,	299.7,	-0.3,	0.7,	-0.3,	0,	17.1,	-0.8,	-0.26,	-0.21,	-0.55,	-6.8,	0.3,	-0.3,	-0.2,	-0.3,	-25.7,	110.9,	1.283,	-0.0007,	0.0008,	-0.0002,	-0.0001,	0.6339,	-0.0262,	-0.0205,	-0.0023,	-0.0107,	2.5819,	0.2559,	0.2376,	0.0368,	0.15,	286.7,	0.7,	1.2,	-0.2,	0.2],
	[8.5,	99.5,	98816,	124,	3,	14,	32,	299,	-0.5,	0.6,	-0.3,	-0.1,	16.82,	-0.64,	-0.2,	-0.25,	-0.54,	-6.9,	-0.1,	-0.3,	-0.2,	-0.3,	-21.78,	187.4,	1.283,	-0.001,	0.0007,	-0.0003,	-0.0002,	0.6293,	-0.0256,	-0.021,	-0.0012,	-0.0104,	2.6065,	0.2727,	0.2773,	0.0053,	0.1559,	286.3,	0.6,	1.2,	-0.2,	0.2],
	[8.5,	100.5,	100902,	148,	14,	16,	37,	300.7,	-0.8,	0.2,	-0.1,	-0.2,	17.36,	-0.59,	0.01,	-0.36,	-0.49,	-8.3,	0.1,	-0.1,	-0.1,	-0.1,	-19.02,	2.34,	1.283,	-0.0012,	0.0006,	-0.0003,	-0.0002,	0.6381,	-0.0247,	-0.0218,	0.0001,	-0.0102,	2.553,	0.2531,	0.3252,	-0.043,	0.158,	287.3,	0.5,	1.2,	-0.3,	0.1],
	[8.5,	101.5,	100929,	151,	20,	16,	41,	301.1,	-0.8,	0.1,	-0.1,	-0.3,	17.25,	-0.65,	0.02,	-0.33,	-0.41,	-9.4,	0,	-0.1,	0,	-0.1,	-16.61,	0,	1.283,	-0.0011,	0.0006,	-0.0003,	-0.0003,	0.6379,	-0.0242,	-0.0223,	0.0012,	-0.0099,	2.5507,	0.257,	0.3563,	-0.0557,	0.1721,	287.3,	0.4,	1.1,	-0.3,	0.1],
	[8.5,	102.5,	100933,	151,	24,	17,	43,	301.1,	-0.8,	0.1,	-0.1,	-0.4,	17.22,	-0.72,	-0.02,	-0.33,	-0.39,	-9.3,	0,	-0.1,	0,	-0.1,	-13.38,	0,	1.283,	-0.001,	0.0007,	-0.0003,	-0.0003,	0.6378,	-0.0235,	-0.0226,	0.0023,	-0.0094,	2.5482,	0.2655,	0.3668,	-0.066,	0.1716,	287.3,	0.4,	1.1,	-0.4,	0],
	[8.5,	103.5,	100937,	151,	28,	17,	44,	300.9,	-0.9,	0,	-0.2,	-0.5,	17.24,	-0.75,	-0.05,	-0.32,	-0.37,	-9.3,	0,	-0.1,	0,	0,	-9.72,	0,	1.283,	-0.0009,	0.0007,	-0.0003,	-0.0003,	0.6378,	-0.0226,	-0.0228,	0.0034,	-0.0088,	2.5536,	0.2768,	0.3788,	-0.0754,	0.1721,	287.3,	0.3,	1.1,	-0.4,	0],
	[8.5,	104.5,	100916,	151,	31,	17,	45,	300.7,	-1,	0,	-0.3,	-0.5,	17.29,	-0.77,	-0.07,	-0.3,	-0.37,	-9.1,	0.2,	-0.1,	0.1,	0,	-6.13,	2.14,	1.2829,	-0.0009,	0.0007,	-0.0003,	-0.0003,	0.6379,	-0.0212,	-0.0228,	0.0046,	-0.0083,	2.5655,	0.283,	0.3882,	-0.0835,	0.1701,	287.2,	0.2,	1.1,	-0.5,	-0.1],
	[8.5,	105.5,	100921,	151,	34,	16,	44,	300.5,	-1.1,	-0.1,	-0.3,	-0.6,	17.33,	-0.8,	-0.1,	-0.29,	-0.39,	-8.8,	0.3,	0,	0.1,	0.1,	-2.74,	1.84,	1.2829,	-0.0008,	0.0007,	-0.0003,	-0.0003,	0.6381,	-0.0195,	-0.0229,	0.0059,	-0.0079,	2.5779,	0.2804,	0.3982,	-0.0965,	0.1664,	287.2,	0.1,	1.1,	-0.6,	-0.1],
	[8.5,	106.5,	100942,	152,	38,	14,	44,	300.4,	-1.3,	-0.2,	-0.4,	-0.6,	17.38,	-0.79,	-0.09,	-0.28,	-0.4,	-8.9,	0.3,	0,	0.1,	0.1,	1.08,	0,	1.2829,	-0.0008,	0.0007,	-0.0003,	-0.0004,	0.6384,	-0.0179,	-0.0232,	0.0072,	-0.0077,	2.5914,	0.2755,	0.4141,	-0.114,	0.1645,	287.1,	0,	1.1,	-0.6,	-0.1],
	[7.5,	96.5,	100923,	99,	-10,	2,	34,	301.1,	-0.3,	0.5,	-0.1,	0,	17.37,	-0.61,	0.02,	-0.16,	-0.39,	-9.3,	0,	-0.1,	0,	-0.1,	-32.8,	0,	1.2829,	-0.0009,	0.0007,	-0.0001,	-0.0001,	0.643,	-0.021,	-0.0187,	-0.0025,	-0.0115,	2.5062,	0.166,	0.2433,	0.0302,	0.1904,	287.1,	0.5,	1.1,	-0.1,	0.3],
	[7.5,	97.5,	100919,	99,	-10,	5,	34,	301.2,	-0.4,	0.5,	-0.1,	0,	17.39,	-0.59,	0.01,	-0.19,	-0.41,	-9.2,	0.1,	-0.1,	-0.1,	-0.1,	-28.31,	0,	1.283,	-0.0008,	0.0007,	-0.0001,	-0.0001,	0.6427,	-0.0206,	-0.0184,	-0.0022,	-0.0112,	2.5166,	0.1803,	0.2426,	0.0297,	0.1834,	287.1,	0.5,	1.1,	-0.1,	0.2],
	[7.5,	98.5,	100844,	97,	-11,	9,	34,	300.9,	-0.4,	0.6,	-0.2,	0,	17.4,	-0.61,	-0.11,	-0.2,	-0.51,	-8.4,	0.2,	-0.1,	-0.1,	-0.1,	-23.66,	6.34,	1.283,	-0.0007,	0.0008,	-0.0002,	-0.0001,	0.6419,	-0.0199,	-0.0182,	-0.0018,	-0.0107,	2.5241,	0.1852,	0.2287,	0.024,	0.1609,	287.1,	0.5,	1.1,	-0.1,	0.2],
	[7.5,	99.5,	99735,	102,	-5,	12,	34,	299.9,	-0.5,	0.6,	-0.2,	0,	17.09,	-0.58,	-0.16,	-0.2,	-0.56,	-7,	0.2,	-0.2,	-0.2,	-0.3,	-19.43,	104.68,	1.2831,	-0.0009,	0.0008,	-0.0002,	-0.0001,	0.6362,	-0.0193,	-0.0186,	-0.0009,	-0.0103,	2.5596,	0.1917,	0.2439,	0.0083,	0.153,	286.7,	0.4,	1.1,	-0.2,	0.2],
	[7.5,	100.5,	100700,	124,	4,	15,	38,	300.5,	-0.8,	0.3,	-0.2,	-0.1,	17.28,	-0.49,	-0.04,	-0.27,	-0.54,	-7.8,	0,	-0.2,	-0.1,	-0.2,	-16.21,	19.73,	1.283,	-0.0011,	0.0007,	-0.0003,	-0.0002,	0.6399,	-0.0185,	-0.0192,	0.0003,	-0.0102,	2.5292,	0.1929,	0.2839,	-0.0311,	0.1515,	287.1,	0.2,	1.1,	-0.3,	0.1],
	[7.5,	101.5,	100929,	136,	14,	17,	43,	300.8,	-0.9,	0.1,	-0.1,	-0.3,	17.32,	-0.55,	0.03,	-0.29,	-0.47,	-8.8,	-0.1,	-0.1,	-0.1,	-0.1,	-13.8,	0,	1.2831,	-0.0012,	0.0007,	-0.0003,	-0.0002,	0.6405,	-0.018,	-0.0199,	0.0015,	-0.0098,	2.5406,	0.1884,	0.3285,	-0.0604,	0.1642,	287.2,	0.1,	1,	-0.3,	0.1],
	[7.5,	102.5,	100935,	137,	20,	17,	47,	300.9,	-0.9,	0,	-0.2,	-0.4,	17.29,	-0.63,	0.02,	-0.28,	-0.43,	-9.2,	0,	-0.1,	0,	-0.1,	-10.9,	0,	1.283,	-0.0011,	0.0006,	-0.0004,	-0.0003,	0.6403,	-0.0174,	-0.0204,	0.0027,	-0.0093,	2.5455,	0.1918,	0.3504,	-0.075,	0.1705,	287.2,	0.1,	1,	-0.4,	0],
	[7.5,	103.5,	100940,	136,	24,	18,	49,	300.8,	-1,	-0.1,	-0.2,	-0.5,	17.25,	-0.61,	0,	-0.26,	-0.4,	-9.3,	0.1,	-0.1,	0,	0,	-7.54,	0,	1.283,	-0.001,	0.0007,	-0.0004,	-0.0003,	0.6405,	-0.0166,	-0.0204,	0.0038,	-0.0085,	2.5421,	0.2125,	0.3652,	-0.0827,	0.1734,	287.1,	0,	1,	-0.4,	-0.1],
	[7.5,	104.5,	100942,	134,	26,	17,	50,	300.7,	-1.1,	-0.2,	-0.3,	-0.6,	17.31,	-0.61,	0.01,	-0.25,	-0.38,	-9.2,	0.2,	0,	0,	0.1,	-3.95,	0,	1.2829,	-0.001,	0.0007,	-0.0004,	-0.0003,	0.6407,	-0.0156,	-0.0204,	0.0049,	-0.0079,	2.5529,	0.2245,	0.379,	-0.0915,	0.1728,	287.1,	0,	1,	-0.5,	-0.1],
	[7.5,	105.5,	100943,	132,	28,	16,	50,	300.5,	-1.3,	-0.2,	-0.3,	-0.6,	17.35,	-0.6,	0.01,	-0.23,	-0.37,	-9.2,	0.3,	0.1,	0.1,	0.1,	-0.07,	0,	1.2829,	-0.0009,	0.0007,	-0.0004,	-0.0004,	0.6411,	-0.0144,	-0.0205,	0.006,	-0.0074,	2.554,	0.2229,	0.3802,	-0.1005,	0.1655,	287,	-0.1,	1,	-0.5,	-0.1],
	[7.5,	106.5,	100941,	129,	30,	13,	48,	300.5,	-1.3,	-0.2,	-0.4,	-0.6,	17.37,	-0.58,	-0.01,	-0.21,	-0.38,	-9.1,	0.3,	0.1,	0.1,	0.1,	3.95,	0,	1.2828,	-0.0009,	0.0006,	-0.0004,	-0.0004,	0.6414,	-0.0134,	-0.0206,	0.007,	-0.0072,	2.5424,	0.2067,	0.3673,	-0.1098,	0.1523,	287,	-0.1,	0.9,	-0.6,	-0.1],
	[6.5,	96.5,	100921,	88,	-14,	0,	32,	301.1,	-0.5,	0.4,	0,	0,	17.47,	-0.48,	0.05,	-0.12,	-0.39,	-9.1,	0,	-0.1,	0,	-0.1,	-31.65,	0,	1.2829,	-0.0011,	0.0006,	-0.0001,	-0.0001,	0.6468,	-0.0153,	-0.0161,	-0.0021,	-0.0107,	2.4827,	0.075,	0.2034,	0.0236,	0.1706,	286.9,	0.3,	1,	0,	0.2],
	[6.5,	97.5,	100915,	88,	-13,	2,	33,	301.2,	-0.5,	0.4,	0,	0,	17.46,	-0.51,	0.04,	-0.15,	-0.4,	-9.2,	0,	-0.1,	-0.1,	-0.1,	-26.78,	0,	1.283,	-0.001,	0.0007,	-0.0001,	-0.0001,	0.6463,	-0.015,	-0.0158,	-0.0019,	-0.0106,	2.4814,	0.0879,	0.2064,	0.0198,	0.1691,	287,	0.3,	1,	-0.1,	0.2],
	[6.5,	98.5,	100912,	84,	-13,	5,	34,	301.3,	-0.4,	0.4,	-0.1,	0,	17.46,	-0.52,	0.03,	-0.19,	-0.44,	-9.2,	0.1,	-0.1,	-0.1,	0,	-21.7,	0,	1.283,	-0.0008,	0.0007,	-0.0001,	-0.0001,	0.6456,	-0.0143,	-0.0155,	-0.0016,	-0.0102,	2.4784,	0.111,	0.2103,	0.0148,	0.1594,	287.1,	0.3,	1,	-0.1,	0.2],
	[6.5,	99.5,	100723,	83,	-14,	7,	35,	300.9,	-0.4,	0.5,	-0.1,	0,	17.34,	-0.54,	-0.06,	-0.19,	-0.53,	-8.5,	0.1,	-0.1,	-0.1,	-0.1,	-17.42,	16.55,	1.2831,	-0.0007,	0.0008,	-0.0002,	-0.0001,	0.6442,	-0.0134,	-0.0155,	-0.0011,	-0.0098,	2.473,	0.1172,	0.2005,	0.0139,	0.1368,	287,	0.2,	1,	-0.1,	0.2],
	[6.5,	100.5,	99058,	89,	-8,	10,	35,	299.4,	-0.7,	0.5,	-0.2,	-0.1,	16.65,	-0.46,	-0.18,	-0.16,	-0.53,	-6.6,	0.3,	-0.1,	-0.1,	-0.1,	-13.55,	165.83,	1.2831,	-0.0009,	0.0008,	-0.0002,	-0.0002,	0.6362,	-0.013,	-0.0158,	-0.0002,	-0.0095,	2.4782,	0.1287,	0.2046,	0.0012,	0.1343,	286.2,	0.1,	1,	-0.2,	0.1],
	[6.5,	101.5,	99384,	113,	5,	15,	39,	299.5,	-1.1,	0.2,	-0.2,	-0.2,	16.59,	-0.24,	0,	-0.21,	-0.53,	-6.4,	0,	-0.1,	-0.2,	-0.1,	-10.45,	138.19,	1.2831,	-0.0012,	0.0007,	-0.0004,	-0.0002,	0.6369,	-0.0129,	-0.0168,	0.001,	-0.0092,	2.4604,	0.1669,	0.2694,	-0.0407,	0.1375,	286.4,	-0.1,	0.9,	-0.3,	0.1],
	[6.5,	102.5,	100889,	130,	17,	18,	48,	300.6,	-1.1,	-0.1,	-0.2,	-0.3,	17.21,	-0.43,	0.11,	-0.24,	-0.5,	-8.4,	-0.2,	-0.1,	-0.1,	0,	-7.82,	4.75,	1.283,	-0.0013,	0.0005,	-0.0004,	-0.0003,	0.6432,	-0.0123,	-0.0177,	0.0023,	-0.0088,	2.4787,	0.1336,	0.3008,	-0.0678,	0.1417,	287,	-0.2,	0.9,	-0.3,	0],
	[6.5,	103.5,	100940,	125,	20,	17,	52,	300.7,	-1,	-0.1,	-0.2,	-0.5,	17.36,	-0.59,	0.04,	-0.21,	-0.42,	-9.1,	0,	-0.1,	0,	0,	-5.16,	0,	1.2829,	-0.0011,	0.0006,	-0.0004,	-0.0003,	0.6435,	-0.0115,	-0.0177,	0.0034,	-0.0082,	2.5109,	0.1247,	0.3074,	-0.0718,	0.1595,	286.9,	-0.2,	0.8,	-0.4,	-0.1],
	[6.5,	104.5,	100941,	119,	21,	16,	53,	300.6,	-1.1,	-0.2,	-0.3,	-0.6,	17.34,	-0.54,	0.05,	-0.2,	-0.4,	-9.2,	0.2,	0,	0,	0,	-1.95,	0,	1.2829,	-0.001,	0.0006,	-0.0004,	-0.0004,	0.6437,	-0.0108,	-0.0176,	0.0044,	-0.0075,	2.5087,	0.1439,	0.3177,	-0.0775,	0.1589,	286.9,	-0.2,	0.8,	-0.4,	-0.1],
	[6.5,	105.5,	100940,	113,	21,	14,	52,	300.6,	-1.2,	-0.2,	-0.4,	-0.6,	17.38,	-0.5,	0.05,	-0.18,	-0.38,	-9.1,	0.3,	0,	0.1,	0.1,	2.37,	0,	1.2828,	-0.001,	0.0006,	-0.0004,	-0.0004,	0.644,	-0.0101,	-0.0175,	0.0053,	-0.007,	2.5069,	0.1461,	0.3162,	-0.0868,	0.1496,	286.9,	-0.2,	0.8,	-0.5,	-0.1],
	[6.5,	106.5,	100937,	108,	22,	11,	50,	300.6,	-1.3,	-0.2,	-0.4,	-0.5,	17.4,	-0.45,	0.03,	-0.15,	-0.38,	-9.1,	0.3,	0,	0.1,	0.1,	6.73,	0,	1.2828,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6445,	-0.0095,	-0.0174,	0.0059,	-0.0068,	2.493,	0.1359,	0.302,	-0.094,	0.1366,	286.9,	-0.2,	0.8,	-0.5,	-0.1],
	[5.5,	96.5,	99579,	71,	-14,	-3,	29,	299.7,	-0.6,	0.3,	-0.1,	0,	17.01,	-0.27,	0.1,	-0.08,	-0.34,	-7.3,	0,	-0.2,	-0.1,	-0.2,	-29.21,	123.08,	1.2827,	-0.0012,	0.0006,	-0.0001,	-0.0002,	0.645,	-0.0107,	-0.0133,	-0.0021,	-0.0091,	2.4102,	0.0185,	0.1547,	0.0235,	0.1311,	286.1,	0.1,	0.9,	0,	0.2],
	[5.5,	97.5,	100821,	81,	-13,	1,	31,	300.8,	-0.6,	0.3,	-0.1,	0,	17.35,	-0.38,	0.08,	-0.11,	-0.36,	-8.2,	0,	-0.1,	0,	0,	-24.49,	10.89,	1.2828,	-0.0011,	0.0006,	-0.0002,	-0.0002,	0.6498,	-0.0099,	-0.0131,	-0.0018,	-0.0094,	2.3756,	0.0299,	0.1612,	0.0129,	0.1292,	286.8,	0,	0.9,	0,	0.2],
	[5.5,	98.5,	100915,	76,	-14,	3,	34,	301.1,	-0.4,	0.3,	-0.1,	0,	17.67,	-0.48,	0.08,	-0.14,	-0.38,	-8.9,	0,	0,	-0.1,	0,	-19.64,	0,	1.2828,	-0.0009,	0.0006,	-0.0002,	-0.0002,	0.6493,	-0.0093,	-0.0127,	-0.0015,	-0.0093,	2.4225,	0.0448,	0.1708,	0.0061,	0.1286,	286.9,	0,	0.8,	-0.1,	0.2],
	[5.5,	99.5,	100908,	68,	-17,	2,	34,	301.1,	-0.3,	0.4,	-0.1,	0,	17.6,	-0.48,	0.05,	-0.16,	-0.44,	-8.7,	0.1,	0,	0,	0,	-15.42,	0,	1.2829,	-0.0007,	0.0007,	-0.0001,	-0.0002,	0.6486,	-0.0085,	-0.0124,	-0.0012,	-0.009,	2.4162,	0.0452,	0.155,	0.0105,	0.1111,	287,	0.1,	0.9,	-0.1,	0.2],
	[5.5,	100.5,	98290,	65,	-14,	4,	33,	299.1,	-0.4,	0.5,	-0.1,	0,	16.63,	-0.51,	-0.16,	-0.15,	-0.51,	-6.8,	0.3,	-0.1,	0,	-0.1,	-11.01,	234.81,	1.283,	-0.0009,	0.0007,	-0.0002,	-0.0003,	0.6361,	-0.0084,	-0.0125,	-0.0006,	-0.0087,	2.4371,	0.0243,	0.1242,	0.0091,	0.0992,	285.8,	0,	0.8,	-0.1,	0.1],
	[5.5,	101.5,	94927,	78,	1,	8,	34,	297,	-0.9,	0.3,	-0.3,	-0.1,	15.42,	-0.09,	-0.06,	-0.12,	-0.51,	-5.6,	0.1,	-0.1,	-0.1,	-0.2,	-6.84,	541.32,	1.2831,	-0.0012,	0.0006,	-0.0003,	-0.0003,	0.6199,	-0.0089,	-0.0136,	0.0004,	-0.0082,	2.4875,	0.1303,	0.1986,	-0.0136,	0.1202,	284.3,	-0.1,	0.8,	-0.2,	0],
	[5.5,	102.5,	98941,	113,	14,	16,	43,	299.1,	-1.3,	-0.1,	-0.3,	-0.2,	16.35,	0.01,	0.21,	-0.16,	-0.56,	-6.1,	-0.4,	-0.1,	-0.1,	-0.1,	-3.88,	178.82,	1.283,	-0.0014,	0.0004,	-0.0004,	-0.0003,	0.6379,	-0.0087,	-0.0147,	0.0017,	-0.008,	2.3877,	0.1457,	0.2597,	-0.0482,	0.1105,	286,	-0.3,	0.7,	-0.3,	0],
	[5.5,	103.5,	100869,	116,	18,	16,	51,	300.4,	-1,	-0.2,	-0.3,	-0.4,	17.36,	-0.49,	0.07,	-0.18,	-0.49,	-8.5,	-0.1,	-0.1,	-0.1,	0,	-1.86,	6.56,	1.2827,	-0.0012,	0.0004,	-0.0004,	-0.0003,	0.6463,	-0.0073,	-0.0148,	0.0028,	-0.0078,	2.4428,	0.0646,	0.2423,	-0.0592,	0.1279,	286.8,	-0.4,	0.7,	-0.4,	-0.1],
	[5.5,	104.5,	100938,	105,	16,	14,	53,	300.6,	-1,	-0.1,	-0.3,	-0.5,	17.41,	-0.53,	0.04,	-0.16,	-0.4,	-9.1,	0.1,	0,	0,	0,	0.68,	0,	1.2828,	-0.0011,	0.0005,	-0.0004,	-0.0004,	0.6466,	-0.0066,	-0.0146,	0.0037,	-0.0073,	2.4618,	0.0691,	0.2502,	-0.0649,	0.1436,	286.8,	-0.3,	0.7,	-0.4,	-0.1],
	[5.5,	105.5,	100936,	96,	15,	12,	52,	300.7,	-1.1,	-0.2,	-0.4,	-0.5,	17.41,	-0.43,	0.07,	-0.13,	-0.38,	-9.2,	0.3,	0,	0,	0.1,	4.67,	0,	1.2828,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6468,	-0.0062,	-0.0142,	0.0042,	-0.0068,	2.4588,	0.0805,	0.2526,	-0.0712,	0.1357,	286.8,	-0.3,	0.7,	-0.4,	-0.1],
	[5.5,	106.5,	100932,	90,	14,	9,	50,	300.6,	-1.2,	-0.2,	-0.4,	-0.5,	17.42,	-0.35,	0.07,	-0.1,	-0.37,	-9.2,	0.3,	0,	0.1,	0.1,	9.12,	0,	1.2828,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6474,	-0.0058,	-0.0139,	0.0046,	-0.0065,	2.4485,	0.0776,	0.2418,	-0.0757,	0.1243,	286.8,	-0.3,	0.7,	-0.4,	-0.1],
	[4.5,	96.5,	92266,	36,	-13,	-6,	23,	295.2,	-0.4,	0.3,	-0.1,	-0.1,	15.2,	0,	0.18,	-0.09,	-0.32,	-5.4,	-0.2,	0,	-0.1,	-0.2,	-27.12,	792,	1.2827,	-0.0011,	0.0005,	-0.0001,	-0.0004,	0.6168,	-0.0092,	-0.0115,	-0.0021,	-0.0074,	2.5494,	0.0029,	0.14,	0.0347,	0.0971,	282.6,	0.2,	0.8,	0,	0.1],
	[4.5,	97.5,	95312,	54,	-10,	-2,	26,	297.2,	-0.6,	0.3,	-0.1,	-0.1,	15.89,	-0.14,	0.03,	-0.07,	-0.34,	-6,	-0.1,	0,	-0.1,	-0.1,	-21.86,	508.88,	1.2826,	-0.001,	0.0006,	-0.0002,	-0.0003,	0.6283,	-0.0068,	-0.0107,	-0.002,	-0.0078,	2.3961,	0.0204,	0.1239,	0.026,	0.094,	284.2,	0,	0.7,	0,	0.1],
	[4.5,	98.5,	100897,	70,	-15,	2,	33,	300.8,	-0.4,	0.3,	-0.1,	-0.1,	17.55,	-0.32,	0.08,	-0.1,	-0.34,	-8.1,	0,	0,	0,	0,	-17.95,	4.35,	1.2826,	-0.0009,	0.0006,	-0.0002,	-0.0003,	0.6526,	-0.0052,	-0.0102,	-0.0014,	-0.0082,	2.3098,	0.0226,	0.1326,	0.0067,	0.094,	286.8,	-0.1,	0.7,	-0.1,	0.1],
	[4.5,	99.5,	100906,	64,	-17,	0,	35,	301,	-0.3,	0.3,	-0.1,	-0.1,	18.07,	-0.32,	0.17,	-0.1,	-0.33,	-8.6,	0,	0,	0,	0,	-13.49,	0,	1.2826,	-0.0008,	0.0006,	-0.0002,	-0.0003,	0.6513,	-0.0046,	-0.0099,	-0.0011,	-0.0081,	2.4235,	0.0117,	0.1376,	0.0059,	0.0968,	286.9,	-0.1,	0.7,	-0.1,	0.1],
	[4.5,	100.5,	99961,	57,	-16,	1,	34,	300.3,	-0.3,	0.4,	-0.1,	-0.1,	17.45,	-0.31,	0.08,	-0.12,	-0.41,	-6.9,	-0.2,	0,	0,	-0.2,	-8.95,	85.88,	1.2827,	-0.0009,	0.0005,	-0.0002,	-0.0003,	0.6459,	-0.0043,	-0.0097,	-0.0006,	-0.0079,	2.3867,	-0.015,	0.0994,	0.0064,	0.078,	286.5,	-0.1,	0.7,	-0.1,	0.1],
	[4.5,	101.5,	95271,	62,	-4,	5,	33,	297.3,	-0.6,	0.3,	-0.2,	-0.1,	15.76,	-0.18,	-0.04,	-0.11,	-0.47,	-5.4,	-0.3,	-0.2,	-0.1,	-0.2,	-4.51,	509.08,	1.283,	-0.0011,	0.0005,	-0.0003,	-0.0004,	0.6238,	-0.0052,	-0.0106,	0.0002,	-0.0078,	2.449,	0.0203,	0.1101,	-0.0003,	0.0878,	284.4,	-0.2,	0.7,	-0.2,	0],
	[4.5,	102.5,	96241,	81,	6,	10,	38,	297.6,	-1.1,	0,	-0.3,	-0.2,	15.81,	-0.03,	0.07,	-0.11,	-0.54,	-5.4,	0,	-0.1,	-0.1,	-0.3,	-1.02,	420.68,	1.283,	-0.0012,	0.0005,	-0.0004,	-0.0004,	0.6279,	-0.0046,	-0.0111,	0.0012,	-0.0075,	2.4134,	0.0931,	0.1783,	-0.022,	0.0939,	284.8,	-0.3,	0.6,	-0.3,	0],
	[4.5,	103.5,	100571,	99,	12,	14,	49,	299.8,	-1,	-0.1,	-0.3,	-0.3,	17.27,	-0.38,	0.06,	-0.17,	-0.56,	-7.4,	-0.1,	-0.2,	-0.1,	-0.1,	1.16,	32.93,	1.2826,	-0.0012,	0.0004,	-0.0004,	-0.0003,	0.6474,	-0.0036,	-0.0115,	0.0022,	-0.0074,	2.3834,	0.0309,	0.1788,	-0.0478,	0.0956,	286.5,	-0.5,	0.6,	-0.3,	-0.1],
	[4.5,	104.5,	100935,	90,	10,	12,	53,	300.5,	-0.9,	-0.1,	-0.3,	-0.4,	17.5,	-0.5,	0.02,	-0.15,	-0.42,	-8.8,	0,	-0.1,	-0.1,	0,	3.62,	0,	1.2827,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6491,	-0.0029,	-0.0113,	0.0029,	-0.0071,	2.4258,	0.0135,	0.1907,	-0.0588,	0.125,	286.7,	-0.4,	0.6,	-0.3,	-0.1],
	[4.5,	105.5,	100932,	81,	8,	9,	52,	300.7,	-1,	-0.1,	-0.3,	-0.5,	17.42,	-0.37,	0.07,	-0.11,	-0.38,	-9.3,	0.2,	0,	0,	0,	7.24,	0,	1.2827,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6495,	-0.0026,	-0.0109,	0.0033,	-0.0067,	2.416,	0.0274,	0.1976,	-0.0595,	0.1228,	286.7,	-0.4,	0.6,	-0.3,	-0.1],
	[4.5,	106.5,	100929,	74,	8,	6,	50,	300.7,	-1.1,	-0.1,	-0.3,	-0.4,	17.42,	-0.27,	0.1,	-0.08,	-0.35,	-9.1,	0.2,	0,	0.1,	0,	11.62,	0,	1.2827,	-0.001,	0.0005,	-0.0004,	-0.0004,	0.6502,	-0.0024,	-0.0104,	0.0034,	-0.0063,	2.4088,	0.0267,	0.1892,	-0.0607,	0.1132,	286.6,	-0.4,	0.6,	-0.3,	-0.1]]

	numData = len(data)
	plat = dlat*180.0/pi
	plon = dlon*180.0/pi
	indx=[0,0,0,0]
	flat = floor(plat-0.5)-4
	flon =106-floor(plon-0.5)

	# transform to polar distance in degrees
	ppod = (-dlat + pi/2)*180/pi;
	# find the index (line in the grid file) of the nearest point
	# changed for the 1 degree grid (GP)
	ipod = floor((ppod+1));
	ilon = floor((plon+1));
	# normalized (to one) differences, can be positive or negative
	# changed for the 1 degree grid (GP)
	diffpod = (ppod - (ipod - 0.5));
	difflon = (plon - (ilon - 0.5));

	indx[0] = (numData-flat*11-1-flon) # LL
	indx[1] = indx[0] - 11                     # UL
	indx[2] = indx[0] + 1                      # LR
	indx[3] = indx[1] + 1                      # UR

	# change the reference epoch to January 1 2000
	dmjd1 = dmjd-51544.5;

	#print(dmjd1)
	# mean gravity in m/s**2
	gm = 9.80665;
	# molar mass of dry air in kg/mol
	dMtr = 28.965*0.001;
	# universal gas constant in J/K/mol
	Rg = 8.3143;
	# factors for amplitudes
	if (it==1):# then  constant parameters
		cosfy = 0;
		coshy = 0;
		sinfy = 0;
		sinhy = 0;
	else:
		cosfy = cos(dmjd1/365.25*2*pi);
		coshy = cos(dmjd1/365.25*4*pi);
		sinfy = sin(dmjd1/365.25*2*pi);
		sinhy = sin(dmjd1/365.25*4*pi);

	# initialization
	coorgrid = zeros([numData,2]);
	pgrid = zeros([numData, 5]);
	Tgrid = zeros([numData, 5]);
	Qgrid = zeros([numData, 5]);
	dTgrid = zeros([numData, 5]);
	u = zeros([numData, 1]);
	Hs = zeros([numData, 1]);
	ahgrid = zeros([numData, 5]);
	awgrid = zeros([numData, 5]);
	lagrid = zeros([numData, 5]);
	Tmgrid = zeros([numData, 5]);

	# initialization of new vectors
	p =  0.0
	T =  0.0
	dT = 0.0
	Tm = 0.0
	e =  0.0
	ah = 0.0
	aw = 0.0
	la = 0.0
	undu = 0.0


	# loop over grid points
	# for the 1 degree grid (GP)
	for n in range(0,numData,1):  # Thailand AreaGrid
		# read mean values and amplitudes
		vec = data[n]
		pgrid[n]  = vec[2:7];          # pressure in Pascal
		Tgrid[n]  = vec[7:12];         # temperature in Kelvin
		Qgrid[n]  = array(vec[12:17])/1000;  # specific humidity in kg/kg
		dTgrid[n] = array(vec[17:22])/1000;  # temperature lapse rate in Kelvin/m
		u[n]        = vec[22];           # geoid undulation in m
		Hs[n]         = vec[23];           # orthometric grid height in m
		ahgrid[n] = array(vec[24:29])/1000;  # hydrostatic mapping function coefficient, dimensionless
		awgrid[n] = array(vec[29:34])/1000;  # wet mapping function coefficient, dimensionless
		lagrid[n] = vec[34:39];    	   # water vapor decrease factor, dimensionless
		Tmgrid[n] = vec[39:44];    # mean temperature in Kelvin
	undul = [0.0,0.0,0.0,0.0]
	Ql =  [0.0,0.0,0.0,0.0]
	dTl =  [0.0,0.0,0.0,0.0]
	Tl =  [0.0,0.0,0.0,0.0]
	pl 	=  [0.0,0.0,0.0,0.0]
	ahl =  [0.0,0.0,0.0,0.0]
	awl =  [0.0,0.0,0.0,0.0]
	lal =  [0.0,0.0,0.0,0.0]
	Tml =  [0.0,0.0,0.0,0.0]
	el  = [0.0,0.0,0.0,0.0]
	for l in range(0,4):
		#transforming ellipsoidal height to orthometric height:
		#Hortho = -N + Hell
		undul[l] = u[indx[l]]
		hgt = hell if is_msl==1 else hell-undul[l]
		#pressure, temperature at the height of the grid
		T0 = Tgrid[indx[l],0] + Tgrid[indx[l],1]*cosfy + Tgrid[indx[l],2]*sinfy + Tgrid[indx[l],3]*coshy + Tgrid[indx[l],4]*sinhy
		p0 = pgrid[indx[l],0] + pgrid[indx[l],1]*cosfy + pgrid[indx[l],2]*sinfy + pgrid[indx[l],3]*coshy + pgrid[indx[l],4]*sinhy

		#humidity
		Ql[l]= Qgrid[indx[l],0] + Qgrid[indx[l],1]*cosfy + Qgrid[indx[l],2]*sinfy + Qgrid[indx[l],3]*coshy + Qgrid[indx[l],4]*sinhy

		#reduction = stationheight - gridheight
		Hs1 = Hs[indx[l]]
		redh = hgt - Hs1

		#lapse rate of the temperature in degree / m
		dTl[l] = dTgrid[indx[l],0] + dTgrid[indx[l],1]*cosfy + dTgrid[indx[l],2]*sinfy + dTgrid[indx[l],3]*coshy + dTgrid[indx[l],4]*sinhy

		#temperature reduction to station height
		Tl[l] = T0 + dTl[l]*redh - 273.15

		#virtual temperature
		Tv = T0*(1+0.6077*Ql[l])
		c = gm*dMtr/(Rg*Tv)

		#pressure in hPa
		pl[l] = (p0*exp(-c*redh))/100


		#hydrostatic coefficient ah
		ahl[l] = ahgrid[indx[l],0] + ahgrid[indx[l],1]*cosfy + ahgrid[indx[l],2]*sinfy + ahgrid[indx[l],3]*coshy + ahgrid[indx[l],4]*sinhy;

		#wet coefficient aw
		awl[l] = awgrid[indx[l],0] + awgrid[indx[l],1]*cosfy + awgrid[indx[l],2]*sinfy + awgrid[indx[l],3]*coshy + awgrid[indx[l],4]*sinhy;

		#water vapor decrease factor la - added by GP
		lal[l] = lagrid[indx[l],0] + lagrid[indx[l],1]*cosfy + lagrid[indx[l],2]*sinfy + lagrid[indx[l],3]*coshy + lagrid[indx[l],4]*sinhy;

		#mean temperature of the water vapor Tm - added by GP
		Tml[l] = Tmgrid[indx[l],0] + Tmgrid[indx[l],1]*cosfy + Tmgrid[indx[l],2]*sinfy + Tmgrid[indx[l],3]*coshy + Tmgrid[indx[l],4]*sinhy;

		#water vapor pressure in hPa - changed by GP
		e0 = Ql[l]*p0/(0.622+0.378*Ql[l])/100; # on the grid
		el[l] = e0*pow((100*pl[l]/p0),(lal[l]+1));  # on the station height - [14] Askne and Nordius, 1987


	dnpod1 = fabs(diffpod); # distance nearer point
	dnpod2 = 1 - dnpod1;  # distance to distant point
	dnlon1 = fabs(difflon);
	dnlon2 = 1 - dnlon1;
	#print(plon)
	#print(ilon)
	# pressure
	R1 = dnpod2*pl[0]+dnpod1*pl[1];
	R2 = dnpod2*pl[2]+dnpod1*pl[3];
	p = dnlon2*R1+dnlon1*R2;

	# temperature
	R1 = dnpod2*Tl[0]+dnpod1*Tl[1];
	R2 = dnpod2*Tl[2]+dnpod1*Tl[3];
	T = dnlon2*R1+dnlon1*R2;

	# temperature in degree per km
	R1 = dnpod2*dTl[0]+dnpod1*dTl[1];
	R2 = dnpod2*dTl[2]+dnpod1*dTl[3];
	dT = (dnlon2*R1+dnlon1*R2)*1000;

	# water vapor pressure in hPa - changed by GP
	R1 = dnpod2*el[0]+dnpod1*el[1];
	R2 = dnpod2*el[2]+dnpod1*el[3];
	e = dnlon2*R1+dnlon1*R2;

	# hydrostatic
	R1 = dnpod2*ahl[0]+dnpod1*ahl[1];
	R2 = dnpod2*ahl[2]+dnpod1*ahl[3];
	ah = dnlon2*R1+dnlon1*R2;

	# wet
	R1 = dnpod2*awl[0]+dnpod1*awl[1];
	R2 = dnpod2*awl[2]+dnpod1*awl[3];
	aw = dnlon2*R1+dnlon1*R2;

	# undulation
	R1 = dnpod2*undul[0]+dnpod1*undul[1];
	R2 = dnpod2*undul[2]+dnpod1*undul[3];
	undu = dnlon2*R1+dnlon1*R2;

	# water vapor decrease factor la - added by GP
	R1 = dnpod2*lal[0]+dnpod1*lal[1];
	R2 = dnpod2*lal[2]+dnpod1*lal[3];
	la = dnlon2*R1+dnlon1*R2;

	# mean temperature of the water vapor Tm - added by GP
	R1 = dnpod2*Tml[0]+dnpod1*Tml[1];
	R2 = dnpod2*Tml[2]+dnpod1*Tml[3];
	Tm = dnlon2*R1+dnlon1*R2;

	hort = hgt if is_msl==1 else (hell - undu)[0]
	return (p,T[0],dT,Tm,e,ah,aw,la,hort)


def gpt3_THA (dmjd,dlat,dlon,hell,nstat,it,is_msl:0):
	#function [p,T,dT,Tm,e,ah,aw,la,undu] = gpt2_1w_THA (dmjd,dlat,dlon,hell,nstat,it)

	# (c) Department of Geodesy and Geoinformation, Vienna University of
	# Technology, 2013
	#
	# The copyright in this document is vested in the Department of Geodesy and
	# Geoinformation (GEO), Vienna University of Technology, Austria. This document
	# may only be reproduced in whole or in part, or stored in a retrieval
	# system, or transmitted in any form, or by any means electronic,
	# mechanical, photocopying or otherwise, either with the prior permission
	# of GEO or in accordance with the terms of ESTEC Contract No.
	# 4000107329/12/NL/LvH.
	# ---
	#
	# This subroutine determines pressure, temperature, temperature lapse rate,
	# mean temperature of the water vapor, water vapor pressure, hydrostatic
	# and wet mapping function coefficients ah and aw, water vapour decrease
	# factor and geoid undulation for specific sites near the Earth surface.
	# It is based on a 1 x 1 degree external grid file ('gpt2_1wA.grd') with mean
	# values as well as sine and cosine amplitudes for the annual and
	# semiannual variation of the coefficients.
	#
	# c Reference:
	# J. B?hm, G. M?ller, M. Schindelegger, G. Pain, R. Weber, Development of an
	# improved blind model for slant delays in the troposphere (GPT2w),
	# GPS Solutions, 2014, doi:10.1007/s10291-014-0403-7
	#
	# input parameters:
	#
	# dmjd:  modified Julian date (scalar, only one epoch per call is possible)
	# dlat:  ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
	# dlon:  longitude in radians [-pi:pi] or [0:2pi] (vector)
	# hell:  ellipsoidal height in m (vector)
	# nstat: number of stations in dlat, dlon, and hell
	#        maximum possible: not relevant for Matlab version
	# it:    case 1: no time variation but static quantities
	#        case 0: with time variation (annual and semiannual terms)
	#
	# output parameters:
	#
	# p:    pressure in hPa (vector of length nstat)
	# T:    temperature in degrees Celsius (vector of length nstat)
	# dT:   temperature lapse rate in degrees per km (vector of length nstat)
	# Tm:   mean temperature of the water vapor in degrees Kelvin (vector of length nstat)
	# e:    water vapor pressure in hPa (vector of length nstat)
	# ah:   hydrostatic mapping function coefficient at zero height (VMF1)
	#       (vector of length nstat)
	# aw:   wet mapping function coefficient (VMF1) (vector of length nstat)
	# la:   water vapor decrease factor (vector of length nstat)
	# undu: geoid undulation in m (vector of length nstat)
	# is_msl : 0: ellipsoidal height, 1 : orthometric height (default:0)

	# The hydrostatic mapping function coefficients have to be used with the
	# height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
	# coefficients refer to zero height.
	#
	# Example 1 (Vienna, 2 August 2012, with time variation,grid file 'gpt2_1wA.grd):
	#
	#dmjd = 56141.0;
	#dlat = 10.200*pi/180.0;
	#dlon = 100.400*pi/180.0;
	#hell = 156.0;
	#nstat = 1;
	#it = 0;
	#
	# output:
	# p = 1002.788 hPa
	# T = 22.060 deg Celsius
	# dT = -6.230 deg / km
	# Tm = 281.304 K
	# e = 16.742 hPa
	# ah = 0.0012646
	# aw = 0.0005752
	# la = 2.6530
	# undu = 45.76 m
	#
	# Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):
	#
	# dmjd = 56141.d0
	# dlat(1) = 48.20d0*pi/180.d0
	# dlon(1) = 16.37d0*pi/180.d0
	# hell(1) = 156.d0
	# nstat = 1
	# it = 1
	#
	# output:
	# p = 1003.709 hPa
	# T = 11.79 deg Celsius
	# dT = -5.49 deg / km
	# Tm = 273.22 K
	# e = 10.26 hPa
	# ah = 0.0012396
	# aw = 0.0005753
	# la = 2.6358
	# undu = 45.76 m
	#
	#
	# Klemens Lagler, 2 August 2012
	# Johannes Boehm, 6 August 2012, revision
	# Klemens Lagler, 21 August 2012, epoch change to January 1 2000
	# Johannes Boehm, 23 August 2012, adding possibility to determine constant field
	# Johannes Boehm, 27 December 2012, reference added
	# Johannes Boehm, 10 January 2013, correction for dlat = -90 degrees
	#                                  (problem found by Changyong He)
	# Johannes Boehm, 21 May 2013, bug with dmjd removed (input parameter dmjd was replaced
	#                 unintentionally; problem found by Dennis Ferguson)
	# Gregory Pain,   17 June 2013, adding water vapor decrease factor la
	# Gregory Pain,   21 June 2013, using the 1 degree grid : better for calculating zenith wet delays (la)
	# Gregory Pain,   01 July 2013, adding mean temperature of the water vapor Tm
	# Gregory Pain,   30 July 2013, changing the method to calculate the water vapor partial pressure (e)
	# Gregory Pain,   31 July 2013, correction for (dlat = -90 degrees, dlon = 360 degrees)
	# Johannes Boehm, 27 December 2013, copyright notice added
	# Johannes Boehm, 25 August 2014, default input file changed to
	#                 gpt2_1wA.grd (slightly different humidity values)
	# Johannes Boehm, 25 August 2014, reference changed to Boehm et al. in GPS
	#                 Solutions
	# ---
	#  lat    lon   p:a0    A1   B1   A2   B2  T:a0    A1   B1   A2   B2  Q:a0    A1    B1    A2    B2 dT:a0    A1   B1   A2   B2    undu       Hs   h:a0      A1      B1      A2      B2    w:a0      A1      B1      A2      B2  lam:a0      A1      B1      A2      B2    Tm:a0    A1   B1   A2   B2
	data=[
	[21.5,96.5,92024,400,-7,-13,-45,296.3,-1.5,1.5,-2,-0.3,13.28,-4.34,-3.01,0.23,-0.7,-5.8,-0.6,-0.9,0.3,-0.1,-45.35,801.33,1.27211,-0.00664,-0.00172,-0.00062,-0.00053,0.56213,-0.075,-0.03436,-0.00322,0.00284,3.0122,1.0544,0.2471,0.3619,-0.0198,284.7,0.3,1,-0.8,-0.4,-10.54,-8.62,-6.62,-0.86,-0.05,6.51,0.55,-0.46,0.07,-0.25,-1.7,0.57,6,-0.54,5.65,14.42,2.95,6.42,-0.13,-7.08],
	[21.5,97.5,88076,354,-14,-29,-54,294.2,-1.7,1.4,-2,-0.2,12.55,-4.14,-2.93,0.41,-0.66,-6.4,-1.6,-1,0.2,-0.1,-42.37,1181.94,1.27197,-0.00717,-0.00217,-0.00041,-0.00043,0.54456,-0.07268,-0.02943,-0.00437,0.00073,3.1104,1.0644,0.1526,0.4327,-0.0338,282.9,0.2,0.7,-0.7,-0.3,-10.91,-9.28,-7.07,-0.82,-0.1,4.27,0.39,-0.49,-0.16,-0.22,-3.63,0.12,1.12,1.19,4.55,5.58,4.58,10.31,-0.4,-5.45],
	[21.5,98.5,89055,373,-32,-22,-62,294.6,-1.6,1.4,-1.9,0,12.6,-4.07,-2.81,0.55,-0.55,-5.9,-1.7,-0.8,0.1,-0.1,-40.37,1086.33,1.27195,-0.00708,-0.00188,-0.00058,-0.00029,0.55053,-0.07237,-0.02644,-0.00529,-0.00069,2.8987,0.89,0.0937,0.3752,-0.0106,283.3,0.2,0.7,-0.6,-0.2,-10.8,-9.26,-7.29,-0.78,-0.13,4.19,0.77,-0.4,0.03,-0.16,-2.63,0.2,-1.74,1.42,3.69,1.08,2.97,7.95,0.47,-3.25],
	[21.5,99.5,85734,340,-32,-43,-67,292.7,-1.9,1.2,-1.8,0.2,11.99,-3.75,-2.14,0.33,-0.47,-6,-1.6,-0.8,0.2,-0.1,-38.26,1416.63,1.27185,-0.00755,-0.00237,-0.00025,-0.00024,0.53386,-0.07123,-0.02625,-0.00372,-0.00228,3.1182,1.0291,0.1683,0.3873,-0.004,281.7,-0.1,0.6,-0.6,0,-10.96,-9.62,-7.71,-0.78,-0.08,4.17,0.99,-0.27,-0.38,-0.16,-1.2,0.42,-3.16,1.14,2.27,2.28,6.22,9.46,-1.65,-3.67],
	[21.5,100.5,89119,394,-47,-34,-77,293.8,-2,0.8,-1.5,0.2,12.72,-3.75,-1.55,0.1,-0.28,-4.7,-1.1,0.1,-0.1,0.1,-36.51,1083.01,1.27203,-0.00735,-0.00212,-0.00037,-0.0001,0.55111,-0.06942,-0.02551,-0.00369,-0.00366,2.836,0.7846,0.2264,0.2133,0.0494,283.1,-0.1,0.8,-0.7,0.1,-11.33,-9.24,-8.1,-0.94,-0.14,4.64,0.76,0.05,-0.45,-0.12,-1.16,0.67,-5.24,1.44,1.63,-2.37,5.09,8.54,-0.58,-0.63],
	[21.5,101.5,89464,413,-45,-38,-81,293.5,-2.5,0.3,-1.4,0.2,12.97,-3.77,-1.15,-0.09,-0.17,-4.8,-1.1,0.4,0.1,0.2,-34.78,1051.32,1.27195,-0.00749,-0.00235,-0.00029,-0.00007,0.55202,-0.06859,-0.02468,-0.00248,-0.00407,2.8511,0.7368,0.2313,0.1612,0.051,283.1,-0.4,0.7,-0.7,0.2,-11.46,-9.22,-8.3,-1.05,-0.16,5.01,1.18,0.36,-0.37,-0.05,-0.32,0.85,-4.73,0.38,0.93,-2.13,5.28,9.71,-0.53,-0.99],
	[21.5,102.5,89645,433,-40,-40,-84,293.1,-3,-0.1,-1.4,0.1,13.05,-3.91,-0.97,-0.24,-0.06,-4.9,-0.6,0.4,0.3,0.3,-32.9,1036.26,1.27197,-0.00766,-0.00257,-0.00023,-0.00008,0.55276,-0.06729,-0.0226,-0.00177,-0.00484,2.8503,0.66,0.1983,0.1124,0.0446,283,-0.7,0.5,-0.8,0.2,-11.32,-9.04,-8.23,-1.12,-0.15,5.35,0.84,0.46,-0.32,-0.14,0.55,-0.59,-3.94,0.05,0.4,-3.59,4.58,6.93,-0.02,-0.62],
	[21.5,103.5,90829,483,-35,-42,-86,293.2,-4,-0.4,-1.3,0,13.22,-4.18,-0.95,-0.34,0.01,-4.2,-0.8,0.9,-0.3,0.6,-31.16,925.87,1.27201,-0.00794,-0.00263,-0.00021,-0.00011,0.55914,-0.06505,-0.01986,-0.00124,-0.0055,2.7757,0.5392,0.1497,0.0672,0.0319,283.1,-1.2,0.3,-0.8,0.2,-11.3,-8.84,-8.36,-1.28,-0.12,5.34,1.76,0.04,-0.52,-0.19,-0.03,-1.26,-3.4,0.54,-0.46,-4.53,5.62,6.01,-0.88,-1.22],
	[21.5,104.5,94263,600,-26,-41,-86,294.5,-5.1,-0.9,-1.3,-0.3,13.54,-4.52,-0.93,-0.35,0.1,-5.8,-0.1,0.4,0,0.4,-29.74,608.36,1.27204,-0.00838,-0.00259,-0.00026,-0.00017,0.57717,-0.06043,-0.01668,-0.00165,-0.00692,2.5629,0.3807,0.112,0.038,0.037,284.1,-1.9,0.1,-0.8,0.2,-11.95,-8.73,-8.43,-1.44,0.01,5.27,1.53,0.36,-0.52,-0.12,-2.9,-1.86,-4.53,2.07,-0.86,-10.96,7.69,4.32,-0.8,-2.74],
	[21.5,105.5,99391,768,12,-24,-77,296.7,-5.8,-1.3,-1.2,-0.6,14.31,-5.12,-1.07,-0.37,0.19,-6.2,-0.1,0.1,-0.2,0.1,-28.45,148.95,1.27186,-0.00863,-0.00269,-0.00042,-0.00035,0.60357,-0.05618,-0.01512,-0.00179,-0.00651,2.3658,0.2184,0.0559,0.0097,0.0292,285.8,-2.6,-0.1,-0.8,0.1,-12,-8.35,-8.17,-1.33,0.11,5.91,1.13,1.94,0.22,0.11,-3.55,-4.47,-5.31,3.98,-0.52,-10.37,5.79,2.59,-1,0.31],
	[21.5,106.5,98030,745,28,-28,-77,295.4,-5.9,-1.7,-1.1,-0.7,13.96,-5.24,-1.09,-0.36,0.22,-5.7,0,0.3,0.1,0.3,-26.23,269.42,1.27221,-0.00871,-0.00288,-0.00041,-0.00035,0.59681,-0.05519,-0.01519,-0.00204,-0.00624,2.4573,0.1899,0.0371,0.0053,0.0279,285.1,-2.6,-0.2,-0.7,0.1,-11.36,-8.08,-8.4,-1.33,0.12,6.55,0.65,1.67,0.29,0.03,-0.81,-5.7,-3.71,3.87,-0.48,-5.81,6.59,4.9,-2.28,-0.89],
	[20.5,96.5,90470,359,0,-16,-41,295.5,-1.3,1.5,-1.9,-0.3,13.05,-4.33,-3.17,0.32,-0.84,-6,-0.5,-0.7,0.1,-0.1,-45.2,950.84,1.27359,-0.00576,-0.00133,-0.00056,-0.00056,0.55725,-0.07504,-0.036,-0.00192,0.00227,3.1985,1.1815,0.2994,0.3984,0.0161,284.1,0.7,1.1,-0.8,-0.4,-8.79,-7.06,-6.08,-1.01,-0.59,6.7,0.89,-0.16,-0.01,-0.21,0.91,0.46,3.97,-3.15,2.73,15.19,2.4,9.14,-2.27,-6.42],
	[20.5,97.5,89160,351,-14,-18,-49,295.1,-1.2,1.7,-1.9,-0.1,12.85,-4.28,-3.23,0.59,-0.75,-6.5,-1.5,-0.9,0.1,-0.1,-42.11,1076.81,1.27364,-0.00599,-0.00133,-0.00057,-0.00043,0.55225,-0.07334,-0.03003,-0.00458,0.00008,3.0484,0.9883,0.1091,0.4323,-0.0131,283.6,0.6,1,-0.7,-0.3,-8.94,-7.54,-6.36,-0.93,-0.56,4.41,0.59,-0.28,-0.17,-0.19,1.6,0.97,3.08,-1.7,1.92,5.16,5.07,12.84,-1.47,-5.02],
	[20.5,98.5,90350,372,-32,-10,-58,295.7,-1.3,1.7,-1.9,0.1,13.14,-4.03,-2.89,0.6,-0.65,-4.4,-0.3,0.3,-0.2,0.1,-40.04,961.59,1.27359,-0.00594,-0.00111,-0.00069,-0.00027,0.55872,-0.07306,-0.02721,-0.00506,-0.00158,2.8931,0.8865,0.0894,0.3737,0.0107,284.2,0.6,1,-0.7,-0.1,-9.25,-7.68,-6.61,-0.92,-0.57,4.47,0.99,-0.08,-0.04,-0.24,0.2,0.19,0.41,0.3,1.93,2.3,3.65,12.29,-0.77,-2.43],
	[20.5,99.5,90348,385,-35,-15,-67,295.3,-1.3,1.3,-1.6,0.2,13.3,-3.79,-1.85,0.19,-0.43,-4.9,-1.1,0.1,-0.3,0.1,-37.89,962.84,1.27348,-0.00608,-0.00138,-0.00055,-0.00017,0.55795,-0.07183,-0.02731,-0.00331,-0.00366,2.88,0.8449,0.2241,0.2528,0.0578,284.1,0.5,1.1,-0.7,0,-9.48,-7.64,-6.94,-1.02,-0.6,4.79,0.95,0.12,-0.4,-0.21,0.26,-1.41,-2.22,0.79,1.12,1.46,8.67,14.13,-2.05,-2.2],
	[20.5,100.5,92013,413,-35,-16,-74,295.5,-1.5,0.8,-1.4,0.2,13.77,-3.68,-1.12,-0.14,-0.21,-5,-0.8,0.2,-0.2,0.1,-36.01,806.06,1.27346,-0.00605,-0.00159,-0.00046,-0.00012,0.56619,-0.07005,-0.02645,-0.00234,-0.00489,2.77,0.7394,0.2674,0.1509,0.0765,284.7,0.3,1.1,-0.8,0.2,-9.75,-7.25,-7.26,-1.17,-0.75,5.37,0.87,0.45,-0.29,-0.09,1.04,-1.73,-2.2,0.08,-0.04,-1.75,4.1,8.7,0.13,0.58],
	[20.5,101.5,90648,406,-27,-24,-75,294.3,-2,0.5,-1.3,0.2,13.6,-3.7,-0.83,-0.23,-0.08,-4.5,-0.9,0.4,-0.4,0.2,-34.28,937.84,1.27339,-0.00632,-0.00194,-0.00032,-0.00017,0.55855,-0.06925,-0.02445,-0.00158,-0.00514,2.8687,0.7338,0.2671,0.1347,0.0715,283.9,0.1,0.9,-0.8,0.2,-9.94,-7.21,-7.45,-1.17,-0.73,5.02,1,0.47,-0.43,-0.09,-0.28,-1.94,-1.54,-0.17,-0.47,-2.21,7.12,10.07,-0.23,-0.35],
	[20.5,102.5,90406,416,-23,-31,-78,293.8,-2.4,0.2,-1.3,0.2,13.4,-3.81,-0.74,-0.3,0.01,-4.1,-0.6,0.7,-0.5,0.3,-32.45,964.1,1.27337,-0.00654,-0.00212,-0.00025,-0.00019,0.55678,-0.06683,-0.02148,-0.00137,-0.00541,2.839,0.6411,0.2169,0.102,0.0512,283.6,-0.2,0.7,-0.8,0.2,-10.03,-7.02,-7.45,-1.16,-0.62,5.52,0.7,0.33,-0.21,-0.21,-1.52,-3.18,-1.82,0.39,-1.28,-3.14,5.78,6.48,-0.62,0.45],
	[20.5,103.5,89625,434,-21,-41,-81,292.9,-3.5,0,-1.3,0.1,13.12,-3.93,-0.77,-0.36,0.05,-4.8,-0.1,0.5,0.3,0.4,-30.7,1042.08,1.27338,-0.007,-0.00218,-0.00013,-0.0002,0.55301,-0.06548,-0.01929,-0.00133,-0.00527,2.8775,0.5945,0.1719,0.0885,0.029,283,-0.6,0.6,-0.7,0.2,-10.34,-7.15,-7.4,-1.19,-0.47,5.37,1.7,0.25,-0.39,-0.34,-2.14,-3.67,-2.14,1.9,-1.81,-3.38,5.55,5.99,-2.39,-0.75],
	[20.5,104.5,93540,559,-15,-37,-84,294.4,-4.8,-0.5,-1.1,-0.3,13.55,-4.15,-0.81,-0.46,0.1,-6.2,0.3,0.2,0.1,0.2,-29.12,675.78,1.27345,-0.00755,-0.00207,-0.00017,-0.00022,0.57458,-0.06124,-0.01611,-0.00252,-0.00581,2.6167,0.434,0.109,0.0462,0.024,284.2,-1.4,0.4,-0.7,0.2,-10.61,-6.82,-7.48,-1.29,-0.39,5.92,1.83,0.5,-0.39,-0.25,-2.21,-4.88,-2.03,3.72,-2.08,-9.15,7.72,4.15,-2.19,-3.11],
	[20.5,105.5,98607,730,19,-27,-77,296.5,-5.7,-1.2,-0.9,-0.6,14.46,-4.65,-0.92,-0.54,0.14,-6.1,0.2,0.5,-0.2,0.2,-27.16,217.44,1.2734,-0.00801,-0.00215,-0.00021,-0.00038,0.60072,-0.05713,-0.01468,-0.00319,-0.0055,2.4526,0.282,0.0699,0.0148,0.0218,285.8,-2.2,0.1,-0.7,0,-10.5,-6.31,-7.49,-1.46,-0.37,5.07,0.75,1.82,0.21,0.06,-1.79,-7.31,-0.8,4.98,-1.25,-12.67,7.58,0.95,-1.35,-0.95],
	[20.5,106.5,100517,797,58,-19,-70,297.1,-5.7,-2,-0.6,-0.8,15.2,-5.1,-1.03,-0.44,0.06,-6.2,-0.4,0.9,-0.4,0.2,-24.72,48.55,1.27332,-0.008,-0.00236,-0.00027,-0.00046,0.60894,-0.05481,-0.01431,-0.00362,-0.00464,2.5273,0.2205,0.0496,0.0305,0.003,286.4,-2.3,-0.2,-0.6,0,-9.85,-5.94,-7.27,-1.54,-0.32,6.12,0.25,1.89,0.38,0.03,-0.03,-9.71,1.31,5.12,-0.15,-8.6,8.35,3,-1.94,-0.9],
	[19.5,96.5,93985,365,-5,6,-30,297.4,-0.6,1.7,-1.7,-0.3,14.01,-4.46,-3.29,0.35,-0.99,-5.9,-0.5,-0.6,0,-0.2,-45.43,619.89,1.27494,-0.00435,-0.00037,-0.00081,-0.0005,0.57501,-0.07638,-0.03742,-0.00135,0.00069,3.0356,1.0005,0.3268,0.2593,0.0336,285.8,1.2,1.5,-0.9,-0.4,-7.75,-5.38,-5.48,-0.94,-0.96,7.1,0.81,-0.15,0.1,-0.2,0.95,-2.1,-0.77,-1.89,2.56,15.92,2.43,9.83,-2.73,-5.52],
	[19.5,97.5,92903,360,-20,5,-41,297.4,-0.7,1.8,-1.8,-0.1,13.88,-4.38,-3.42,0.71,-0.86,-5.9,0.6,-0.2,0.3,-0.1,-42.16,721.04,1.27502,-0.00453,-0.00028,-0.00088,-0.00037,0.5716,-0.07459,-0.03165,-0.00376,-0.00078,2.9015,0.8702,0.1294,0.3547,0.0182,285.4,1.1,1.4,-0.7,-0.3,-8.09,-5.5,-5.74,-1.01,-1.02,4.5,1.09,-0.27,-0.14,-0.36,0.69,-1.4,-1.24,-0.66,1.21,7.36,5.98,15.2,-3.15,-4.64],
	[19.5,98.5,91040,364,-26,-2,-52,296.1,-1,1.8,-1.8,0.1,13.63,-4.02,-2.8,0.58,-0.69,-5.9,-1.4,-0.2,-0.1,-0.1,-39.92,896.24,1.27484,-0.00502,-0.0006,-0.00071,-0.00024,0.5629,-0.0741,-0.02855,-0.00369,-0.00199,2.9506,0.8879,0.1172,0.3639,0.0236,284.6,1,1.3,-0.7,-0.1,-7.99,-5.89,-5.91,-0.95,-0.95,4.63,1.13,0.11,-0.08,-0.28,-0.32,-2.01,-2.72,0.49,0.38,4.1,5.85,16.39,-2.49,-1.91],
	[19.5,99.5,93323,398,-33,2,-62,297,-1,1.4,-1.5,0.2,14.31,-3.7,-1.62,0.06,-0.38,-5.1,-0.5,0.3,-0.3,0.1,-37.87,681.21,1.2747,-0.00497,-0.00072,-0.00066,-0.00014,0.57383,-0.07249,-0.02764,-0.00247,-0.00412,2.8019,0.7788,0.253,0.1919,0.0824,285.6,0.9,1.4,-0.8,0.1,-8.22,-5.49,-5.98,-1,-0.99,5.44,0.75,0.56,-0.35,-0.18,0.66,-2.39,-2.5,0.75,-0.56,-0.47,9.19,15.38,-1.93,-0.17],
	[19.5,100.5,93111,403,-23,-4,-65,296.3,-1.2,0.9,-1.3,0.1,14.38,-3.54,-0.84,-0.3,-0.16,-4.6,-0.9,0.3,-0.5,0.1,-35.72,702.61,1.27465,-0.00512,-0.00119,-0.00046,-0.0002,0.57248,-0.07174,-0.02652,-0.0013,-0.00468,2.8053,0.7428,0.3057,0.1212,0.0861,285.4,0.7,1.3,-0.8,0.1,-8.51,-5.33,-6.16,-1.1,-0.97,5.79,0.76,0.7,-0.22,-0.07,1,-2.33,-1.5,0.71,-1.09,-0.2,4.84,9.22,0.01,1.13],
	[19.5,101.5,91478,393,-10,-15,-65,295.1,-1.3,0.6,-1.2,0.1,13.96,-3.54,-0.63,-0.33,-0.02,-5.5,-1.7,0.2,-0.3,0.2,-33.79,859.61,1.27453,-0.00538,-0.00156,-0.00032,-0.00028,0.56385,-0.07086,-0.02415,-0.00092,-0.00465,2.8608,0.7191,0.289,0.1157,0.0736,284.5,0.5,1.1,-0.8,0.1,-8.48,-5.35,-6.29,-1.14,-0.86,4.98,0.68,0.6,-0.43,-0.16,0.74,-2.75,-1.27,1.21,-1.68,-1.86,7.87,9.49,-0.64,0.73],
	[19.5,102.5,89194,375,-2,-29,-66,293.5,-1.9,0.4,-1.3,0,13.39,-3.58,-0.52,-0.39,0.06,-5.4,-1.2,0.2,0.1,0.2,-31.83,1081.76,1.27452,-0.00576,-0.00182,-0.00015,-0.00032,0.55147,-0.06866,-0.02154,-0.00113,-0.00491,2.9833,0.7083,0.2778,0.1073,0.0576,283.3,0.3,0.9,-0.7,0.1,-8.32,-5.46,-6.43,-1.21,-0.79,5.43,0.55,0.4,-0.14,-0.23,0.5,-2.33,-0.43,1.9,-1.94,-1.83,6.44,6.61,-1.75,0.62],
	[19.5,103.5,88469,388,-1,-37,-71,292.5,-2.9,0.3,-1.3,0,13.06,-3.62,-0.6,-0.39,0.1,-5,-0.2,0.4,0.2,0.4,-29.78,1154.42,1.27472,-0.00611,-0.00186,-0.00007,-0.0003,0.54775,-0.06607,-0.0193,-0.00209,-0.00493,3.0197,0.6868,0.2109,0.1217,0.0322,282.7,-0.1,0.7,-0.7,0.2,-8.28,-5.35,-6.65,-1.28,-0.8,5.28,1.33,0.34,-0.28,-0.39,1.94,-3.38,0.48,2.26,-1.4,-5.17,6.65,5.02,-3.28,-1.02],
	[19.5,104.5,93712,532,2,-25,-78,294.6,-4.5,-0.1,-1.1,-0.2,13.82,-3.81,-0.76,-0.54,0.09,-5.9,0.4,0.2,0,0.3,-28.05,659.43,1.27464,-0.00669,-0.00164,-0.00017,-0.00028,0.57502,-0.0608,-0.01621,-0.00398,-0.00491,2.6791,0.4871,0.1141,0.0633,0.0165,284.5,-1,0.6,-0.7,0.1,-8.08,-4.52,-6.75,-1.49,-0.88,5.89,2.33,0.59,-0.33,-0.39,3.83,-6.34,2.36,2.87,-0.91,-9,9,2.3,-3.25,-3.7],
	[19.5,105.5,99832,727,35,-14,-72,297.4,-5.4,-1.1,-0.6,-0.6,15.01,-4.17,-0.89,-0.63,0.02,-6.4,-0.3,0.8,-0.4,0.3,-25.6,109.01,1.27444,-0.00741,-0.00166,-0.00013,-0.0004,0.60619,-0.05611,-0.01499,-0.00483,-0.00464,2.4741,0.3244,0.0839,0.0321,0.0098,286.5,-1.9,0.3,-0.6,0,-8.13,-4.08,-6.36,-1.55,-0.68,5.32,0.86,1.59,0.16,-0.03,2.45,-11.7,4.55,4.22,0.35,-13.35,10.49,-0.86,-1.79,-2.37],
	[19.5,106.5,101035,773,69,-12,-66,297.9,-5.3,-2.1,-0.2,-0.7,15.95,-4.51,-0.88,-0.48,-0.08,-6.4,-0.4,2.4,-1,0.2,-23.07,1.38,1.27436,-0.00743,-0.00187,-0.00012,-0.00044,0.61066,-0.05481,-0.01535,-0.00544,-0.00416,2.6379,0.2881,0.1115,0.059,0.0078,287,-2,0,-0.5,-0.1,-8.33,-3.71,-6.1,-1.58,-0.62,5.62,-0.28,1.79,0.55,0.04,-0.6,-13.01,5.71,4.53,1.32,-7.39,7.55,1.4,-1.47,-1.38],
	[18.5,96.5,96843,358,-9,21,-19,298.6,-0.2,1.6,-1.5,-0.3,14.94,-4.35,-3.05,0.16,-1.09,-5.8,-1.2,-0.3,-0.5,0,-45.59,357.98,1.27605,-0.00322,0.00027,-0.00088,-0.00038,0.58978,-0.07666,-0.03863,-0.00002,-0.00053,2.969,0.9243,0.3851,0.1486,0.0559,287,1.6,1.7,-0.9,-0.4,-5.97,-3.76,-4.81,-0.85,-1.21,7.89,0.11,-0.1,0.25,-0.08,1.05,-5.09,-5.65,-0.03,2.44,14.7,2.12,10.78,-3.4,-3.81],
	[18.5,97.5,94165,343,-15,14,-32,297.6,-0.3,1.8,-1.6,-0.1,14.49,-4.3,-3.09,0.48,-0.96,-5.3,0.1,-0.2,-0.1,0,-42.33,604.14,1.27598,-0.00358,0.00008,-0.00085,-0.0003,0.57865,-0.07594,-0.03391,-0.00147,-0.00095,2.9404,0.8571,0.2005,0.2771,0.0177,285.9,1.5,1.6,-0.8,-0.3,-6.44,-3.69,-5.01,-0.91,-1.28,4.87,1.28,-0.1,-0.07,-0.43,2.02,-3.75,-4.62,0.13,0.53,9.93,6.12,16.15,-4.41,-3.51],
	[18.5,98.5,93698,362,-25,12,-46,297.4,-0.6,1.7,-1.6,0,14.6,-3.98,-2.52,0.43,-0.66,-4.8,0.2,0.4,-0.3,0.1,-39.91,646.86,1.27584,-0.00401,-0.00007,-0.00078,-0.00019,0.57672,-0.07413,-0.029,-0.00256,-0.00234,2.8967,0.8067,0.166,0.2841,0.0465,285.8,1.3,1.5,-0.8,-0.1,-6.91,-3.9,-5.2,-0.91,-1.21,4.8,1.48,0.19,-0.24,-0.48,1.98,-5.04,-4.81,1.26,-0.74,4.07,9.42,20.41,-4.41,-0.28],
	[18.5,99.5,95392,394,-28,14,-55,298.1,-0.7,1.4,-1.4,0.1,15.04,-3.57,-1.39,-0.1,-0.29,-5.2,-0.1,0,0,0.1,-37.83,490.38,1.27573,-0.00415,-0.00032,-0.00064,-0.00018,0.58473,-0.0728,-0.02735,-0.00178,-0.00364,2.7741,0.7349,0.2657,0.1536,0.0861,286.6,1.2,1.6,-0.8,0.1,-7.17,-3.65,-5.24,-0.95,-1.21,5.85,0.37,0.7,-0.22,-0.13,4.33,-6.4,-3.29,1.64,-1.67,-2.8,7.92,12.21,-1.28,0.96],
	[18.5,100.5,94894,398,-14,8,-56,297.4,-0.8,0.9,-1.2,0,14.92,-3.41,-0.7,-0.44,-0.09,-4.7,-0.1,-0.1,0.1,0,-35.42,537.64,1.27565,-0.00432,-0.00079,-0.00045,-0.00027,0.58224,-0.07217,-0.02597,-0.00074,-0.00387,2.7687,0.7014,0.3067,0.093,0.0817,286.2,1,1.4,-0.9,0.1,-7.09,-3.51,-5.17,-1.04,-1.2,6.02,0.81,0.88,-0.24,-0.14,5.32,-6.52,-1.07,0.82,-1.88,1.1,7.46,11.65,-0.68,1.65],
	[18.5,101.5,94829,409,-3,3,-58,297,-1.2,0.6,-1.1,-0.1,14.63,-3.45,-0.58,-0.48,-0.01,-5.1,-0.1,0,-0.1,-0.1,-33.24,546.59,1.27548,-0.00451,-0.00109,-0.00038,-0.00035,0.58232,-0.06902,-0.02235,-0.0011,-0.00342,2.6825,0.6063,0.2332,0.0762,0.0409,285.9,0.8,1.2,-0.8,0.1,-7.09,-3.7,-5.12,-1.13,-1.05,5.37,0.72,0.74,-0.33,-0.29,5.43,-8.89,-3.07,1.37,-2.14,-2.86,9.35,10.56,-1.3,2.13],
	[18.5,102.5,95299,427,6,2,-61,297.1,-1.7,0.4,-1.2,-0.1,14.32,-3.63,-0.63,-0.48,0.07,-5.6,0.2,0.1,-0.1,0.1,-31.6,506.14,1.27536,-0.00465,-0.00123,-0.00038,-0.00039,0.58494,-0.06548,-0.01968,-0.00183,-0.00273,2.5761,0.4998,0.1619,0.0621,0.0144,285.9,0.4,1,-0.8,0.1,-6.45,-3.5,-5.1,-1.23,-0.99,5.41,0.61,0.69,-0.12,-0.25,8.63,-9.52,-1.54,1.01,-2.08,-4.45,6.36,5.08,-1.78,0.66],
	[18.5,103.5,94919,444,11,0,-63,296.4,-2.4,0.5,-1.4,-0.2,14.06,-3.75,-0.73,-0.44,0.12,-5.8,0.1,-0.1,0.2,0.1,-29.85,543.47,1.27537,-0.0049,-0.00124,-0.00041,-0.00041,0.58188,-0.06311,-0.01847,-0.00276,-0.00273,2.5986,0.4744,0.1232,0.075,0.0084,285.5,0.1,0.9,-0.8,0,-5.36,-2.41,-5.38,-1.33,-1.16,5.74,1.28,0.62,-0.17,-0.3,9.13,-7.34,2.36,0.21,-1.56,-7.18,6.95,3.01,-3.24,-1.44],
	[18.5,104.5,93414,469,17,-9,-65,294.7,-3.6,0.2,-1.3,-0.3,13.96,-3.56,-0.75,-0.48,0.08,-5.4,0.4,0.3,0.1,0.5,-27.25,684.34,1.27541,-0.00555,-0.00129,-0.00028,-0.00036,0.57313,-0.06147,-0.01799,-0.0046,-0.00403,2.7876,0.5636,0.1328,0.1003,0.0192,284.6,-0.4,0.8,-0.7,0.1,-5.73,-1.62,-5.72,-1.5,-1.29,6.31,2.24,0.62,-0.23,-0.44,2.37,-3.53,3.72,1.08,-1.37,-3.46,7.73,1.07,-4.84,-4.11],
	[18.5,105.5,98221,617,34,-3,-66,296.9,-4.7,-0.4,-0.9,-0.4,14.89,-3.44,-0.85,-0.7,-0.03,-6.4,0.1,0.5,-0.1,0.1,-23.81,250.64,1.27524,-0.0064,-0.0012,-0.00022,-0.00034,0.59792,-0.05653,-0.01642,-0.00638,-0.0047,2.5683,0.4474,0.1104,0.047,0.0207,286.2,-1.2,0.5,-0.6,0,-7.02,-0.81,-5.55,-1.71,-1.26,6.33,2.33,1.02,-0.13,-0.41,-6.43,-3.77,2.83,3.29,0.14,-8.05,14.63,-0.08,-3.72,-4.07],
	[18.5,106.5,100967,708,66,-1,-60,298.5,-4.8,-1.4,-0.3,-0.6,15.99,-3.62,-0.77,-0.65,-0.14,-6.7,-0.4,1.9,-0.8,0.1,-21.14,6.47,1.2752,-0.00674,-0.00138,-0.00009,-0.00036,0.61082,-0.05351,-0.01616,-0.00717,-0.00452,2.6076,0.3678,0.1449,0.051,0.0268,287.3,-1.6,0.3,-0.5,-0.1,-8.34,-1.8,-5.1,-1.53,-0.97,6.12,-0.37,1.67,0.59,0.02,-9.36,-9.14,2.39,4.32,1.65,-7.16,7.89,0,-0.59,-1.47],
	[17.5,96.5,99787,346,-11,31,-7,299.7,0,1.4,-1.3,-0.5,15.93,-4.02,-2.38,-0.25,-0.99,-4.4,1.1,0.3,0.2,0.5,-45.58,94.3,1.27661,-0.00228,0.00065,-0.00083,-0.0003,0.6058,-0.07563,-0.03851,0.00116,-0.00184,2.9188,0.8974,0.491,0.0327,0.1255,288.2,1.9,1.9,-0.9,-0.3,-4.11,-2.35,-4.23,-0.67,-1.37,8.33,-0.65,-0.02,0.2,-0.05,-0.64,-9.15,-9.27,2.45,2.83,11.47,1.3,12.33,-3.31,-1.8],
	[17.5,97.5,96771,324,-11,22,-20,298.7,0.4,1.6,-1.3,-0.1,15.15,-4.14,-2.35,-0.07,-0.92,-5.6,-0.8,-0.3,-0.3,0,-42.3,366.33,1.2766,-0.00256,0.00036,-0.00075,-0.00025,0.59243,-0.07576,-0.03528,0.00059,-0.00143,2.8643,0.7869,0.3359,0.1077,0.0492,287,1.8,1.8,-0.8,-0.2,-4.86,-1.91,-4.3,-0.75,-1.49,6.14,0.76,-0.01,0.02,-0.38,2.66,-7.58,-5.84,1.14,1.07,11.07,6.38,14.88,-4.74,-1.59],
	[17.5,98.5,93918,335,-10,12,-35,297.2,-0.1,1.6,-1.3,0,14.65,-3.71,-1.82,-0.04,-0.6,-5.8,-1.3,-0.2,-0.4,0.1,-39.56,627.99,1.27656,-0.00337,-0.00007,-0.00059,-0.00023,0.57866,-0.07453,-0.03027,-0.00102,-0.00261,2.9172,0.808,0.272,0.1837,0.0486,285.8,1.6,1.6,-0.8,-0.1,-5.64,-1.96,-4.36,-0.82,-1.45,5.21,1.84,0.35,-0.3,-0.64,1.66,-6.73,-3.53,1.31,-0.19,4.96,11.41,22,-6.2,0.92],
	[17.5,99.5,97494,382,-22,23,-45,299.2,-0.4,1.3,-1.2,0,15.59,-3.28,-1.08,-0.29,-0.2,-5.2,-0.9,0,0,0,-37.35,300.53,1.2765,-0.00344,-0.00009,-0.00056,-0.00024,0.59622,-0.07094,-0.02695,-0.00169,-0.00363,2.707,0.6952,0.2648,0.1232,0.0778,287.4,1.4,1.7,-0.8,0.1,-6.14,-1.84,-4.47,-0.81,-1.4,5.91,0.29,0.65,-0.1,-0.14,1.13,-9.84,-4.06,2.15,-1.34,-5.38,9.45,8.6,-0.27,1.85],
	[17.5,100.5,95371,375,-5,13,-46,297.8,-0.7,0.9,-1.1,-0.1,15.1,-3.22,-0.48,-0.58,0.01,-5.3,-0.9,0.3,-0.6,0.2,-34.86,494.78,1.27648,-0.00383,-0.00056,-0.00038,-0.00036,0.58595,-0.07075,-0.02619,-0.00059,-0.00335,2.7863,0.6743,0.3336,0.0682,0.0879,286.4,1.2,1.5,-0.8,0.1,-6.08,-2.03,-4.33,-0.97,-1.29,6.19,0.78,0.88,-0.18,-0.17,0.18,-11.17,-1.49,2.02,-1.24,3.38,9.29,12.26,-1.52,1.63],
	[17.5,101.5,95447,397,3,9,-51,297.5,-1.3,0.6,-1,-0.1,14.85,-3.26,-0.37,-0.63,0.06,-5.2,-0.2,0.4,-0.5,0.2,-32.45,490.45,1.27634,-0.00413,-0.00081,-0.00032,-0.00039,0.58711,-0.06555,-0.02079,-0.0018,-0.00294,2.6854,0.5662,0.2491,0.0525,0.0393,286.2,0.8,1.3,-0.8,0,-5.8,-2.19,-4.29,-0.96,-1.18,5,1.17,0.69,-0.32,-0.52,-1.7,-12.14,-1.83,3.35,-1.04,-3.34,11.19,12.67,-2.48,2.42],
	[17.5,102.5,98747,453,5,19,-56,299.3,-1.9,0.6,-1.1,-0.1,14.93,-3.53,-0.54,-0.66,0.08,-6.3,-0.5,0.1,-0.2,0.2,-30.96,194.34,1.27623,-0.00417,-0.00077,-0.00039,-0.00041,0.60433,-0.06052,-0.01739,-0.00288,-0.00153,2.4335,0.403,0.1392,0.0286,0.0026,287.3,0.4,1.1,-0.8,0,-5.62,-2.14,-4.2,-0.94,-1.08,5.14,0.15,0.75,0.01,-0.19,-2.24,-13.72,-2.1,4.89,-0.38,-8.01,4.76,1.29,-1.28,-0.23],
	[17.5,103.5,98962,464,10,24,-54,299.1,-2.2,0.7,-1.3,-0.2,14.83,-3.71,-0.79,-0.58,0.1,-5.8,-0.1,0.2,-0.3,0.4,-29.22,175.46,1.27623,-0.00414,-0.00068,-0.00051,-0.00044,0.60457,-0.05951,-0.01775,-0.00356,-0.00127,2.4383,0.3802,0.094,0.0438,0.0019,287.2,0.2,1,-0.8,0,-5.95,-1.84,-4.18,-0.97,-1.07,5.95,0.42,0.89,-0.03,-0.14,-4.27,-12.61,-1.07,5.51,0.79,-3.97,4.79,1.52,-3.3,-2.08],
	[17.5,104.5,97962,471,21,21,-52,297.8,-2.8,0.6,-1.3,-0.3,14.78,-3.56,-0.83,-0.55,0.06,-6.2,0,0,-0.1,0.4,-26.22,267.43,1.27604,-0.00438,-0.00073,-0.00052,-0.00042,0.59866,-0.05904,-0.01778,-0.00533,-0.00209,2.527,0.4325,0.0967,0.0709,0.0117,286.7,0,0.9,-0.7,0,-5.61,-0.65,-4.28,-1.12,-1.22,7.24,1.54,0.94,-0.01,-0.26,-3.51,-10.55,-0.2,4.13,0.42,-0.3,5.43,-0.19,-4.9,-4.2],
	[17.5,105.5,96102,488,34,10,-54,296.1,-3.5,0.2,-1.1,-0.4,14.85,-3.17,-0.8,-0.6,-0.06,-5.5,0.3,0.5,-0.4,0.4,-22.25,436.92,1.27611,-0.00501,-0.00085,-0.00035,-0.00034,0.58857,-0.05778,-0.01795,-0.00733,-0.00388,2.7269,0.5366,0.1341,0.094,0.0326,285.8,-0.4,0.8,-0.6,0,-6.02,0.41,-4.4,-1.39,-1.41,7.79,2.67,0.73,-0.09,-0.54,-3.22,-5.72,2.05,2.2,-0.52,-1.81,12.39,-0.1,-4.49,-4.83],
	[17.5,106.5,99148,582,54,9,-53,297.8,-4.2,-0.4,-0.6,-0.5,15.46,-2.87,-0.73,-0.74,-0.18,-6.3,0.2,0.8,-0.6,0.2,-19.16,166.2,1.2761,-0.0057,-0.00091,-0.00021,-0.0003,0.60329,-0.05273,-0.017,-0.00854,-0.00478,2.5923,0.4676,0.1497,0.0538,0.0351,286.8,-1,0.6,-0.5,-0.1,-7.56,0.57,-4.17,-1.47,-1.39,7.41,1.47,1.28,0.28,-0.29,-9.02,-5.29,2.43,2.64,0.29,-5.93,12.73,0.13,-1.1,-2.11],
	[16.5,96.5,100783,330,3,26,2,300.3,-0.2,0.7,-0.9,-0.7,16.94,-2.85,-0.98,-0.73,-0.52,-5.2,2.1,0.5,0.4,1.1,-45.22,7.43,1.27713,-0.00189,0.00056,-0.0006,-0.00029,0.61332,-0.07448,-0.03789,0.00128,-0.00322,3.139,1.1744,0.816,-0.0402,0.3351,288.5,2,1.9,-0.8,-0.3,-3.93,-1.29,-3.5,-0.59,-1.33,8.39,-0.64,0.09,0.19,-0.06,-1.91,-12.87,-6.22,3.1,4.48,8.94,0.3,11.88,-2.57,-0.11],
	[16.5,97.5,99893,305,-4,27,-5,299.8,0.2,1.1,-1,-0.3,16.46,-3.28,-1.24,-0.59,-0.6,-4.9,1,0,0.5,0.2,-41.8,87.86,1.277,-0.00175,0.00047,-0.00058,-0.00026,0.60971,-0.07436,-0.03572,0.00084,-0.0031,2.9101,0.916,0.5447,0.0176,0.1769,288.1,2.1,1.8,-0.7,-0.2,-4.41,-1.15,-3.57,-0.6,-1.4,7.76,-0.3,0.1,-0.02,-0.22,-2.78,-13.48,-6.79,4.02,3.02,13,8.58,16.13,-6.14,0.68],
	[16.5,98.5,95436,306,-1,15,-23,297.7,0.3,1.2,-1,-0.1,15.14,-3.38,-1.2,-0.4,-0.48,-4.7,-1,0.2,-0.5,0.1,-38.59,490.04,1.27715,-0.00267,0.00003,-0.00047,-0.0003,0.58936,-0.07349,-0.03018,-0.00149,-0.00342,2.8593,0.7606,0.3448,0.1,0.0659,286.3,1.8,1.7,-0.7,-0.1,-5.15,-1.03,-3.58,-0.68,-1.38,5.98,2.12,0.49,-0.31,-0.76,-5.32,-12.12,-6.81,5.09,1.44,9.11,12.16,24.25,-7.94,3.15],
	[16.5,99.5,97942,354,-12,24,-36,299.4,-0.2,1.2,-1,0,15.61,-2.95,-0.82,-0.42,-0.11,-6.3,-0.7,-0.1,0,0.1,-36.26,262.11,1.27715,-0.00304,-0.00004,-0.00044,-0.0003,0.60186,-0.06823,-0.02509,-0.00282,-0.00336,2.6468,0.6531,0.2439,0.1153,0.0632,287.4,1.5,1.6,-0.7,0,-5.8,-1.02,-3.66,-0.69,-1.31,5.67,0.89,0.49,0.02,-0.39,-6.41,-14.13,-8.8,6.25,0.07,-7.78,12.64,7.81,0.1,2.23],
	[16.5,100.5,95808,353,1,17,-37,298.1,-0.6,0.9,-0.9,-0.1,15.18,-2.95,-0.37,-0.68,0.06,-6,-1.1,0.2,-0.4,0.3,-33.82,455.54,1.2772,-0.00348,-0.00038,-0.00031,-0.0004,0.59089,-0.06735,-0.02514,-0.00174,-0.00236,2.7561,0.6258,0.309,0.0571,0.0753,286.5,1.3,1.5,-0.7,0,-5.24,-1.1,-3.56,-0.77,-1.19,6.1,0.24,0.84,-0.1,-0.12,-1.61,-13.89,-5.79,5.38,0.41,1.86,9.43,8.56,-1.27,1.29],
	[16.5,101.5,97868,385,5,19,-42,298.9,-1.2,0.7,-0.9,-0.1,15.41,-3.04,-0.26,-0.79,0.07,-6.1,-0.1,0,-0.1,0.2,-31.42,271.77,1.27698,-0.00354,-0.00049,-0.00029,-0.00042,0.60218,-0.06197,-0.02035,-0.00261,-0.00151,2.5762,0.4844,0.2278,0.025,0.029,287.1,0.9,1.3,-0.7,0,-5.15,-0.35,-3.54,-0.85,-1.31,4.85,1.71,0.75,-0.36,-0.65,1.54,-13.32,-0.76,4.19,0.96,-3.14,11.17,12.06,-2.79,2.37],
	[16.5,102.5,98587,417,8,21,-48,299.3,-1.8,0.7,-0.9,-0.1,15.13,-3.16,-0.41,-0.78,0.05,-5.9,-0.1,0.3,-0.3,0.3,-29.93,207.82,1.27689,-0.00383,-0.00048,-0.00029,-0.00038,0.60547,-0.05842,-0.01745,-0.00388,-0.00075,2.4839,0.4113,0.1597,0.0231,0.0049,287.3,0.5,1.2,-0.7,0,-4.92,-0.17,-3.53,-0.91,-1.35,4.83,0.08,0.66,0.09,-0.23,-0.57,-12.83,0.97,3.92,1.58,-5.95,4.39,0.64,-1.95,-1.27],
	[16.5,103.5,98749,423,12,24,-46,299.2,-2.1,0.8,-1,-0.1,14.91,-3.21,-0.68,-0.67,0.01,-5.9,0.1,0.3,-0.3,0.5,-27.92,192.95,1.27695,-0.00377,-0.00037,-0.00036,-0.00037,0.60635,-0.05802,-0.018,-0.00484,-0.00081,2.4574,0.4045,0.1076,0.052,0.0007,287.3,0.4,1.1,-0.6,0,-4.65,-0.14,-3.47,-0.98,-1.34,6.21,0.16,0.9,-0.02,-0.12,-3.45,-11.62,0.87,4.05,1.58,-0.95,2.69,1.33,-3.37,-2.27],
	[16.5,104.5,98715,423,21,27,-43,298.8,-2.2,0.8,-1.1,-0.2,15,-3.27,-0.8,-0.59,-0.04,-5.8,0,0.3,-0.3,0.5,-25.1,197.36,1.27685,-0.00367,-0.00037,-0.00043,-0.00038,0.60509,-0.05671,-0.01764,-0.00614,-0.0015,2.4863,0.4093,0.0903,0.077,0.0054,287.1,0.3,1.1,-0.6,0,-5.18,-0.05,-3.42,-0.98,-1.32,7.27,0.36,0.95,0.13,-0.13,-5.5,-9.66,0.89,3.53,1.15,-1.79,3.99,-0.92,-2.99,-3.42],
	[16.5,105.5,98489,431,33,28,-40,298.2,-2.4,0.6,-1.1,-0.3,15.21,-3.18,-0.8,-0.59,-0.11,-6,-0.1,0.2,-0.3,0.3,-21.37,220.51,1.27674,-0.00379,-0.0004,-0.00046,-0.00037,0.60289,-0.05578,-0.01838,-0.00709,-0.00285,2.5407,0.427,0.1138,0.0858,0.0214,286.9,0.1,1,-0.6,0,-5.68,0.93,-3.36,-1.04,-1.46,8.95,1.51,0.8,0.13,-0.27,-7.21,-7.51,0.74,2.67,0.33,1.68,8.01,-1.8,-3.69,-4.35],
	[16.5,106.5,95336,428,49,11,-44,296.1,-3,0.3,-1,-0.3,14.86,-2.71,-0.72,-0.55,-0.22,-5.8,0,0.3,-0.2,0.4,-17.37,507.83,1.27703,-0.00441,-0.00061,-0.00029,-0.0003,0.58675,-0.05346,-0.01871,-0.00872,-0.00458,2.7548,0.5394,0.1612,0.1241,0.036,285.6,-0.2,0.9,-0.4,0,-6.81,0.6,-3.42,-1.07,-1.5,8.89,2.57,0.9,0.17,-0.51,-11.31,-4.53,0.92,2.7,-0.17,1.46,12.68,1.14,-2.67,-2.81],
	[15.5,96.5,100888,310,18,18,11,300.5,-0.4,0.1,-0.5,-0.8,17.34,-2.01,-0.15,-0.78,-0.23,-7.3,1.6,0.6,0.4,0.9,-44.67,0,1.27769,-0.0018,0.00038,-0.00043,-0.00031,0.617,-0.07242,-0.03659,-0.00044,-0.00446,3.1986,1.2626,0.9218,-0.0226,0.4054,288.3,2,1.7,-0.5,-0.2,-3.7,0,-2.81,-0.52,-1.26,8.51,-0.54,0.09,0.21,-0.07,-5.68,-17.7,-10.05,6.3,5.73,9.62,0.1,11.63,-3.34,0.9],
	[15.5,97.5,100096,280,11,20,6,299.8,0.2,0.4,-0.6,-0.5,16.93,-2.28,-0.39,-0.71,-0.34,-5.8,0.6,-0.2,0.7,0.3,-41,71.56,1.27749,-0.00161,0.00029,-0.00039,-0.00035,0.61438,-0.07221,-0.03469,-0.00107,-0.00469,2.9904,1.0521,0.6537,0.0363,0.2413,288,2.1,1.7,-0.5,-0.1,-4.14,0.01,-2.86,-0.52,-1.28,8.56,-0.65,0.29,-0.1,-0.08,-6.52,-16.86,-10.11,7.26,4.87,14.36,9.56,17.68,-7.44,2.17],
	[15.5,98.5,95582,277,8,13,-12,297.4,0.4,0.8,-0.7,-0.2,15.32,-2.99,-0.61,-0.65,-0.31,-4.7,-1,0.5,-0.6,0.2,-37.49,478.4,1.27767,-0.00235,0.00001,-0.00032,-0.00037,0.59463,-0.07026,-0.02875,-0.00335,-0.00431,2.8503,0.7326,0.4102,0.0583,0.1005,286.1,1.8,1.6,-0.5,0,-4.08,-0.4,-2.83,-0.63,-1.14,6.48,2.06,0.61,-0.29,-0.78,-5.27,-14.72,-6.82,6.34,4.24,10.23,12.61,23.43,-8.11,4.67],
	[15.5,99.5,98585,328,-5,24,-27,299.7,-0.2,1.1,-0.9,0,15.53,-2.7,-0.58,-0.55,-0.07,-7,-0.8,-0.2,0.1,0.2,-35.11,206.63,1.27767,-0.00274,-0.00001,-0.00031,-0.00034,0.60924,-0.06373,-0.02274,-0.00469,-0.0034,2.5365,0.5681,0.2173,0.0973,0.0564,287.5,1.5,1.6,-0.6,0.1,-4.34,-0.13,-2.92,-0.63,-1.17,5.57,1.24,0.39,0.12,-0.55,-4.13,-15.75,-5.69,6.23,4.14,-9.94,13.02,7.76,0.03,2.12],
	[15.5,100.5,99407,349,-3,28,-29,300.1,-0.8,0.9,-0.8,-0.1,15.75,-2.67,-0.37,-0.79,0.05,-6,-0.7,0.3,-0.4,0.4,-32.75,130.99,1.27773,-0.00297,-0.00011,-0.00026,-0.00041,0.61219,-0.06041,-0.02242,-0.00381,-0.00222,2.5272,0.506,0.2179,0.0444,0.053,287.8,1.2,1.5,-0.6,0,-4.47,0.02,-2.83,-0.7,-1.2,6.23,-0.29,0.87,-0.09,0.01,-3.62,-15.47,-4.59,6.25,4.24,-1.07,9.9,5.8,-0.61,1.64],
	[15.5,101.5,98458,351,8,23,-30,299.3,-1,0.7,-0.7,-0.1,15.56,-2.78,-0.12,-0.94,0.05,-5.4,0.1,0.5,-0.6,0.4,-30.32,217.87,1.27756,-0.0031,-0.00028,-0.00021,-0.00044,0.60723,-0.05799,-0.02048,-0.00378,-0.00134,2.5582,0.4477,0.2366,0.0071,0.0357,287.3,1,1.4,-0.6,0,-4.82,0.31,-2.74,-0.76,-1.29,5.34,1.19,0.77,-0.24,-0.52,-4.33,-14.62,-2.2,5.7,3.79,-2.25,9.47,7.96,-1.7,2.23],
	[15.5,102.5,99131,372,11,23,-36,299.7,-1.5,0.8,-0.8,0,15.39,-2.86,-0.3,-0.92,-0.03,-5.4,0.2,0.6,-0.6,0.4,-28.64,158.8,1.2774,-0.00327,-0.00023,-0.0002,-0.00037,0.61036,-0.0553,-0.0182,-0.00442,-0.00063,2.4776,0.388,0.1762,0.009,0.0076,287.5,0.7,1.3,-0.6,0,-5.64,0.55,-2.73,-0.78,-1.33,4.99,0.43,0.62,0.08,-0.32,-6.04,-13.62,-0.18,4.95,3.33,-2.99,3.95,0.69,-1.76,-1.12],
	[15.5,103.5,99439,376,14,26,-36,299.8,-1.6,0.9,-0.8,0,15.26,-2.92,-0.58,-0.81,-0.13,-5.9,0.3,0.4,-0.3,0.4,-26.58,130.95,1.27733,-0.0032,-0.00012,-0.00023,-0.00033,0.6126,-0.05529,-0.0185,-0.00511,-0.00065,2.443,0.3742,0.124,0.0348,-0.0066,287.6,0.6,1.2,-0.5,0,-6.14,0.43,-2.7,-0.76,-1.27,6.33,-0.02,0.86,0.05,-0.11,-7.43,-13.47,0.86,3.94,2.73,0.43,1.3,0.14,-2.48,-1.79],
	[15.5,104.5,99441,369,21,28,-34,299.7,-1.5,0.9,-0.9,0,15.33,-3.07,-0.74,-0.71,-0.22,-5.9,0.2,0.2,-0.2,0.3,-23.62,131.27,1.27735,-0.00299,-0.00009,-0.00029,-0.00032,0.61257,-0.05467,-0.01833,-0.00588,-0.00134,2.4527,0.36,0.1013,0.06,-0.0114,287.5,0.6,1.2,-0.5,0,-6.25,0.72,-2.61,-0.76,-1.29,7.88,0.14,0.85,0.08,-0.12,-8.64,-13.13,2.25,2.7,2.36,0.4,2.81,-1.2,-2.07,-3.22],
	[15.5,105.5,98236,351,33,26,-30,298.7,-1.3,0.9,-1,-0.1,15.21,-3.1,-0.79,-0.6,-0.28,-6,-0.1,0.1,-0.2,0.1,-19.83,241.12,1.27739,-0.00284,-0.00013,-0.00034,-0.00035,0.60586,-0.05374,-0.0188,-0.00669,-0.00278,2.5055,0.3659,0.1093,0.0882,-0.0027,287,0.6,1.1,-0.5,-0.1,-6.34,0.98,-2.59,-0.7,-1.32,8.87,0.12,0.81,0.15,-0.07,-9.86,-11.32,2.51,1.44,1.78,3.04,5.51,-1.67,-2.61,-3.39],
	[15.5,106.5,93010,320,52,8,-31,295.4,-1.7,0.8,-1,-0.2,14.47,-2.75,-0.73,-0.47,-0.3,-5.4,-0.3,0.1,-0.2,0.1,-15.32,721.58,1.27756,-0.00332,-0.00042,-0.00023,-0.00037,0.57877,-0.05282,-0.01982,-0.00786,-0.00451,2.8258,0.5074,0.1623,0.1604,0.0205,284.8,0.4,1,-0.4,0,-6.09,0.72,-2.65,-0.66,-1.33,8.99,1.78,0.79,0.11,-0.4,-7.83,-7.62,2.91,-0.12,0.92,3.3,8.98,0.41,-2.55,-2.62],
	[14.5,96.5,100902,284,24,14,18,300.7,-0.3,0,-0.3,-0.7,17.24,-1.88,-0.02,-0.7,-0.18,-8.5,0.8,0.2,0.3,0.5,-44.11,0,1.27811,-0.0017,0.00028,-0.00029,-0.00034,0.62106,-0.06823,-0.03443,-0.00249,-0.00629,3.0399,1.0733,0.7751,0.011,0.344,288.1,1.9,1.6,-0.4,-0.1,-3.19,0.24,-2.27,-0.45,-0.97,9,-0.61,0.05,0.19,-0.07,-9.8,-19.93,-12.72,7.98,6.84,10.12,0.08,11.45,-4.22,1.64],
	[14.5,97.5,100561,258,18,17,16,300.2,0.1,0.2,-0.4,-0.6,17.08,-1.84,-0.03,-0.69,-0.15,-7.8,0.4,0,0.3,0.2,-40.18,31.51,1.27784,-0.00149,0.00016,-0.00023,-0.00044,0.62119,-0.06786,-0.03251,-0.00316,-0.0064,2.8915,0.9716,0.5943,0.0666,0.2405,287.9,2,1.6,-0.4,-0.1,-3.53,0.27,-2.23,-0.49,-0.97,8.92,-0.74,0.34,-0.07,-0.03,-8.38,-19.1,-10.61,7.18,5.86,13.76,8.73,16.42,-6.95,3.47],
	[14.5,98.5,97178,255,12,15,-1,298.1,0.4,0.6,-0.6,-0.3,15.81,-2.62,-0.15,-0.8,-0.13,-5.2,-0.9,0.6,-0.5,0.5,-36.58,334.03,1.27794,-0.00193,0.00002,-0.0002,-0.00044,0.60686,-0.06468,-0.02741,-0.00491,-0.00553,2.7624,0.6722,0.4256,0.041,0.14,286.5,1.8,1.5,-0.4,0,-3.71,0.48,-2.17,-0.55,-1.02,6.86,2.02,0.6,-0.26,-0.81,-5.25,-16.16,-4.67,4.86,5.58,8.51,13.35,19.87,-6.82,5.21],
	[14.5,99.5,98841,296,2,23,-16,299.7,-0.2,0.9,-0.7,-0.1,15.56,-2.54,-0.32,-0.71,-0.05,-6.5,-0.1,0,0.1,0.1,-34.15,184.07,1.27804,-0.00242,0.00003,-0.00019,-0.00038,0.61405,-0.0592,-0.02281,-0.00585,-0.00438,2.5021,0.5053,0.2421,0.0674,0.0746,287.3,1.4,1.5,-0.5,0.1,-4.08,0.93,-2.18,-0.57,-1.17,5.5,1.15,0.35,0.22,-0.54,-4.91,-15.08,-0.36,4.16,6.12,-9.13,11.08,5.65,-0.13,1.75],
	[14.5,100.5,100088,316,4,27,-17,300.4,-0.8,0.8,-0.6,-0.1,15.91,-2.42,-0.17,-0.89,0.03,-5.5,0.9,-0.2,0.1,-0.1,-31.72,70.9,1.27805,-0.0026,-0.00002,-0.00016,-0.00043,0.61926,-0.05639,-0.02254,-0.00532,-0.0034,2.4888,0.4631,0.2373,0.0295,0.0697,287.9,1.2,1.5,-0.5,0.1,-4.41,1.17,-2.21,-0.6,-1.24,6.69,-0.38,0.84,-0.08,0.06,-7.03,-14.08,1.36,4.33,5.96,-0.06,9.62,4.81,-0.51,2.46],
	[14.5,101.5,97747,301,17,19,-17,298.8,-0.8,0.6,-0.5,-0.1,15.54,-2.48,0.03,-0.96,0.01,-5.5,0.1,0.6,-0.7,0.4,-28.85,281.53,1.27798,-0.00269,-0.00018,-0.00013,-0.00047,0.60736,-0.05539,-0.02132,-0.00543,-0.00251,2.5951,0.4442,0.2739,0.0119,0.0575,286.9,1.1,1.4,-0.4,0,-4.87,0.91,-2.12,-0.65,-1.17,5.97,0.78,0.72,-0.21,-0.41,-8.09,-11.46,2.72,4.09,5.15,-1.01,7.55,5.37,-1.18,2.35],
	[14.5,102.5,98477,317,18,21,-23,299.2,-1,0.7,-0.5,0,15.56,-2.61,-0.2,-0.95,-0.13,-5.3,0.2,0.6,-0.6,0.4,-26.39,216.41,1.27796,-0.00275,-0.00011,-0.00011,-0.00039,0.61027,-0.05334,-0.01963,-0.0055,-0.00158,2.5351,0.3818,0.2062,0.0145,0.0186,287.2,1,1.3,-0.4,0,-5.69,1.05,-2.07,-0.64,-1.2,5.83,0.56,0.57,0.09,-0.36,-7.93,-10.99,3.85,3.44,5.19,-2.75,3.68,-0.04,-0.69,-0.59],
	[14.5,103.5,99484,322,19,24,-25,299.9,-1,0.9,-0.6,0.1,15.61,-2.74,-0.5,-0.93,-0.29,-6,0.2,0.4,-0.4,0.2,-23.99,126.4,1.27793,-0.0026,0.00003,-0.0001,-0.0003,0.61585,-0.05262,-0.01989,-0.00547,-0.00146,2.4655,0.3366,0.1501,0.0205,-0.0072,287.6,0.9,1.3,-0.4,0,-6.44,1.24,-1.95,-0.64,-1.23,6.59,-0.03,0.74,0.03,-0.13,-7.23,-12.22,4.82,2.16,5.28,0.78,2.07,-1.18,-1.05,-1.49],
	[14.5,104.5,99471,313,24,25,-23,299.9,-0.7,1,-0.6,0.1,15.65,-2.92,-0.7,-0.85,-0.43,-6.1,-0.1,0.1,-0.3,0,-21.01,127.66,1.27794,-0.00235,0.0001,-0.00013,-0.00027,0.6162,-0.05185,-0.01993,-0.00559,-0.00186,2.4636,0.3072,0.1206,0.036,-0.0233,287.6,0.9,1.3,-0.4,0,-6.68,0.97,-1.86,-0.6,-1.15,7.98,0.13,0.78,0.04,-0.16,-7.03,-12.93,5.3,0.39,4.83,1.04,2.48,-1.73,-0.4,-2.64],
	[14.5,105.5,99112,302,32,26,-19,299.6,-0.4,1.1,-0.7,0,15.66,-3.1,-0.85,-0.75,-0.49,-5.9,-0.3,-0.1,-0.3,-0.1,-17.27,161,1.27795,-0.00209,0.00013,-0.00021,-0.0003,0.61389,-0.05023,-0.01976,-0.00582,-0.00273,2.4769,0.2852,0.1076,0.0525,-0.0197,287.4,0.8,1.2,-0.4,-0.1,-6.57,0.53,-1.76,-0.55,-1.05,8.88,-0.1,0.76,0.03,-0.05,-4.4,-12.25,5.94,-2.11,4.24,1.45,3.82,-2.29,-0.5,-3.03],
	[14.5,106.5,97020,296,46,23,-15,298.1,-0.5,1,-0.9,-0.2,15.22,-3.03,-0.92,-0.55,-0.39,-6,-0.7,-0.2,-0.2,0.1,-12.68,351.14,1.27791,-0.00222,0.00004,-0.00029,-0.00038,0.60162,-0.0483,-0.02042,-0.00627,-0.00401,2.5563,0.3086,0.1249,0.0891,0.0172,286.5,0.6,1.2,-0.4,-0.1,-6.31,0.81,-1.6,-0.48,-1.06,8.95,0.88,0.78,0.07,-0.19,0.51,-9.61,7.67,-4.72,4,0.21,6.59,-2.13,-0.49,-2.88],
	[13.5,96.5,100912,257,25,12,24,300.8,-0.2,0.1,-0.3,-0.6,17.14,-1.81,-0.02,-0.64,-0.19,-9.1,0.2,-0.1,0.2,0.2,-43.54,0,1.27835,-0.0016,0.00021,-0.00016,-0.00037,0.62561,-0.06331,-0.03216,-0.00384,-0.00788,2.8702,0.8667,0.6097,0.0376,0.2775,288,1.8,1.6,-0.3,0,-2.98,0.94,-1.64,-0.44,-0.84,9.93,-0.5,0.07,0.17,-0.1,-10.2,-22.41,-11.55,6.01,5.81,9.63,-0.47,9.62,-4.08,2.47],
	[13.5,97.5,100921,238,20,15,23,300.6,0.1,0.2,-0.3,-0.6,17.19,-1.62,0.09,-0.63,-0.08,-8.5,0.3,-0.1,0.2,0.2,-39.42,0,1.27819,-0.00141,0.00012,-0.00012,-0.00049,0.6272,-0.0625,-0.03029,-0.00444,-0.00758,2.7962,0.8507,0.5118,0.086,0.2208,287.9,1.9,1.5,-0.3,0,-3.52,0.67,-1.64,-0.47,-0.81,8.9,-0.66,0.27,0,-0.03,-7.35,-20.72,-7.78,4.61,5.22,10.12,6.9,12.38,-4.93,4.34],
	[13.5,98.5,98084,228,15,15,10,298.5,0.3,0.5,-0.5,-0.4,16.21,-2.23,0.03,-0.8,-0.06,-5.8,0.4,0.1,0.1,0.3,-35.9,252.98,1.27819,-0.00162,0.00008,-0.00013,-0.0005,0.61541,-0.05973,-0.02697,-0.00569,-0.00667,2.7367,0.6518,0.4201,0.0541,0.161,286.7,1.7,1.5,-0.3,0,-3.65,0.86,-1.6,-0.55,-0.9,7.47,1.35,0.57,-0.22,-0.63,-6.08,-16.32,-2.19,2.95,5.54,8.61,11.24,14.27,-5.42,4.54],
	[13.5,99.5,99039,258,10,21,-3,299.5,-0.1,0.7,-0.6,-0.2,15.94,-2.26,-0.08,-0.78,-0.04,-5.9,0.3,0.1,0.1,0.2,-33.05,167.27,1.27827,-0.00203,0.00008,-0.00012,-0.00044,0.61909,-0.05536,-0.02377,-0.00666,-0.00563,2.5463,0.5007,0.2959,0.0567,0.107,287.2,1.4,1.5,-0.3,0.1,-4.15,1.21,-1.65,-0.5,-1,5.57,1.35,0.38,0.2,-0.6,-7.68,-14.87,0.87,2.52,6.23,-7.8,9.99,3.94,-0.4,1.59],
	[13.5,100.5,100813,282,14,25,-4,300.7,-0.7,0.5,-0.4,-0.2,16.56,-2.15,0.08,-0.93,0.03,-6,0.9,-0.1,0.1,0.1,-30.6,7.85,1.27817,-0.00218,0.00003,-0.00007,-0.00046,0.62651,-0.05268,-0.02314,-0.00643,-0.00468,2.5368,0.4585,0.2959,0.0248,0.1056,288,1.3,1.4,-0.3,0.1,-4.13,1.15,-1.59,-0.49,-0.98,6.75,-0.47,0.73,0,0.09,-7.48,-14.47,2.17,2.11,6.88,0.33,8.19,3.36,0.04,2.67],
	[13.5,101.5,99852,270,23,21,-3,300,-0.6,0.4,-0.3,-0.2,16.51,-2.28,0.15,-1,0.01,-6.1,0.6,0.3,-0.1,0.1,-27.38,93.67,1.27803,-0.00211,-0.00009,-0.00004,-0.00052,0.62218,-0.05192,-0.02171,-0.00658,-0.00388,2.5741,0.4119,0.2931,0.0095,0.0874,287.6,1.2,1.3,-0.3,0,-4.24,1.34,-1.42,-0.6,-1.04,6.92,0.59,0.7,-0.2,-0.34,-8.22,-12.25,2.75,1.52,7.07,1.24,5.71,5.09,-1.05,2.39],
	[13.5,102.5,99791,275,23,21,-9,299.9,-0.5,0.6,-0.4,-0.1,16.38,-2.41,-0.09,-1.04,-0.19,-5.9,0.4,0.4,-0.2,0.2,-24.23,99.47,1.27811,-0.00213,-0.00001,0,-0.00041,0.62053,-0.04975,-0.01984,-0.00664,-0.00308,2.5442,0.3585,0.2346,0.0051,0.0434,287.6,1.1,1.3,-0.3,0,-5.03,1.42,-1.4,-0.52,-1.08,6.59,0.68,0.46,0.02,-0.44,-9.83,-11.75,2.11,1.23,6.98,-3.88,4.67,-0.25,0.1,-0.59],
	[13.5,103.5,100553,282,21,24,-13,300.3,-0.5,0.8,-0.3,0.1,16.36,-2.61,-0.36,-1.1,-0.44,-5.4,0.7,0.2,-0.2,0.1,-21.55,31.62,1.27813,-0.00205,0.00016,0.00004,-0.00027,0.6232,-0.04843,-0.02033,-0.00568,-0.00257,2.49,0.2983,0.1895,-0.0092,0.0059,287.9,1,1.3,-0.4,0,-5.39,1.2,-1.32,-0.45,-0.99,6.9,0.12,0.6,0.04,-0.19,-7.47,-12.66,2.85,-0.26,6.77,-1.97,2.62,-3.09,1.15,-1.3],
	[13.5,104.5,100136,276,26,23,-12,300.1,-0.2,0.9,-0.3,0.1,16.2,-2.84,-0.63,-1.04,-0.6,-5.2,0.6,0.1,-0.3,-0.1,-18.5,68.12,1.27815,-0.00188,0.00024,0.00002,-0.00022,0.62116,-0.04801,-0.02135,-0.00456,-0.00263,2.485,0.2553,0.1542,-0.0062,-0.0172,287.8,1,1.3,-0.4,0,-5.38,1.24,-1.19,-0.42,-1,7.98,0.14,0.7,-0.05,-0.17,-4.96,-12.23,4.5,-2.28,6.11,1.12,2.25,-2.18,1.4,-1.91],
	[13.5,105.5,100001,272,32,26,-8,299.9,-0.1,1,-0.5,0,16.24,-3.04,-0.85,-0.92,-0.62,-5.4,-0.1,-0.1,0,-0.1,-14.76,80.93,1.27819,-0.00164,0.0003,-0.00009,-0.00026,0.61986,-0.04616,-0.02143,-0.00384,-0.00321,2.5,0.2408,0.1321,0.0089,-0.0091,287.7,0.9,1.3,-0.4,-0.1,-5.4,1.41,-1.1,-0.34,-1.05,8.88,0.05,0.64,0,-0.05,-5.06,-11.18,4.72,-3.59,5.52,-1.08,2.26,-2.97,1.79,-2.74],
	[13.5,106.5,99519,276,40,30,-2,299.4,-0.2,1,-0.8,-0.3,15.96,-3.07,-0.99,-0.68,-0.48,-5.8,-0.8,-0.1,-0.3,0.1,-10.24,125.7,1.2782,-0.00161,0.00035,-0.00028,-0.00037,0.61549,-0.04335,-0.02171,-0.00347,-0.00446,2.495,0.2473,0.1366,0.0364,0.0331,287.5,0.7,1.3,-0.5,-0.1,-5.4,1.33,-0.96,-0.31,-0.98,9.21,0.37,0.63,0.03,-0.08,-5.38,-10.35,4.39,-5.17,5.65,-1.87,5.13,-3.74,1.95,-3.54],
	[12.5,96.5,100919,229,22,11,29,300.9,-0.1,0.2,-0.2,-0.5,17.17,-1.61,-0.03,-0.55,-0.21,-9.2,0,-0.2,0,0.1,-42.75,0,1.27843,-0.00152,0.00022,-0.0001,-0.00041,0.62953,-0.05739,-0.02998,-0.00475,-0.00897,2.7611,0.7159,0.4975,0.058,0.2373,287.8,1.6,1.5,-0.2,0,-2.89,1.18,-1.05,-0.38,-0.65,10.91,-0.41,0.09,0.13,-0.09,-8.96,-24.72,-8.5,3.67,5.47,8.68,-1.18,6.86,-3.32,3.06],
	[12.5,97.5,100922,217,19,14,28,300.8,0,0.2,-0.3,-0.5,17.28,-1.43,0.1,-0.56,-0.12,-8.8,0.3,-0.1,0.2,0.1,-38.13,0,1.2785,-0.00137,0.00017,-0.00007,-0.00051,0.63071,-0.0565,-0.02857,-0.00497,-0.00836,2.7418,0.7374,0.4527,0.0911,0.2047,287.8,1.7,1.5,-0.2,0,-3.59,1.13,-1.01,-0.38,-0.63,9.04,-0.42,0.2,0.06,-0.08,-7.68,-23.35,-5.47,2.38,5.83,6.88,4.36,7.98,-2.87,3.96],
	[12.5,98.5,99523,206,16,16,21,299.4,0.2,0.4,-0.4,-0.4,16.76,-1.79,0.04,-0.68,-0.11,-6.2,0,0,-0.1,0.3,-34.72,124.9,1.27831,-0.00135,0.00018,-0.00008,-0.00054,0.62564,-0.05448,-0.02644,-0.00574,-0.00763,2.6841,0.6161,0.3837,0.0745,0.1601,287.1,1.6,1.5,-0.2,0.1,-3.79,1.24,-1,-0.43,-0.69,7.86,0.7,0.53,-0.17,-0.43,-5.24,-20.09,-1.32,0.57,6.74,9.14,8.63,9.86,-4.21,3.62],
	[12.5,99.5,97987,222,19,17,9,298.6,0,0.6,-0.5,-0.3,16.15,-1.97,0,-0.75,-0.07,-6.2,0,0.1,-0.1,0.2,-31.5,261.26,1.27844,-0.00178,0.00014,-0.00009,-0.00048,0.61747,-0.05167,-0.0244,-0.00653,-0.00694,2.6701,0.52,0.3529,0.0553,0.1516,286.5,1.4,1.4,-0.3,0.1,-3.89,1.43,-1.04,-0.42,-0.8,5.88,1.29,0.47,0.1,-0.56,-4.56,-16.83,1.81,-0.77,7.18,-6.46,8.16,2.5,0.07,1.43],
	[12.5,100.5,100877,250,23,21,9,300.8,-0.6,0.3,-0.4,-0.3,17.12,-1.63,0.2,-0.77,0,-7.4,0.6,0.1,0.1,0.1,-29.05,3.2,1.27826,-0.00191,0.00007,-0.00005,-0.00049,0.62973,-0.04892,-0.02335,-0.00634,-0.00617,2.6273,0.5164,0.3439,0.0552,0.144,287.8,1.2,1.4,-0.2,0.1,-4.34,1.66,-1.11,-0.35,-0.91,6.42,-0.29,0.53,0.07,0.02,-6,-15.65,3.04,-1.74,7.88,-0.19,6.11,2.31,0.77,2.29],
	[12.5,101.5,100420,238,30,19,12,300.5,-0.5,0.2,-0.2,-0.3,17.14,-1.67,0.17,-0.79,-0.03,-7.4,0.8,0,0.2,0.1,-26,43.45,1.27808,-0.0018,-0.00001,-0.00002,-0.00056,0.62826,-0.04925,-0.02266,-0.00596,-0.00554,2.6371,0.4854,0.3202,0.0485,0.1225,287.7,1.2,1.3,-0.2,0,-4.28,1.69,-1.02,-0.38,-0.91,7.39,0.07,0.64,-0.16,-0.16,-6.25,-13.37,4.11,-2.66,7.61,3.16,3.78,4.38,-0.19,2.51],
	[12.5,102.5,98576,226,31,16,6,299.2,-0.1,0.4,-0.2,-0.2,16.47,-1.98,0.01,-0.9,-0.17,-6.4,0.1,0.4,-0.5,0.1,-22.43,208.82,1.27818,-0.00181,0.00001,0.00001,-0.0005,0.61921,-0.04822,-0.02115,-0.00598,-0.00501,2.631,0.4044,0.2855,0.0206,0.0908,286.9,1.2,1.3,-0.3,0,-4.74,2.08,-0.92,-0.37,-1.01,7.49,0.45,0.45,-0.14,-0.36,-10.41,-11.41,3.62,-2.73,6.66,0.52,4.86,1.59,-0.07,0.14],
	[12.5,103.5,98902,234,30,17,-1,299.4,0,0.7,-0.2,0,16.25,-2.23,-0.25,-0.95,-0.42,-5.5,-0.4,0.2,-0.5,-0.1,-19.42,180.36,1.27815,-0.00181,0.00013,0.00006,-0.00034,0.61909,-0.04561,-0.02078,-0.00526,-0.00444,2.553,0.3208,0.2291,0.0045,0.0453,287.1,1.1,1.3,-0.3,0,-5.34,2.11,-0.86,-0.24,-0.99,6.67,0.62,0.48,0.02,-0.35,-14.64,-11.62,1.75,-2.12,6.51,-4.41,3.12,-3.48,2.26,-1.12],
	[12.5,104.5,100229,245,30,20,-3,300.1,0,0.8,-0.2,0.1,16.49,-2.57,-0.53,-1,-0.65,-5,0.2,0,-0.2,0,-16.42,60.86,1.2781,-0.00164,0.00027,0.00009,-0.00022,0.62418,-0.04386,-0.02208,-0.00359,-0.00426,2.497,0.2461,0.1843,-0.0127,0.01,287.7,0.9,1.3,-0.4,0,-5.13,1.78,-0.78,-0.21,-0.92,7.76,0.02,0.6,-0.07,-0.07,-12.42,-13.23,2.02,-2.45,6.92,-1.5,2.66,-3.49,2.5,-1.54],
	[12.5,105.5,100358,242,34,22,2,299.9,0,0.9,-0.4,-0.1,16.59,-2.84,-0.77,-0.93,-0.67,-4.7,0.1,-0.2,-0.3,-0.1,-12.67,49.03,1.27813,-0.00139,0.00035,0.00001,-0.00027,0.62433,-0.0424,-0.0229,-0.00238,-0.00483,2.5021,0.2083,0.1579,-0.0138,0.0122,287.7,0.8,1.3,-0.4,0,-4.77,1.92,-0.67,-0.18,-0.98,8.78,-0.18,0.52,-0.09,-0.04,-11.24,-15.12,2.15,-2.77,6.36,0.53,2.55,-2.06,2.35,-1.98],
	[12.5,106.5,99161,233,40,22,6,299.2,-0.1,1,-0.7,-0.3,16.07,-2.82,-0.9,-0.68,-0.53,-5.7,-0.7,-0.1,-0.2,0.1,-8.22,157.46,1.27826,-0.0013,0.00041,-0.00017,-0.00037,0.6178,-0.04012,-0.02312,-0.00165,-0.00601,2.501,0.2236,0.1617,0.0135,0.0517,287.2,0.7,1.3,-0.5,-0.1,-4.63,2.47,-0.46,-0.12,-1.03,9.21,-0.07,0.4,-0.04,-0.07,-14.62,-15.17,0.79,-2.5,5.03,0.92,4.72,-2.84,2.72,-3.27],
	[11.5,96.5,100923,201,17,10,32,300.9,-0.1,0.3,-0.2,-0.4,17.24,-1.38,-0.02,-0.46,-0.25,-9.2,0,-0.2,0,0,-41.47,0,1.27847,-0.00147,0.00027,-0.00009,-0.00042,0.63289,-0.05063,-0.02799,-0.00533,-0.0098,2.6919,0.5934,0.4292,0.0624,0.2193,287.7,1.5,1.5,-0.1,0.1,-3.52,1.65,-0.58,-0.27,-0.49,11.48,-0.3,0.07,0.12,-0.11,-7.86,-25.34,-7.42,1.69,5.9,7.35,-1.37,4.38,-2.27,2.86],
	[11.5,97.5,100924,194,15,14,31,300.8,-0.1,0.3,-0.2,-0.4,17.33,-1.27,0.06,-0.5,-0.2,-9,0.2,-0.1,0.1,0,-36.3,0,1.27868,-0.00135,0.00025,-0.00006,-0.00048,0.63412,-0.04989,-0.02713,-0.00522,-0.0093,2.6827,0.6101,0.4053,0.0785,0.1958,287.6,1.5,1.4,-0.2,0.1,-3.82,1.4,-0.58,-0.26,-0.44,9.45,-0.23,0.17,0.07,-0.13,-6.72,-25.13,-4.94,0.45,6.59,5.72,2.94,4.77,-1.28,3.04],
	[11.5,98.5,100056,186,13,16,27,299.8,0,0.5,-0.3,-0.3,17.08,-1.47,-0.03,-0.58,-0.23,-7.1,-0.5,-0.2,0.4,0,-32.7,77.21,1.27847,-0.0013,0.00029,-0.00008,-0.0005,0.63104,-0.04837,-0.02573,-0.00527,-0.00885,2.6589,0.543,0.3565,0.0657,0.1626,287.2,1.4,1.4,-0.2,0.1,-4.13,1.26,-0.56,-0.33,-0.45,8.08,0.44,0.48,-0.08,-0.32,-5.33,-23.68,-1.69,-1.5,7.43,6.98,6.65,5.87,-1.95,2.45],
	[11.5,99.5,99565,200,18,18,21,299.4,-0.1,0.5,-0.4,-0.3,16.98,-1.48,0,-0.61,-0.18,-6.4,-0.2,0,-0.2,0.2,-29.56,120.55,1.27846,-0.00154,0.00027,-0.0001,-0.00047,0.62797,-0.04623,-0.02449,-0.00531,-0.00828,2.6724,0.5014,0.3517,0.0535,0.1642,287,1.2,1.4,-0.2,0.1,-4.34,1.52,-0.51,-0.41,-0.56,6.09,1.16,0.54,0.08,-0.45,-4.96,-22.06,1.02,-3.64,7.99,-5.47,6.68,0.52,1.59,1.19],
	[11.5,100.5,100921,218,26,19,21,300.8,-0.4,0.3,-0.3,-0.4,17.45,-1.22,0.15,-0.57,-0.1,-8.5,0.4,0,0.1,0.1,-27.05,0,1.27829,-0.00173,0.00019,-0.0001,-0.00049,0.63349,-0.04502,-0.02404,-0.00474,-0.0076,2.6518,0.4991,0.3583,0.0535,0.1656,287.7,1.1,1.3,-0.2,0.1,-4.47,1.84,-0.46,-0.4,-0.69,6.4,-0.16,0.44,0.07,-0.03,-6.83,-21.35,2.32,-4.95,8.27,1.14,4.07,2.01,1.15,2.41],
	[11.5,101.5,100920,213,32,17,23,301,-0.4,0.2,-0.2,-0.4,17.37,-1.26,0.1,-0.57,-0.12,-8.9,0.4,0,0.1,0,-24.19,0,1.27811,-0.00166,0.00011,-0.00006,-0.00056,0.6339,-0.04574,-0.02394,-0.00389,-0.00682,2.609,0.464,0.3262,0.0474,0.1414,287.7,1.1,1.3,-0.2,0,-4.26,1.64,-0.36,-0.44,-0.64,7.21,0.03,0.55,-0.09,-0.13,-7.01,-19.31,3.63,-5.58,8.33,2.59,1.54,2.17,1.13,2.53],
	[11.5,102.5,100521,199,34,15,22,300.4,-0.1,0.2,-0.1,-0.3,17.36,-1.35,0.04,-0.64,-0.19,-8.1,0.3,-0.1,0.2,0,-20.55,36.32,1.27806,-0.00148,0.00009,0,-0.00056,0.63211,-0.04519,-0.02329,-0.00347,-0.00635,2.6068,0.4418,0.2956,0.0413,0.1166,287.5,1.2,1.3,-0.3,0,-3.99,1.66,-0.32,-0.51,-0.72,8.11,-0.21,0.41,-0.19,-0.17,-4.78,-16.81,5.5,-6.66,7.99,3.96,3.55,2.29,-0.08,0.94],
	[11.5,103.5,97761,194,36,13,13,298.4,0.1,0.5,-0.1,-0.2,16.42,-1.95,-0.18,-0.8,-0.35,-5.7,-0.1,0,-0.4,-0.1,-17.14,282.16,1.27821,-0.00159,0.00013,0,-0.00045,0.6181,-0.04354,-0.0233,-0.00299,-0.00654,2.6403,0.3248,0.2773,-0.0047,0.0972,286.4,1,1.3,-0.3,0,-4.13,1.93,-0.32,-0.44,-0.82,7.03,0.72,0.47,-0.13,-0.36,-6.43,-14.52,4.66,-6.07,6.72,-1.36,3.41,-1.61,1.46,-0.84],
	[11.5,104.5,99999,212,34,17,10,299.7,0,0.8,-0.3,0,16.75,-2.21,-0.51,-0.79,-0.59,-5.4,0,-0.1,0.1,0,-14.16,82.06,1.27805,-0.00148,0.00028,0.00004,-0.00029,0.62707,-0.03996,-0.02357,-0.00208,-0.00662,2.5125,0.2406,0.202,-0.0066,0.049,287.4,0.9,1.3,-0.4,0,-4.59,2.12,-0.32,-0.31,-0.92,7.2,0.29,0.47,-0.06,-0.12,-10.46,-16.77,3.05,-5.98,6.62,-3.57,5.08,-3.71,2.86,-1.69],
	[11.5,105.5,100766,214,35,18,13,300,0,0.9,-0.4,-0.1,16.8,-2.42,-0.76,-0.74,-0.66,-4.8,0.1,-0.4,0,0,-10.51,13.3,1.27805,-0.0013,0.00041,-0.00001,-0.00029,0.62999,-0.03752,-0.02419,-0.00091,-0.00698,2.4617,0.1842,0.1689,-0.0157,0.0403,287.7,0.7,1.3,-0.4,0,-4.35,1.9,-0.16,-0.36,-0.85,8.65,-0.35,0.37,-0.09,-0.01,-10.88,-19.83,3.17,-7.29,6.67,1.72,4.71,-1.95,3.35,-1.22],
	[11.5,106.5,100020,198,39,17,16,299.7,0.1,0.9,-0.5,-0.3,16.49,-2.47,-0.85,-0.62,-0.58,-5.7,-0.2,-0.4,0.2,0,-6.05,80.37,1.27826,-0.00104,0.00051,-0.00011,-0.00038,0.62628,-0.03538,-0.02449,0.00027,-0.00752,2.4588,0.1798,0.1792,-0.0182,0.061,287.4,0.7,1.3,-0.5,0,-4.02,1.64,-0.08,-0.39,-0.75,9.01,-0.43,0.25,-0.09,-0.08,-8.14,-20.28,4.54,-9.29,5.74,0.87,4.22,-3.15,3.93,-1.97],
	[10.5,96.5,100925,172,10,9,34,301,-0.2,0.4,-0.2,-0.3,17.3,-1.14,-0.01,-0.38,-0.3,-9.3,-0.1,-0.1,0,0,-39.66,0,1.27856,-0.00142,0.00034,-0.0001,-0.00038,0.63588,-0.04381,-0.02583,-0.00533,-0.01071,2.6409,0.486,0.3799,0.0558,0.213,287.5,1.3,1.4,-0.1,0.2,-3.95,1.97,-0.16,-0.24,-0.37,11.69,-0.23,0.03,0.13,-0.13,-7.49,-26.24,-7.22,-0.7,5.74,5.54,-1.03,2.74,-1.14,2.37],
	[10.5,97.5,100926,168,9,12,33,300.9,-0.2,0.4,-0.2,-0.3,17.4,-1.05,0.04,-0.4,-0.28,-9.2,0.1,-0.1,0,0,-34.36,0,1.27877,-0.00134,0.00035,-0.00009,-0.00041,0.63696,-0.04318,-0.02549,-0.00496,-0.01036,2.6409,0.4989,0.3694,0.0602,0.1968,287.5,1.2,1.4,-0.1,0.2,-4.24,2.05,-0.07,-0.25,-0.38,9.82,-0.15,0.11,0.06,-0.14,-6.14,-26.59,-5.74,-1.71,6.25,4.91,2.88,2.36,-0.04,2.24],
	[10.5,98.5,100004,165,9,15,29,299.8,-0.2,0.5,-0.3,-0.2,17.22,-1.13,-0.08,-0.46,-0.35,-7.2,-0.4,-0.2,0,-0.1,-30.54,81.63,1.27857,-0.00136,0.00041,-0.00012,-0.00042,0.63316,-0.04175,-0.02486,-0.00439,-0.01017,2.6525,0.4653,0.3397,0.0469,0.171,287,1.2,1.4,-0.2,0.1,-4.25,1.89,-0.04,-0.27,-0.39,8.1,0.56,0.45,0.01,-0.29,-4.61,-25.91,-3.38,-3.08,6.78,5.09,5.75,2.69,0.38,1.52],
	[10.5,99.5,100657,180,14,18,28,300.2,-0.3,0.5,-0.3,-0.2,17.51,-1.06,-0.02,-0.49,-0.31,-7.6,0.4,-0.1,0,0,-27.44,24.01,1.27837,-0.00154,0.00041,-0.00017,-0.00041,0.6354,-0.03972,-0.0243,-0.0037,-0.00975,2.6411,0.4299,0.3418,0.0255,0.1671,287.3,1.1,1.3,-0.2,0.1,-4.33,1.84,-0.01,-0.33,-0.44,6.29,0.92,0.54,0.14,-0.31,-4.37,-25.54,-0.97,-4.79,7.37,-4.45,5.76,-1.57,3.21,0.92],
	[10.5,100.5,100926,191,23,18,29,300.8,-0.4,0.3,-0.2,-0.3,17.52,-1,0.1,-0.47,-0.22,-9,0.2,0,0,0,-24.85,0,1.27823,-0.00169,0.00033,-0.00017,-0.00045,0.63641,-0.03936,-0.02453,-0.0029,-0.00908,2.6253,0.4104,0.364,0.0108,0.1759,287.5,1,1.3,-0.2,0.1,-4.56,1.89,0,-0.35,-0.49,6.38,-0.11,0.42,0.04,-0.03,-4.45,-24.24,0.94,-5.72,7.69,1.6,2.35,0.97,1.55,2.4],
	[10.5,101.5,100925,189,29,16,30,301.1,-0.4,0.2,-0.1,-0.3,17.29,-1.11,0.02,-0.48,-0.24,-9.3,0.1,-0.1,0,0,-21.95,0,1.27808,-0.00163,0.00027,-0.00014,-0.0005,0.63655,-0.03977,-0.02447,-0.00215,-0.00844,2.5621,0.3723,0.3327,0.0048,0.154,287.5,1,1.3,-0.2,0,-4.47,1.97,0.02,-0.32,-0.53,6.92,0.11,0.45,-0.06,-0.13,-3.49,-22.62,2.52,-6.12,7.9,1.4,-0.47,0.26,1.81,2.23],
	[10.5,102.5,100928,184,32,16,31,301,-0.3,0.2,-0.1,-0.4,17.33,-1.13,-0.02,-0.5,-0.24,-8.9,0.2,-0.1,0,0,-18.41,0,1.27803,-0.00149,0.00023,-0.00009,-0.00053,0.63572,-0.03925,-0.02437,-0.0011,-0.00808,2.5556,0.3763,0.3169,0.0045,0.1416,287.6,0.9,1.3,-0.3,0,-4.25,1.88,0.11,-0.35,-0.55,7.9,-0.25,0.37,-0.13,-0.07,-2.49,-20.59,4.16,-7.05,7.89,1.93,1.17,0.63,1.35,1.46],
	[10.5,103.5,100518,179,34,16,28,300.5,-0.2,0.4,-0.1,-0.3,17.32,-1.31,-0.15,-0.55,-0.33,-8,0.3,-0.1,0.1,-0.1,-14.81,36.54,1.27804,-0.00136,0.00025,-0.00006,-0.00048,0.6334,-0.03804,-0.02458,-0.00019,-0.00814,2.5735,0.3544,0.2929,0.0005,0.1253,287.4,0.9,1.3,-0.3,0,-4.16,1.92,0.2,-0.4,-0.62,7.64,0.36,0.43,-0.21,-0.23,-1.9,-18.88,4.88,-8.27,7.14,1.97,2.83,-0.48,1.35,-0.08],
	[10.5,104.5,100446,185,33,17,23,300.1,-0.2,0.7,-0.3,-0.2,17.08,-1.66,-0.46,-0.55,-0.51,-6.4,0.5,-0.3,0.2,-0.1,-11.64,42.5,1.27798,-0.00136,0.00037,-0.00009,-0.00037,0.63236,-0.03543,-0.02483,0.00072,-0.00828,2.5228,0.2657,0.2386,-0.0122,0.0937,287.4,0.7,1.3,-0.4,0,-4.5,2.26,0.27,-0.38,-0.75,7.27,0.4,0.35,-0.11,-0.2,-4.02,-18.32,3.83,-8.9,5.91,-0.89,6.31,-2.32,3.05,-1.2],
	[10.5,105.5,100874,189,35,17,24,300.1,-0.2,0.8,-0.4,-0.2,16.96,-1.85,-0.69,-0.5,-0.61,-5.5,0.2,-0.6,0.3,-0.1,-8.04,4.52,1.27798,-0.00129,0.00051,-0.00015,-0.00034,0.63352,-0.03236,-0.02489,0.00198,-0.00848,2.4673,0.2053,0.2138,-0.0285,0.0794,287.5,0.6,1.3,-0.5,0,-4.63,2.24,0.34,-0.35,-0.71,8.3,-0.12,0.26,-0.03,-0.09,-4.66,-17.28,4.01,-9.7,5.29,-0.13,4.94,-3.57,5.05,-1.04],
	[10.5,106.5,100807,186,43,15,28,300.1,-0.3,0.6,-0.4,-0.5,17.1,-1.77,-0.58,-0.49,-0.54,-6.4,0.2,-0.5,0.3,0.1,-3.91,10.77,1.27813,-0.00111,0.00056,-0.00019,-0.0004,0.632,-0.02939,-0.0253,0.00373,-0.0085,2.5327,0.2411,0.2848,-0.058,0.1076,287.5,0.5,1.2,-0.5,-0.1,-3.97,1.89,0.37,-0.36,-0.66,8.47,-0.14,0.24,-0.11,-0.11,-0.52,-15.6,6.38,-11.25,4.84,-4.3,3.33,-5.47,6.18,-0.87],
	[9.5,96.5,100924,144,3,8,35,301,-0.2,0.4,-0.2,-0.2,17.33,-0.92,-0.01,-0.29,-0.34,-9.3,0,-0.1,0,0,-37.38,0,1.2787,-0.00142,0.0004,-0.00011,-0.00033,0.63891,-0.03672,-0.02348,-0.00471,-0.0114,2.5898,0.3758,0.3288,0.0462,0.2066,287.4,1,1.3,-0.1,0.2,-4.29,1.87,0.01,-0.13,-0.22,11.54,-0.22,-0.05,0.13,-0.11,-8.72,-28.45,-8.7,-1.3,4.2,3.34,-0.48,2.16,-0.56,2.21],
	[9.5,97.5,100925,139,0,11,34,300.9,-0.2,0.5,-0.2,-0.2,17.42,-0.87,-0.04,-0.3,-0.37,-9,0,-0.1,0,0,-32.29,0,1.27883,-0.00133,0.00045,-0.00012,-0.00034,0.63966,-0.03601,-0.02313,-0.00433,-0.01107,2.5915,0.3805,0.3109,0.0493,0.1856,287.4,1,1.3,-0.1,0.2,-4.21,1.86,0.06,-0.14,-0.25,10.12,-0.24,0.04,0.07,-0.12,-7.78,-28.57,-8.83,-1.28,4.03,4.75,2.96,1.6,0.42,2.02],
	[9.5,98.5,99593,136,1,14,30,299.5,-0.2,0.7,-0.3,-0.1,17.12,-0.91,-0.19,-0.32,-0.47,-6.9,0.2,-0.3,-0.2,-0.2,-28.14,118.1,1.2787,-0.00137,0.00051,-0.00018,-0.00035,0.63376,-0.03452,-0.02279,-0.00357,-0.01094,2.6257,0.3645,0.2906,0.0395,0.1616,286.8,1,1.3,-0.2,0.2,-4.57,1.94,0.12,-0.17,-0.31,8.2,0.54,0.37,0.03,-0.25,-6.95,-26.68,-7.71,-1.59,3.96,3.12,4.61,0.28,1.85,0.97],
	[9.5,99.5,100496,154,8,17,31,299.9,-0.4,0.5,-0.3,-0.1,17.44,-0.82,-0.11,-0.39,-0.48,-7.1,0.2,-0.1,-0.1,-0.1,-24.73,38.34,1.27843,-0.00158,0.00051,-0.00022,-0.00034,0.6367,-0.03258,-0.02285,-0.00229,-0.01062,2.6044,0.3374,0.3081,0.0014,0.1525,287.2,0.8,1.3,-0.2,0.2,-5.05,2.12,0.22,-0.19,-0.38,6.62,0.83,0.51,0.18,-0.23,-6.24,-25.55,-5.64,-3.04,4.65,-3.79,4.4,-3.77,4.64,0.36],
	[9.5,100.5,100929,168,19,17,34,300.8,-0.6,0.3,-0.1,-0.3,17.49,-0.76,0.04,-0.42,-0.39,-8.8,0.1,-0.1,-0.1,0,-22.11,0,1.27826,-0.00173,0.00042,-0.00022,-0.00039,0.6378,-0.03241,-0.02396,-0.00098,-0.01015,2.5942,0.327,0.3504,-0.0251,0.1671,287.4,0.7,1.2,-0.3,0.1,-5.05,2.03,0.26,-0.23,-0.41,6.13,0.03,0.4,0.06,-0.03,-4.16,-24.55,-2.78,-4.8,5.57,-1.25,0.99,-1.68,2.7,1.8],
	[9.5,101.5,100928,168,25,16,36,301.1,-0.6,0.2,-0.1,-0.3,17.26,-0.86,0,-0.4,-0.35,-9.3,0,-0.1,0,-0.1,-19.36,0,1.27812,-0.00165,0.00039,-0.00021,-0.00045,0.63746,-0.03289,-0.02452,-0.00015,-0.00976,2.5514,0.3131,0.3473,-0.0318,0.1637,287.4,0.7,1.2,-0.3,0.1,-4.85,1.94,0.35,-0.28,-0.43,6.8,0,0.35,-0.05,-0.08,-1.54,-23.8,-0.09,-6.08,6.16,0.66,-1.52,-1.45,2.39,2.01],
	[9.5,102.5,100931,167,29,16,37,301.1,-0.6,0.2,-0.1,-0.4,17.21,-0.93,-0.04,-0.41,-0.33,-9.3,0,-0.1,0,-0.1,-15.92,0,1.2781,-0.00153,0.00038,-0.0002,-0.00048,0.63699,-0.03243,-0.02455,0.00104,-0.00932,2.5395,0.3155,0.3455,-0.0381,0.1592,287.4,0.7,1.2,-0.3,0,-4.69,1.9,0.42,-0.29,-0.46,7.61,-0.02,0.3,-0.09,-0.08,-1.11,-23.02,1.56,-6.93,6.03,1.25,-0.58,-0.94,2.57,1.81],
	[9.5,103.5,100933,167,32,17,37,301,-0.6,0.2,-0.1,-0.4,17.29,-0.97,-0.1,-0.4,-0.34,-9.1,0.2,-0.2,0,-0.1,-12.25,0,1.27811,-0.00141,0.00039,-0.00019,-0.00048,0.63674,-0.03135,-0.02464,0.00229,-0.00898,2.5475,0.3224,0.3395,-0.0377,0.1542,287.4,0.6,1.2,-0.4,0,-4.62,1.91,0.52,-0.32,-0.5,7.74,0.06,0.3,-0.17,-0.1,-0.98,-21.72,2.73,-7.92,5.26,1.95,1.73,-0.9,2.5,1.02],
	[9.5,104.5,100880,167,32,17,36,300.6,-0.6,0.3,-0.3,-0.3,17.26,-1.1,-0.27,-0.37,-0.42,-8.3,0.5,-0.2,0.2,-0.1,-8.84,4.73,1.27803,-0.00133,0.00047,-0.00022,-0.00043,0.63632,-0.02959,-0.02481,0.00346,-0.00869,2.5395,0.295,0.3176,-0.0428,0.1384,287.4,0.5,1.2,-0.5,0,-4.63,2.13,0.67,-0.34,-0.55,7.56,0.13,0.25,-0.13,-0.15,-1.05,-19.5,3.37,-8.74,4.2,1.31,4.64,-1.62,3.49,0.02],
	[9.5,105.5,100698,168,35,17,35,300.2,-0.7,0.4,-0.4,-0.4,17.16,-1.25,-0.4,-0.37,-0.5,-7.3,0.3,-0.4,0.3,-0.1,-5.42,20.88,1.27798,-0.00125,0.00056,-0.00026,-0.0004,0.63503,-0.02714,-0.02486,0.00493,-0.00859,2.5442,0.2654,0.321,-0.0646,0.1308,287.3,0.4,1.2,-0.5,0,-4.46,2.27,0.73,-0.27,-0.59,7.91,0.06,0.22,-0.12,-0.16,-1.5,-17.62,2.8,-8.48,3.28,-1.22,3.43,-3.88,5.08,-0.66],
	[9.5,106.5,100897,174,44,15,37,300.1,-0.9,0.1,-0.4,-0.6,17.37,-1.15,-0.25,-0.38,-0.45,-7.7,0.2,-0.2,0.2,0.1,-1.6,3.63,1.27803,-0.00121,0.00055,-0.00028,-0.00043,0.63502,-0.0243,-0.02521,0.00668,-0.00843,2.611,0.3025,0.3985,-0.0983,0.1544,287.3,0.2,1.2,-0.6,-0.1,-4.56,2.28,0.73,-0.18,-0.61,8.22,0.02,0.24,-0.12,-0.07,-4.33,-17.88,0.26,-7.04,2.62,-4.33,3.22,-4.92,6.17,-0.12],
	[8.5,96.5,100923,119,-4,5,35,301.1,-0.3,0.5,-0.1,-0.1,17.35,-0.76,0,-0.22,-0.37,-9.3,-0.1,-0.1,0,0,-34.84,0,1.27885,-0.00142,0.00045,-0.00009,-0.0003,0.64197,-0.02931,-0.021,-0.0038,-0.01163,2.5427,0.2651,0.2813,0.0376,0.1986,287.2,0.8,1.2,-0.1,0.2,-3.86,1.8,0.18,-0.12,-0.12,11.09,-0.24,-0.09,0.14,-0.07,-8.58,-27.24,-7.34,-1.74,2.94,1.17,-0.26,1.78,-0.44,1.89],
	[8.5,97.5,100921,115,-6,9,35,301.1,-0.3,0.5,-0.2,-0.1,17.4,-0.73,-0.05,-0.23,-0.42,-9,0.1,-0.1,0,0,-30.19,0,1.27882,-0.00133,0.00049,-0.0001,-0.0003,0.6418,-0.02886,-0.02061,-0.00353,-0.01123,2.548,0.2713,0.2648,0.0395,0.1779,287.3,0.8,1.2,-0.1,0.2,-3.7,1.78,0.19,-0.15,-0.16,10.24,-0.18,0.01,0.1,-0.08,-5.78,-26.21,-6.63,-1.85,3.4,3.11,1.89,1.38,0.41,2.04],
	[8.5,98.5,99670,111,-6,12,32,299.7,-0.3,0.7,-0.3,0,17.1,-0.8,-0.26,-0.21,-0.55,-6.8,0.3,-0.3,-0.2,-0.3,-25.7,110.9,1.27879,-0.00133,0.00057,-0.00017,-0.00029,0.63609,-0.02778,-0.02037,-0.00289,-0.01095,2.5819,0.2559,0.2376,0.0368,0.15,286.7,0.7,1.2,-0.2,0.2,-3.81,1.78,0.24,-0.17,-0.21,8.58,0.29,0.24,0.04,-0.18,-3.67,-23.98,-6.15,-1.85,3.21,2.06,3.22,-0.66,2.34,0.95],
	[8.5,99.5,98816,124,3,14,32,299,-0.5,0.6,-0.3,-0.1,16.82,-0.64,-0.2,-0.25,-0.54,-6.9,-0.1,-0.3,-0.2,-0.3,-21.78,187.4,1.27867,-0.00157,0.00054,-0.00024,-0.00032,0.63082,-0.02619,-0.02084,-0.00138,-0.01076,2.6065,0.2727,0.2773,0.0053,0.1559,286.3,0.6,1.2,-0.2,0.2,-4.42,1.75,0.3,-0.17,-0.25,7.17,0.77,0.39,0.19,-0.2,-4.51,-22.53,-6.05,-1.84,2.99,-2.68,2.86,-4.82,5.41,-0.05],
	[8.5,100.5,100902,148,14,16,37,300.7,-0.8,0.2,-0.1,-0.2,17.36,-0.59,0.01,-0.36,-0.49,-8.3,0.1,-0.1,-0.1,-0.1,-19.02,2.34,1.27839,-0.00174,0.00046,-0.00026,-0.00034,0.63865,-0.02514,-0.02217,0.00054,-0.01056,2.553,0.2531,0.3252,-0.043,0.158,287.3,0.5,1.2,-0.3,0.1,-5.19,1.76,0.39,-0.18,-0.28,6.07,0.11,0.32,0.12,-0.01,-4.54,-22.96,-5.16,-2.69,3.57,-4.54,0.44,-4.54,4.2,0.94],
	[8.5,101.5,100929,151,20,16,41,301.1,-0.8,0.1,-0.1,-0.3,17.25,-0.65,0.02,-0.33,-0.41,-9.4,0,-0.1,0,-0.1,-16.61,0,1.27821,-0.00166,0.00047,-0.00027,-0.00041,0.63819,-0.02547,-0.02324,0.00172,-0.01016,2.5507,0.257,0.3563,-0.0557,0.1721,287.3,0.4,1.1,-0.3,0.1,-5.22,1.79,0.46,-0.19,-0.33,6.65,-0.02,0.27,-0.05,-0.03,-3.59,-23.13,-3.75,-3.56,4.12,0.42,-1.59,-2.72,3.06,1.73],
	[8.5,102.5,100933,151,24,17,43,301.1,-0.8,0.1,-0.1,-0.4,17.22,-0.72,-0.02,-0.33,-0.39,-9.3,0,-0.1,0,-0.1,-13.38,0,1.27814,-0.00154,0.0005,-0.00029,-0.00045,0.63816,-0.02519,-0.02355,0.00271,-0.00949,2.5482,0.2655,0.3668,-0.066,0.1716,287.3,0.4,1.1,-0.4,0,-4.98,1.81,0.52,-0.18,-0.37,7.33,0.02,0.25,-0.09,-0.05,-3.32,-22.95,-2.88,-4.11,4.22,0.78,-1.41,-2.19,3.3,2.13],
	[8.5,103.5,100937,151,28,17,44,300.9,-0.9,0,-0.2,-0.5,17.24,-0.75,-0.05,-0.32,-0.37,-9.3,0,-0.1,0,0,-9.72,0,1.27816,-0.00143,0.00051,-0.00031,-0.00046,0.63799,-0.02423,-0.02353,0.00383,-0.00875,2.5536,0.2768,0.3788,-0.0754,0.1721,287.3,0.3,1.1,-0.4,0,-4.8,1.82,0.57,-0.17,-0.39,7.71,-0.01,0.18,-0.16,-0.08,-2.67,-21.96,-2.12,-4.67,3.84,1.66,0.32,-1.61,3.28,1.83],
	[8.5,104.5,100916,151,31,17,45,300.7,-1,0,-0.3,-0.5,17.29,-0.77,-0.07,-0.3,-0.37,-9.1,0.2,-0.1,0.1,0,-6.13,2.14,1.2781,-0.00133,0.00054,-0.00033,-0.00047,0.63804,-0.02269,-0.02344,0.00508,-0.00825,2.5655,0.283,0.3882,-0.0835,0.1701,287.2,0.2,1.1,-0.5,-0.1,-4.83,1.87,0.61,-0.12,-0.4,7.65,-0.07,0.16,-0.17,-0.1,-2.67,-20.65,-2.17,-4.9,3.16,1.65,2.58,-1.56,3.77,0.94],
	[8.5,105.5,100921,151,34,16,44,300.5,-1.1,-0.1,-0.3,-0.6,17.33,-0.8,-0.1,-0.29,-0.39,-8.8,0.3,0,0.1,0.1,-2.74,1.84,1.27798,-0.00126,0.00056,-0.00034,-0.00047,0.63825,-0.02084,-0.0236,0.00662,-0.00803,2.5779,0.2804,0.3982,-0.0965,0.1664,287.2,0.1,1.1,-0.6,-0.1,-4.97,2,0.68,-0.07,-0.43,7.78,-0.01,0.19,-0.2,-0.11,-4.27,-20.08,-3.67,-4.47,2.15,0.48,3.06,-1.94,4.15,0.32],
	[8.5,106.5,100942,152,38,14,44,300.4,-1.3,-0.2,-0.4,-0.6,17.38,-0.79,-0.09,-0.28,-0.4,-8.9,0.3,0,0.1,0.1,1.08,0,1.27799,-0.00125,0.00053,-0.00035,-0.00049,0.63833,-0.01897,-0.02396,0.0081,-0.00789,2.5914,0.2755,0.4141,-0.114,0.1645,287.1,0,1.1,-0.6,-0.1,-4.91,2.18,0.75,-0.02,-0.48,8.34,-0.09,0.2,-0.18,-0.05,-7.99,-20.62,-6.53,-3.36,1.16,0.65,4.24,-1.42,4.49,0.55],
	[7.5,96.5,100923,99,-10,2,34,301.1,-0.3,0.5,-0.1,0,17.37,-0.61,0.02,-0.16,-0.39,-9.3,0,-0.1,0,-0.1,-32.8,0,1.27879,-0.00145,0.00049,-0.00007,-0.00027,0.64555,-0.02188,-0.01867,-0.00297,-0.01176,2.5062,0.166,0.2433,0.0302,0.1904,287.1,0.5,1.1,-0.1,0.3,-3.08,1.33,0.2,-0.05,0.02,10.98,-0.18,-0.09,0.15,-0.04,-10.06,-25.57,-7.54,-1.54,0.77,-0.31,-0.38,0.95,-0.37,0.99],
	[7.5,97.5,100919,99,-10,5,34,301.2,-0.4,0.5,-0.1,0,17.39,-0.59,0.01,-0.19,-0.41,-9.2,0.1,-0.1,-0.1,-0.1,-28.31,0,1.27871,-0.00138,0.0005,-0.00007,-0.00026,0.64448,-0.02173,-0.01812,-0.00277,-0.01139,2.5166,0.1803,0.2426,0.0297,0.1834,287.1,0.5,1.1,-0.1,0.2,-3.24,1.22,0.2,-0.06,0.01,10.42,-0.11,0.01,0.12,-0.05,-8.73,-24.14,-6.69,-1.55,1.81,0.35,-0.43,0.63,0.44,1.45],
	[7.5,98.5,100844,97,-11,9,34,300.9,-0.4,0.6,-0.2,0,17.4,-0.61,-0.11,-0.2,-0.51,-8.4,0.2,-0.1,-0.1,-0.1,-23.66,6.34,1.27874,-0.00128,0.00058,-0.00012,-0.00025,0.64307,-0.0213,-0.01798,-0.00225,-0.0109,2.5241,0.1852,0.2287,0.024,0.1609,287.1,0.5,1.1,-0.1,0.2,-3.67,1.45,0.29,-0.06,-0.06,9.09,-0.06,0.09,0.04,-0.1,-6.82,-21.08,-6.31,-0.98,2.16,0.88,1.05,-0.75,1.88,1.21],
	[7.5,99.5,99735,102,-5,12,34,299.9,-0.5,0.6,-0.2,0,17.09,-0.58,-0.16,-0.2,-0.56,-7,0.2,-0.2,-0.2,-0.3,-19.43,104.68,1.27867,-0.00139,0.00059,-0.0002,-0.00028,0.63699,-0.02005,-0.01857,-0.00116,-0.01063,2.5596,0.1917,0.2439,0.0083,0.153,286.7,0.4,1.1,-0.2,0.2,-4.31,1.95,0.48,-0.01,-0.15,8.06,0.64,0.31,0.13,-0.14,-6.94,-18.79,-7.45,0.74,1.41,-0.4,2.32,-3.57,4.39,0.51],
	[7.5,100.5,100700,124,4,15,38,300.5,-0.8,0.3,-0.2,-0.1,17.28,-0.49,-0.04,-0.27,-0.54,-7.8,0,-0.2,-0.1,-0.2,-16.21,19.73,1.27852,-0.00161,0.00056,-0.00026,-0.00031,0.63996,-0.01856,-0.01954,0.00042,-0.01049,2.5292,0.1929,0.2839,-0.0311,0.1515,287.1,0.2,1.1,-0.3,0.1,-5.19,1.97,0.56,0.03,-0.16,6.69,0.46,0.32,0.16,-0.02,-10.1,-19.49,-9.38,1.75,0.95,-3.78,1.13,-5.08,5,0.46],
	[7.5,101.5,100929,136,14,17,43,300.8,-0.9,0.1,-0.1,-0.3,17.32,-0.55,0.03,-0.29,-0.47,-8.8,-0.1,-0.1,-0.1,-0.1,-13.8,0,1.27828,-0.00167,0.00052,-0.0003,-0.00038,0.64016,-0.01845,-0.02061,0.00189,-0.01006,2.5406,0.1884,0.3285,-0.0604,0.1642,287.2,0.1,1,-0.3,0.1,-5.76,1.41,0.42,-0.06,-0.18,6.57,0.18,0.34,0,0.04,-10.12,-21.12,-9.25,0.31,1.16,-1.1,-0.59,-3.65,4.45,1.47],
	[7.5,102.5,100935,137,20,17,47,300.9,-0.9,0,-0.2,-0.4,17.29,-0.63,0.02,-0.28,-0.43,-9.2,0,-0.1,0,-0.1,-10.9,0,1.2781,-0.00158,0.00052,-0.00033,-0.00043,0.64026,-0.01836,-0.02114,0.00298,-0.00935,2.5455,0.1918,0.3504,-0.075,0.1705,287.2,0.1,1,-0.4,0,-5.57,1.25,0.36,-0.08,-0.26,6.99,-0.04,0.21,-0.14,-0.01,-8.44,-22.27,-8.9,-0.66,1.55,0.45,-1.4,-2.4,3.38,1.91],
	[7.5,103.5,100940,136,24,18,49,300.8,-1,-0.1,-0.2,-0.5,17.25,-0.61,0,-0.26,-0.4,-9.3,0.1,-0.1,0,0,-7.54,0,1.27805,-0.00146,0.00054,-0.00035,-0.00046,0.6404,-0.01767,-0.02106,0.00405,-0.00842,2.5421,0.2125,0.3652,-0.0827,0.1734,287.1,0,1,-0.4,-0.1,-5.04,1.4,0.41,-0.01,-0.29,7.34,-0.12,0.12,-0.2,-0.07,-8.45,-22.09,-9.19,-0.29,1.77,1.27,-0.78,-1.83,3.35,2.19],
	[7.5,104.5,100942,134,26,17,50,300.7,-1.1,-0.2,-0.3,-0.6,17.31,-0.61,0.01,-0.25,-0.38,-9.2,0.2,0,0,0.1,-3.95,0,1.27802,-0.00135,0.00054,-0.00037,-0.00049,0.64047,-0.01651,-0.02083,0.00526,-0.00769,2.5529,0.2245,0.379,-0.0915,0.1728,287.1,0,1,-0.5,-0.1,-4.76,1.56,0.5,0.03,-0.3,7.71,-0.19,0.09,-0.21,-0.08,-8.26,-21.1,-9.54,-0.1,1.28,1.84,1.42,-0.95,3.73,1.73],
	[7.5,105.5,100943,132,28,16,50,300.5,-1.3,-0.2,-0.3,-0.6,17.35,-0.6,0.01,-0.23,-0.37,-9.2,0.3,0.1,0.1,0.1,-0.07,0,1.27801,-0.0013,0.00052,-0.00039,-0.00051,0.64087,-0.01519,-0.02095,0.00664,-0.0073,2.554,0.2229,0.3802,-0.1005,0.1655,287,-0.1,1,-0.5,-0.1,-4.98,1.68,0.58,0.06,-0.32,8.02,-0.15,0.13,-0.23,-0.08,-9.19,-20.46,-10.47,0,0.48,2.82,3,-0.12,3.65,1.34],
	[7.5,106.5,100941,129,30,13,48,300.5,-1.3,-0.2,-0.4,-0.6,17.37,-0.58,-0.01,-0.21,-0.38,-9.1,0.3,0.1,0.1,0.1,3.95,0,1.27797,-0.00129,0.00048,-0.0004,-0.00053,0.64134,-0.01418,-0.02133,0.00775,-0.00717,2.5424,0.2067,0.3673,-0.1098,0.1523,287,-0.1,0.9,-0.6,-0.1,-5.08,1.75,0.66,0.07,-0.31,8.14,-0.21,0.1,-0.21,-0.04,-10.71,-19.98,-11.39,0.35,-0.09,3.5,4.4,0.59,3.53,1.11],
	[6.5,96.5,100921,88,-14,0,32,301.1,-0.5,0.4,0,0,17.47,-0.48,0.05,-0.12,-0.39,-9.1,0,-0.1,0,-0.1,-31.65,0,1.27851,-0.00156,0.0005,-0.00007,-0.00027,0.64818,-0.01578,-0.01655,-0.00229,-0.01137,2.4827,0.075,0.2034,0.0236,0.1706,286.9,0.3,1,0,0.2,-3.06,0.75,0.05,0,0.07,11.37,-0.14,-0.1,0.19,-0.02,-12.6,-23.17,-9.63,-1.73,-3.75,-1.52,-0.76,-0.09,-0.09,-0.2],
	[6.5,97.5,100915,88,-13,2,33,301.2,-0.5,0.4,0,0,17.46,-0.51,0.04,-0.15,-0.4,-9.2,0,-0.1,-0.1,-0.1,-26.78,0,1.27846,-0.00146,0.00051,-0.00007,-0.00026,0.64815,-0.01572,-0.01585,-0.00221,-0.01119,2.4814,0.0879,0.2064,0.0198,0.1691,287,0.3,1,-0.1,0.2,-3.57,0.85,0.08,-0.03,0.02,10.75,-0.15,0,0.11,-0.01,-13.78,-22.28,-9.28,-1.88,-2.91,1.64,-2.36,0.54,0.23,0.8],
	[6.5,98.5,100912,84,-13,5,34,301.3,-0.4,0.4,-0.1,0,17.46,-0.52,0.03,-0.19,-0.44,-9.2,0.1,-0.1,-0.1,0,-21.7,0,1.27867,-0.00128,0.00055,-0.00008,-0.00027,0.64727,-0.01543,-0.01519,-0.00223,-0.01044,2.4784,0.111,0.2103,0.0148,0.1594,287.1,0.3,1,-0.1,0.2,-3.86,0.82,0.06,0,-0.03,9.5,-0.22,0.02,0.01,-0.03,-16.2,-20.87,-10.32,-1.14,-2.85,0.22,-1.43,0.05,0.66,1.44],
	[6.5,99.5,100723,83,-14,7,35,300.9,-0.4,0.5,-0.1,0,17.34,-0.54,-0.06,-0.19,-0.53,-8.5,0.1,-0.1,-0.1,-0.1,-17.42,16.55,1.27864,-0.00118,0.00062,-0.0001,-0.00028,0.64489,-0.0147,-0.01534,-0.00203,-0.00994,2.473,0.1172,0.2005,0.0139,0.1368,287,0.2,1,-0.1,0.2,-4.17,1.07,0.16,0.11,-0.06,8.71,0.05,0.09,-0.01,-0.08,-15.92,-19.72,-12.43,0.94,-3.05,0.41,1.85,-0.98,1.82,1.22],
	[6.5,100.5,99058,89,-8,10,35,299.4,-0.7,0.5,-0.2,-0.1,16.65,-0.46,-0.18,-0.16,-0.53,-6.6,0.3,-0.1,-0.1,-0.1,-13.55,165.83,1.27863,-0.00135,0.00066,-0.0002,-0.00034,0.63712,-0.01356,-0.01596,-0.00116,-0.00972,2.4782,0.1287,0.2046,0.0012,0.1343,286.2,0.1,1,-0.2,0.1,-5.11,1.42,0.32,0.16,-0.09,7.85,0.53,0.24,0.16,-0.05,-17.84,-18.36,-13.71,2.65,-2.81,1.03,1.15,-3.92,4.61,0.24],
	[6.5,101.5,99384,113,5,15,39,299.5,-1.1,0.2,-0.2,-0.2,16.59,-0.24,0,-0.21,-0.53,-6.4,0,-0.1,-0.2,-0.1,-10.45,138.19,1.2784,-0.00167,0.00056,-0.00032,-0.00038,0.63673,-0.01286,-0.01717,0.00065,-0.00941,2.4604,0.1669,0.2694,-0.0407,0.1375,286.4,-0.1,0.9,-0.3,0.1,-6.18,1.14,0.28,0.06,-0.1,6.74,0.81,0.49,0.17,0.07,-18.75,-17.57,-13.13,1.98,-1.85,-4.85,-0.02,-5.23,6.52,0.46],
	[6.5,102.5,100889,130,17,18,48,300.6,-1.1,-0.1,-0.2,-0.3,17.21,-0.43,0.11,-0.24,-0.5,-8.4,-0.2,-0.1,-0.1,0,-7.82,4.75,1.2781,-0.00172,0.00042,-0.00036,-0.00042,0.64231,-0.01245,-0.01828,0.00248,-0.00896,2.4787,0.1336,0.3008,-0.0678,0.1417,287,-0.2,0.9,-0.3,0,-6.48,0.73,0.2,-0.06,-0.14,6.66,-0.19,0.14,-0.19,0.05,-15.35,-18.44,-11.75,0.45,-0.64,-1.29,-0.6,-2.95,3.61,1.27],
	[6.5,103.5,100940,125,20,17,52,300.7,-1,-0.1,-0.2,-0.5,17.36,-0.59,0.04,-0.21,-0.42,-9.1,0,-0.1,0,0,-5.16,0,1.27791,-0.00152,0.00046,-0.00037,-0.00047,0.64348,-0.01227,-0.0183,0.00365,-0.0081,2.5109,0.1247,0.3074,-0.0718,0.1595,286.9,-0.2,0.8,-0.4,-0.1,-5.83,0.93,0.28,0,-0.17,6.5,-0.32,0.05,-0.27,-0.07,-12.86,-19.14,-11.13,0.39,0.19,0.35,-1.07,-1.56,2.82,1.77],
	[6.5,104.5,100941,119,21,16,53,300.6,-1.1,-0.2,-0.3,-0.6,17.34,-0.54,0.05,-0.2,-0.4,-9.2,0.2,0,0,0,-1.95,0,1.27777,-0.00139,0.00047,-0.00038,-0.00051,0.64344,-0.01175,-0.01805,0.00467,-0.00712,2.5087,0.1439,0.3177,-0.0775,0.1589,286.9,-0.2,0.8,-0.4,-0.1,-5.02,1.17,0.38,0.06,-0.16,7.41,-0.4,-0.03,-0.25,-0.07,-11.73,-18.94,-11.61,0.94,0.5,1.49,0.52,-0.35,3.51,2.31],
	[6.5,105.5,100940,113,21,14,52,300.6,-1.2,-0.2,-0.4,-0.6,17.38,-0.5,0.05,-0.18,-0.38,-9.1,0.3,0,0.1,0.1,2.37,0,1.27788,-0.00132,0.00045,-0.00041,-0.00054,0.64369,-0.01086,-0.01785,0.00561,-0.00659,2.5069,0.1461,0.3162,-0.0868,0.1496,286.9,-0.2,0.8,-0.5,-0.1,-4.7,1.37,0.49,0.09,-0.17,8.19,-0.27,0.04,-0.25,-0.07,-11.16,-18.17,-12.11,1.41,0.18,3.53,2.51,0.74,3.29,1.88],
	[6.5,106.5,100937,108,22,11,50,300.6,-1.3,-0.2,-0.4,-0.5,17.4,-0.45,0.03,-0.15,-0.38,-9.1,0.3,0,0.1,0.1,6.73,0,1.27792,-0.00131,0.00042,-0.00043,-0.00056,0.64443,-0.01026,-0.01786,0.00627,-0.00647,2.493,0.1359,0.302,-0.094,0.1366,286.9,-0.2,0.8,-0.5,-0.1,-4.7,1.49,0.56,0.1,-0.19,8,-0.28,0.02,-0.22,-0.05,-10.99,-17.6,-12.76,2.16,-0.16,4.07,3.9,1.36,2.76,1.48],
	[5.5,96.5,99579,71,-14,-3,29,299.7,-0.6,0.3,-0.1,0,17.01,-0.27,0.1,-0.08,-0.34,-7.3,0,-0.2,-0.1,-0.2,-29.21,123.08,1.27872,-0.0016,0.00047,-0.00008,-0.00036,0.64721,-0.01097,-0.01343,-0.00236,-0.00935,2.4102,0.0185,0.1547,0.0235,0.1311,286.1,0.1,0.9,0,0.2,-4.77,0.71,0.08,0.07,0.06,11.35,0.16,0.02,0.23,-0.01,-28.94,-17.24,-12.96,0.42,-8.74,-0.79,-1.5,-1.84,0.55,-1.41],
	[5.5,97.5,100821,81,-13,1,31,300.8,-0.6,0.3,-0.1,0,17.35,-0.38,0.08,-0.11,-0.36,-8.2,0,-0.1,0,0,-24.49,10.89,1.27828,-0.00151,0.00049,-0.00011,-0.00034,0.65118,-0.01053,-0.01332,-0.00185,-0.00989,2.3756,0.0299,0.1612,0.0129,0.1292,286.8,0,0.9,0,0.2,-3.96,0.33,-0.08,-0.02,0.05,10.26,-0.19,0,0.16,0.04,-26.99,-15.2,-9.72,-1.24,-7.15,-0.68,-4.81,-1.53,1.4,-0.67],
	[5.5,98.5,100915,76,-14,3,34,301.1,-0.4,0.3,-0.1,0,17.67,-0.48,0.08,-0.14,-0.38,-8.9,0,0,-0.1,0,-19.64,0,1.27826,-0.00125,0.00055,-0.00012,-0.00035,0.65183,-0.01048,-0.01265,-0.00183,-0.00946,2.4225,0.0448,0.1708,0.0061,0.1286,286.9,0,0.8,-0.1,0.2,-3.53,0.54,-0.01,-0.05,-0.01,9.64,-0.37,-0.04,-0.05,0.02,-21.27,-14.49,-8.04,-2.94,-5.53,-0.64,-1.07,1.55,0.38,1.22],
	[5.5,99.5,100908,68,-17,2,34,301.1,-0.3,0.4,-0.1,0,17.6,-0.48,0.05,-0.16,-0.44,-8.7,0.1,0,0,0,-15.42,0,1.27833,-0.00112,0.00057,-0.00008,-0.00037,0.65041,-0.0097,-0.01174,-0.00222,-0.0088,2.4162,0.0452,0.155,0.0105,0.1111,287,0.1,0.9,-0.1,0.2,-3.82,0.64,-0.03,0.04,-0.06,9.11,-0.16,0.08,-0.14,-0.04,-18.67,-17.84,-10.97,-2.78,-5.35,-1.43,3.76,2.6,-0.15,1.43],
	[5.5,100.5,98290,65,-14,4,33,299.1,-0.4,0.5,-0.1,0,16.63,-0.51,-0.16,-0.15,-0.51,-6.8,0.3,-0.1,0,-0.1,-11.01,234.81,1.27842,-0.00123,0.0006,-0.00013,-0.00042,0.63799,-0.00935,-0.01209,-0.00188,-0.0086,2.4371,0.0243,0.1242,0.0091,0.0992,285.8,0,0.8,-0.1,0.1,-4.25,0.63,-0.03,0.12,-0.05,8.89,0.18,0.1,0.05,-0.07,-16.42,-18.35,-13.45,-0.65,-4.34,4.31,0.3,-2.12,2.75,-0.41],
	[5.5,101.5,94927,78,1,8,34,297,-0.9,0.3,-0.3,-0.1,15.42,-0.09,-0.06,-0.12,-0.51,-5.6,0.1,-0.1,-0.1,-0.2,-6.84,541.32,1.27843,-0.00152,0.0005,-0.00028,-0.00048,0.62042,-0.00934,-0.01361,-0.00024,-0.00835,2.4875,0.1303,0.1986,-0.0136,0.1202,284.3,-0.1,0.8,-0.2,0,-4.84,0.61,0.02,0.14,-0.04,6.94,0.71,0.28,0.16,-0.01,-15.51,-15.29,-12.93,1.53,-2.97,-5.66,-1.9,-6.24,5.76,-0.43],
	[5.5,102.5,98941,113,14,16,43,299.1,-1.3,-0.1,-0.3,-0.2,16.35,0.01,0.21,-0.16,-0.56,-6.1,-0.4,-0.1,-0.1,-0.1,-3.88,178.82,1.27824,-0.00178,0.00036,-0.0004,-0.00045,0.63726,-0.00868,-0.01501,0.00169,-0.00825,2.3877,0.1457,0.2597,-0.0482,0.1105,286,-0.3,0.7,-0.3,0,-5.74,0.74,0.08,0.15,0.02,6.89,-0.12,0.07,-0.16,0.02,-17.53,-14.7,-13.23,2.28,-2.29,-2.83,-0.79,-4.41,4.4,-0.01],
	[5.5,103.5,100869,116,18,16,51,300.4,-1,-0.2,-0.3,-0.4,17.36,-0.49,0.07,-0.18,-0.49,-8.5,-0.1,-0.1,-0.1,0,-1.86,6.56,1.27789,-0.00157,0.00036,-0.00039,-0.00047,0.64664,-0.00779,-0.01535,0.00307,-0.00776,2.4428,0.0646,0.2423,-0.0592,0.1279,286.8,-0.4,0.7,-0.4,-0.1,-6.13,0.81,0.17,0.1,0,5.68,-0.57,-0.09,-0.29,-0.01,-16.04,-14.95,-11.97,1.35,-1.04,-3.05,-0.23,-2.24,2.94,0.86],
	[5.5,104.5,100938,105,16,14,53,300.6,-1,-0.1,-0.3,-0.5,17.41,-0.53,0.04,-0.16,-0.4,-9.1,0.1,0,0,0,0.68,0,1.27765,-0.00137,0.00042,-0.00039,-0.00053,0.64655,-0.00743,-0.01507,0.00387,-0.007,2.4618,0.0691,0.2502,-0.0649,0.1436,286.8,-0.3,0.7,-0.4,-0.1,-5.51,0.99,0.32,0.1,-0.06,6.74,-0.58,-0.13,-0.29,-0.08,-12.47,-15.89,-11.14,1.17,-0.35,0.14,0.47,0.02,2.98,2.16],
	[5.5,105.5,100936,96,15,12,52,300.7,-1.1,-0.2,-0.4,-0.5,17.41,-0.43,0.07,-0.13,-0.38,-9.2,0.3,0,0,0.1,4.67,0,1.27766,-0.0013,0.00042,-0.00042,-0.00057,0.64659,-0.00692,-0.01441,0.0043,-0.00649,2.4588,0.0805,0.2526,-0.0712,0.1357,286.8,-0.3,0.7,-0.4,-0.1,-4.84,1.05,0.39,0.12,-0.06,7.83,-0.36,-0.03,-0.25,-0.04,-11.16,-16.27,-11.77,1.95,-0.25,2.93,2.23,1.52,2.55,2.2],
	[5.5,106.5,100932,90,14,9,50,300.6,-1.2,-0.2,-0.4,-0.5,17.42,-0.35,0.07,-0.1,-0.37,-9.2,0.3,0,0.1,0.1,9.12,0,1.2777,-0.00128,0.00039,-0.00044,-0.00059,0.6473,-0.00654,-0.01402,0.00464,-0.00625,2.4485,0.0776,0.2418,-0.0757,0.1243,286.8,-0.3,0.7,-0.4,-0.1,-4.64,1.07,0.44,0.12,-0.06,8.12,-0.33,-0.02,-0.23,-0.06,-10.45,-16.13,-12.75,2.79,-0.52,3.71,3.51,2.23,1.9,1.75],
	[4.5,96.5,92266,36,-13,-6,23,295.2,-0.4,0.3,-0.1,-0.1,15.2,0,0.18,-0.09,-0.32,-5.4,-0.2,0,-0.1,-0.2,-27.12,792,1.27886,-0.00143,0.00047,-0.00005,-0.00056,0.62214,-0.01031,-0.01145,-0.00274,-0.00703,2.5494,0.0029,0.14,0.0347,0.0971,282.6,0.2,0.8,0,0.1,-0.86,0.63,0.06,0.08,0.02,12.67,0.11,0,0.11,-0.04,-15.1,-11.94,-9.93,1.54,-6.15,8.6,-1.94,-1.72,0.29,-1.1],
	[4.5,97.5,95312,54,-10,-2,26,297.2,-0.6,0.3,-0.1,-0.1,15.89,-0.14,0.03,-0.07,-0.34,-6,-0.1,0,-0.1,-0.1,-21.86,508.88,1.27836,-0.00138,0.00052,-0.00011,-0.0005,0.63088,-0.00769,-0.01101,-0.00222,-0.00801,2.3961,0.0204,0.1239,0.026,0.094,284.2,0,0.7,0,0.1,-2.87,0.5,0.01,0.14,0.03,8.92,0.29,0.16,0.15,0.01,-24.07,-11.72,-9.86,2.13,-6.6,-7.17,-4.67,-4.25,2.39,-2.45],
	[4.5,98.5,100897,70,-15,2,33,300.8,-0.4,0.3,-0.1,-0.1,17.55,-0.32,0.08,-0.1,-0.34,-8.1,0,0,0,0,-17.95,4.35,1.27788,-0.00122,0.00057,-0.00016,-0.00044,0.6542,-0.00633,-0.01083,-0.00124,-0.00847,2.3098,0.0226,0.1326,0.0067,0.094,286.8,-0.1,0.7,-0.1,0.1,-2.93,0.1,-0.1,0.15,0.06,10.14,-0.45,-0.11,-0.09,0.03,-19.92,-13.15,-7.81,1.15,-5.81,-7.75,0.12,0.41,2.24,-0.19],
	[4.5,99.5,100906,64,-17,0,35,301,-0.3,0.3,-0.1,-0.1,18.07,-0.32,0.17,-0.1,-0.33,-8.6,0,0,0,0,-13.49,0,1.27793,-0.00116,0.00052,-0.00013,-0.00048,0.65371,-0.00529,-0.0097,-0.00116,-0.0079,2.4235,0.0117,0.1376,0.0059,0.0968,286.9,-0.1,0.7,-0.1,0.1,-3.44,0.07,-0.05,0.07,0.01,9.5,-0.16,0.12,-0.1,0.05,-12.73,-13.25,-6.22,-1.58,-4.74,-4.3,6.1,5.12,0.6,1.56],
	[4.5,100.5,99961,57,-16,1,34,300.3,-0.3,0.4,-0.1,-0.1,17.45,-0.31,0.08,-0.12,-0.41,-6.9,-0.2,0,0,-0.2,-8.95,85.88,1.27802,-0.00124,0.00046,-0.00012,-0.0005,0.64775,-0.00466,-0.00886,-0.00123,-0.00754,2.3867,-0.015,0.0994,0.0064,0.078,286.5,-0.1,0.7,-0.1,0.1,-3.69,0.31,-0.02,0.06,-0.05,9.38,-0.06,0.07,-0.05,-0.08,-8.37,-15,-9.42,-1.92,-4.8,2.77,1.98,2.04,0.66,-0.87],
	[4.5,101.5,95271,62,-4,5,33,297.3,-0.6,0.3,-0.2,-0.1,15.76,-0.18,-0.04,-0.11,-0.47,-5.4,-0.3,-0.2,-0.1,-0.2,-4.51,509.08,1.27824,-0.00141,0.00041,-0.00023,-0.00055,0.6246,-0.00589,-0.01029,-0.00038,-0.00773,2.449,0.0203,0.1101,-0.0003,0.0878,284.4,-0.2,0.7,-0.2,0,-3.93,0.61,0.12,0.17,0,7.73,0.41,0.18,0.07,-0.04,-7.25,-13.83,-11.63,0.16,-3.5,-3.32,-4.23,-5.58,3.63,-1.48],
	[4.5,102.5,96241,81,6,10,38,297.6,-1.1,0,-0.3,-0.2,15.81,-0.03,0.07,-0.11,-0.54,-5.4,0,-0.1,-0.1,-0.3,-1.02,420.68,1.27821,-0.0015,0.00038,-0.00034,-0.00053,0.62789,-0.00495,-0.01113,0.00092,-0.00778,2.4134,0.0931,0.1783,-0.022,0.0939,284.8,-0.3,0.6,-0.3,0,-4.25,0.9,0.21,0.24,0.02,6.87,0.01,0.01,-0.09,-0.03,-9.48,-12.37,-12.55,1.88,-2.52,-2.97,-1.71,-5.12,4.19,-1.08],
	[4.5,103.5,100571,99,12,14,49,299.8,-1,-0.1,-0.3,-0.3,17.27,-0.38,0.06,-0.17,-0.56,-7.4,-0.1,-0.2,-0.1,-0.1,1.16,32.93,1.27778,-0.00148,0.00033,-0.00037,-0.00049,0.64807,-0.00374,-0.01192,0.0024,-0.0076,2.3834,0.0309,0.1788,-0.0478,0.0956,286.5,-0.5,0.6,-0.3,-0.1,-4.64,1.02,0.29,0.18,-0.01,5.77,-0.6,-0.2,-0.21,0.04,-10.4,-12.89,-12.75,2.04,-1.96,-5.7,0.85,-3.91,4.02,-0.16],
	[4.5,104.5,100935,90,10,12,53,300.5,-0.9,-0.1,-0.3,-0.4,17.5,-0.5,0.02,-0.15,-0.42,-8.8,0,-0.1,-0.1,0,3.62,0,1.27758,-0.00129,0.00042,-0.00039,-0.00054,0.64917,-0.00317,-0.01177,0.00323,-0.00712,2.4258,0.0135,0.1907,-0.0588,0.125,286.7,-0.4,0.6,-0.3,-0.1,-5.02,0.85,0.29,0.11,-0.01,6.39,-0.53,-0.16,-0.29,-0.08,-9.97,-13.9,-11.71,1.43,-1.01,-1.34,1.34,-0.16,2.61,1.5],
	[4.5,105.5,100932,81,8,9,52,300.7,-1,-0.1,-0.3,-0.5,17.42,-0.37,0.07,-0.11,-0.38,-9.3,0.2,0,0,0,7.24,0,1.27751,-0.00124,0.00042,-0.00041,-0.00058,0.6493,-0.00307,-0.01098,0.00338,-0.00656,2.416,0.0274,0.1976,-0.0595,0.1228,286.7,-0.4,0.6,-0.3,-0.1,-5.02,0.78,0.31,0.08,-0.01,7.45,-0.39,-0.1,-0.25,-0.04,-9.46,-14.54,-11.47,1.76,-0.59,2.49,2.44,2.12,1.8,2.31],
	[4.5,106.5,100929,74,8,6,50,300.7,-1.1,-0.1,-0.3,-0.4,17.42,-0.27,0.1,-0.08,-0.35,-9.1,0.2,0,0.1,0,11.62,0,1.27758,-0.00124,0.00038,-0.00042,-0.00061,0.64992,-0.00282,-0.01023,0.00343,-0.00606,2.4088,0.0267,0.1892,-0.0607,0.1132,286.6,-0.4,0.6,-0.3,-0.1,-5.02,0.75,0.33,0.07,0.02,8.29,-0.28,-0.03,-0.23,-0.08,-9.6,-14.73,-12.01,2.44,-0.71,3.54,3.34,3.17,1.13,2.01]
	]

	numData = len(data)
	plat = dlat*180.0/pi
	plon = dlon*180.0/pi
	indx=[0,0,0,0]
	flat = floor(plat-0.5)-4
	flon =106-floor(plon-0.5)

	# transform to polar distance in degrees
	ppod = (-dlat + pi/2)*180/pi;
	# find the index (line in the grid file) of the nearest point
	# changed for the 1 degree grid (GP)
	ipod = floor((ppod+1));
	ilon = floor((plon+1));
	# normalized (to one) differences, can be positive or negative
	# changed for the 1 degree grid (GP)
	diffpod = (ppod - (ipod - 0.5));
	difflon = (plon - (ilon - 0.5));

	indx[0] = (numData-flat*11-1-flon) # LL
	indx[1] = indx[0] - 11                     # UL
	indx[2] = indx[0] + 1                      # LR
	indx[3] = indx[1] + 1                      # UR

	# change the reference epoch to January 1 2000
	dmjd1 = dmjd-51544.5;

	#print(dmjd1)
	# mean gravity in m/s**2
	gm = 9.80665;
	# molar mass of dry air in kg/mol
	dMtr = 28.965*0.001;
	# universal gas constant in J/K/mol
	Rg = 8.3143;
	# factors for amplitudes
	if (it==1):# then  constant parameters
		cosfy = 0;
		coshy = 0;
		sinfy = 0;
		sinhy = 0;
	else:
		cosfy = cos(dmjd1/365.25*2*pi);
		coshy = cos(dmjd1/365.25*4*pi);
		sinfy = sin(dmjd1/365.25*2*pi);
		sinhy = sin(dmjd1/365.25*4*pi);

	# initialization
	coorgrid = zeros([numData,2]);
	pgrid = zeros([numData, 5]);
	Tgrid = zeros([numData, 5]);
	Qgrid = zeros([numData, 5]);
	dTgrid = zeros([numData, 5]);
	u = zeros([numData, 1]);
	Hs = zeros([numData, 1]);
	ahgrid = zeros([numData, 5]);
	awgrid = zeros([numData, 5]);
	lagrid = zeros([numData, 5]);
	Tmgrid = zeros([numData, 5]);

	# initialization of new vectors
	p =  0.0
	T =  0.0
	dT = 0.0
	Tm = 0.0
	e =  0.0
	ah = 0.0
	aw = 0.0
	la = 0.0
	undu = 0.0


	# loop over grid points
	# for the 1 degree grid (GP)
	for n in range(0,numData,1):  # Thailand AreaGrid
		# read mean values and amplitudes
		vec = data[n]
		pgrid[n]  = vec[2:7];          # pressure in Pascal
		Tgrid[n]  = vec[7:12];         # temperature in Kelvin
		Qgrid[n]  = array(vec[12:17])/1000;  # specific humidity in kg/kg
		dTgrid[n] = array(vec[17:22])/1000;  # temperature lapse rate in Kelvin/m
		u[n]        = vec[22];           # geoid undulation in m
		Hs[n]         = vec[23];           # orthometric grid height in m
		ahgrid[n] = array(vec[24:29])/1000;  # hydrostatic mapping function coefficient, dimensionless
		awgrid[n] = array(vec[29:34])/1000;  # wet mapping function coefficient, dimensionless
		lagrid[n] = vec[34:39];    	   # water vapor decrease factor, dimensionless
		Tmgrid[n] = vec[39:44];    # mean temperature in Kelvin
	undul = [0.0,0.0,0.0,0.0]
	Ql =  [0.0,0.0,0.0,0.0]
	dTl =  [0.0,0.0,0.0,0.0]
	Tl =  [0.0,0.0,0.0,0.0]
	pl 	=  [0.0,0.0,0.0,0.0]
	ahl =  [0.0,0.0,0.0,0.0]
	awl =  [0.0,0.0,0.0,0.0]
	lal =  [0.0,0.0,0.0,0.0]
	Tml =  [0.0,0.0,0.0,0.0]
	el  = [0.0,0.0,0.0,0.0]
	for l in range(0,4):
		#transforming ellipsoidal height to orthometric height:
		#Hortho = -N + Hell
		undul[l] = u[indx[l]]
		hgt = hell if is_msl==1 else hell-undul[l]
		#pressure, temperature at the height of the grid
		T0 = Tgrid[indx[l],0] + Tgrid[indx[l],1]*cosfy + Tgrid[indx[l],2]*sinfy + Tgrid[indx[l],3]*coshy + Tgrid[indx[l],4]*sinhy
		p0 = pgrid[indx[l],0] + pgrid[indx[l],1]*cosfy + pgrid[indx[l],2]*sinfy + pgrid[indx[l],3]*coshy + pgrid[indx[l],4]*sinhy

		#humidity
		Ql[l]= Qgrid[indx[l],0] + Qgrid[indx[l],1]*cosfy + Qgrid[indx[l],2]*sinfy + Qgrid[indx[l],3]*coshy + Qgrid[indx[l],4]*sinhy

		#reduction = stationheight - gridheight
		Hs1 = Hs[indx[l]]
		redh = hgt - Hs1

		#lapse rate of the temperature in degree / m
		dTl[l] = dTgrid[indx[l],0] + dTgrid[indx[l],1]*cosfy + dTgrid[indx[l],2]*sinfy + dTgrid[indx[l],3]*coshy + dTgrid[indx[l],4]*sinhy

		#temperature reduction to station height
		Tl[l] = T0 + dTl[l]*redh - 273.15

		#virtual temperature
		Tv = T0*(1+0.6077*Ql[l])
		c = gm*dMtr/(Rg*Tv)

		#pressure in hPa
		pl[l] = (p0*exp(-c*redh))/100


		#hydrostatic coefficient ah
		ahl[l] = ahgrid[indx[l],0] + ahgrid[indx[l],1]*cosfy + ahgrid[indx[l],2]*sinfy + ahgrid[indx[l],3]*coshy + ahgrid[indx[l],4]*sinhy;

		#wet coefficient aw
		awl[l] = awgrid[indx[l],0] + awgrid[indx[l],1]*cosfy + awgrid[indx[l],2]*sinfy + awgrid[indx[l],3]*coshy + awgrid[indx[l],4]*sinhy;

		#water vapor decrease factor la - added by GP
		lal[l] = lagrid[indx[l],0] + lagrid[indx[l],1]*cosfy + lagrid[indx[l],2]*sinfy + lagrid[indx[l],3]*coshy + lagrid[indx[l],4]*sinhy;

		#mean temperature of the water vapor Tm - added by GP
		Tml[l] = Tmgrid[indx[l],0] + Tmgrid[indx[l],1]*cosfy + Tmgrid[indx[l],2]*sinfy + Tmgrid[indx[l],3]*coshy + Tmgrid[indx[l],4]*sinhy;

		#water vapor pressure in hPa - changed by GP
		e0 = Ql[l]*p0/(0.622+0.378*Ql[l])/100; # on the grid
		el[l] = e0*pow((100*pl[l]/p0),(lal[l]+1));  # on the station height - [14] Askne and Nordius, 1987


	dnpod1 = fabs(diffpod); # distance nearer point
	dnpod2 = 1 - dnpod1;  # distance to distant point
	dnlon1 = fabs(difflon);
	dnlon2 = 1 - dnlon1;
	#print(plon)
	#print(ilon)
	# pressure
	R1 = dnpod2*pl[0]+dnpod1*pl[1];
	R2 = dnpod2*pl[2]+dnpod1*pl[3];
	p = dnlon2*R1+dnlon1*R2;

	# temperature
	R1 = dnpod2*Tl[0]+dnpod1*Tl[1];
	R2 = dnpod2*Tl[2]+dnpod1*Tl[3];
	T = dnlon2*R1+dnlon1*R2;

	# temperature in degree per km
	R1 = dnpod2*dTl[0]+dnpod1*dTl[1];
	R2 = dnpod2*dTl[2]+dnpod1*dTl[3];
	dT = (dnlon2*R1+dnlon1*R2)*1000;

	# water vapor pressure in hPa - changed by GP
	R1 = dnpod2*el[0]+dnpod1*el[1];
	R2 = dnpod2*el[2]+dnpod1*el[3];
	e = dnlon2*R1+dnlon1*R2;

	# hydrostatic
	R1 = dnpod2*ahl[0]+dnpod1*ahl[1];
	R2 = dnpod2*ahl[2]+dnpod1*ahl[3];
	ah = dnlon2*R1+dnlon1*R2;

	# wet
	R1 = dnpod2*awl[0]+dnpod1*awl[1];
	R2 = dnpod2*awl[2]+dnpod1*awl[3];
	aw = dnlon2*R1+dnlon1*R2;

	# undulation
	R1 = dnpod2*undul[0]+dnpod1*undul[1];
	R2 = dnpod2*undul[2]+dnpod1*undul[3];
	undu = dnlon2*R1+dnlon1*R2;

	# water vapor decrease factor la - added by GP
	R1 = dnpod2*lal[0]+dnpod1*lal[1];
	R2 = dnpod2*lal[2]+dnpod1*lal[3];
	la = dnlon2*R1+dnlon1*R2;

	# mean temperature of the water vapor Tm - added by GP
	R1 = dnpod2*Tml[0]+dnpod1*Tml[1];
	R2 = dnpod2*Tml[2]+dnpod1*Tml[3];
	Tm = dnlon2*R1+dnlon1*R2;

	hort = hgt if is_msl==1 else (hell - undu)[0]
	return (p,T[0],dT,Tm,e,ah,aw,la,hort)