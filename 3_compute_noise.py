#!/usr/bin/env python

import pickle
from glob import glob
from obspy.signal import PPSD
from obspy import read_inventory
import numpy as np
import pandas as pd

def rms(s, f):
	    # Parseval: the RMS in time domain is the sqrt of the integral of the power spectrum
	return np.sqrt(np.trapz(s, f))

fmin=5
fmax=15


with open('params_BGS_IRIS.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
	datelist, station_list, dataset, nslc_all = pickle.load(f)
with open('params_RASPISHAKE.pkl', 'rb') as f:  # Python 3: open(..., 'wb')
        datelist, station_list2, dataset, nslc_all2 = pickle.load(f)
with open('params_BGS_EIDA.pkl', 'rb') as f:
	 datelist, station_list3, dataset, nslc_all2 = pickle.load(f)

for n in nslc_all2:
	nslc_all.append(n)
for s in station_list2:
	station_list.append(s)
for s in station_list3:
	station_list.append(s)

resp = read_inventory("GB_IRIS.xml", format="STATIONXML")
resp += read_inventory("AM_RASPISHAKE.xml", format="STATIONXML")
resp += read_inventory("GB_BGS_EIDA.xml", format="STATIONXML")

df = pd.DataFrame(data=[])

data = []
for station in station_list:
	df_sta = pd.DataFrame(data=[])
	sta_resp = resp.select(station=station)
	lat = sta_resp[0][0].latitude
	lon = sta_resp[0][0].longitude
	net = sta_resp[0].code
	ppsds = {}
	for day in datelist:
		datestr = day.strftime("%Y-%m-%d")
		print(datestr)
		fn_pattern = "data/{}_*{}*.npz".format(datestr, station)
		for fn in glob(fn_pattern):
			mseedid = fn.replace(".npz", "").split("_")[-1]
			if mseedid not in ppsds:
				ppsds[mseedid] = PPSD.load_npz(fn)#, allow_pickle=True)
			else:
				ppsds[mseedid].add_npz(fn)#, allow_pickle=True)
	displacement_RMS = {}
	for mseedid, ppsd in ppsds.items():
		per = ppsd.period_bin_centers
		displacement_RMS["disp"] = []
		for n_psd, psd in enumerate(ppsd.psd_values):
			RMS = {}
			ix = np.where((per>=1.0/fmax) & (per<=1.0/fmin))
			spec = psd.copy()[ix][::-1]
			f = 1.0/per.copy()[ix][::-1]
			# remove NaNs from the list
			valid = np.where(np.isfinite(spec))[0]
			spec = spec[valid]
			f = f[valid]
			w2f = (2.0 * np.pi * f)
			# The acceleration amplitude spectrum (dB to Power! = divide by 10 and not 20!)
			amp = 10.0**(spec/10.) 

			# velocity spectrum (divide by omega**2)
			vamp = amp / w2f**2

			# displacement spectrum (divide by omega**2)
			damp =  vamp / w2f**2

			RMS = rms(damp, f)
			displacement_RMS["disp"].append(RMS)
		data= pd.DataFrame(data=displacement_RMS, index=pd.DatetimeIndex([d.datetime for d in ppsd.times_processed]), dtype=float)
		rs = data.copy().between_time("5:30", "18:00")
		print(rs)
		rs = rs.resample("1D" ).median().tshift(12, "H")
		rs["lat"] = lat
		rs["lon"] = lon
		rs["net"] = net
		rs["sta"] = station
		df_sta = df_sta.append(rs)
	df_sta_ = pd.concat([df_sta.loc['2020-02-03':'2020-02-06'],
		                 df_sta.loc['2020-02-10':'2020-02-13'],
						 df_sta.loc['2020-02-17':'2020-02-20'],
						 df_sta.loc['2020-02-24':'2020-02-27'],
						 df_sta.loc['2020-03-02':'2020-03-05']])
	try:
		mean_pre_lockdn = df_sta_["disp"].median()
	except:
		continue
	df_sta["disp_pc"] = 100 * (mean_pre_lockdn - df_sta["disp"]) / mean_pre_lockdn
	df = df.append(df_sta)
df.to_csv("BGS_noise.csv")

#		df_tmp = pd.dataframe
#			data.append([station, lat, lon, date, RMS])
#
#df = pd.DataFrame(data=data, columns=["sta", "lat", "lon", "date", "disp"])
#print(df)
