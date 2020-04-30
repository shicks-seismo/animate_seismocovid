#!/usr/bin/env python

import pickle
from obspy import read, read_inventory
import os
from obspy.signal import PPSD
import sys
import multiprocessing
import numpy as np
from functools import partial

def go(nslc_all__):
	print(nslc_all__)
	for nslc in nslc_all__:
		for day in datelist:
			datestr = day.strftime("%Y-%m-%d")
			fn_in = "/DATA/COVID_NOISE/data_{:}/{:}_{:}.mseed".format(dataset, datestr, nslc)
			if day == datelist[-1] :
				continue
			try:
				stall = read(fn_in)
			except:
				continue
			mseedid = [tr.id for tr in stall][0]
			fn_out = "data/{:}_{:}.npz".format(datestr, mseedid)
			if os.path.isfile(fn_out):
				print("%s done already."%fn_out)
				continue
			st = stall.select(id=mseedid)
			st.attach_response(resp)
			ppsd = PPSD(st[0].stats, metadata=resp,
						ppsd_length=1800, overlap=0.5,
						period_smoothing_width_octaves=0.025,
						period_step_octaves=0.0125,
						period_limits=(0.008, 50),
						db_bins=(-200, 20, 0.25))
			ppsd.add(st)
			ppsd.save_npz(fn_out[:-4])
			print(st)
			sys.stdout.flush()
			del st, ppsd
	print("Done")


resp = read_inventory("GB_IRIS.xml", format="STATIONXML")

with open('params_BGS_IRIS.pkl', 'rb') as f:  # Python 3: open(..., 'wb
	datelist, station_list, dataset, nslc_all = pickle.load(f)

nproc = multiprocessing.cpu_count()
nslc_all_splt = np.array_split(nslc_all, nproc)
print(nslc_all)
print(len(nslc_all))
with multiprocessing.Pool(processes=nproc) as pool:
	pool.map(go, nslc_all_splt)




