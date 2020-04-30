#!/usr/bin/env python

from obspy import UTCDateTime 
from obspy.clients.fdsn import Client
import pandas as pd
import os
import pickle

start = UTCDateTime("2020-01-01")
end = UTCDateTime()
data_provider = "IRIS"
network = "GB"
dataset = "BGS_IRIS"

c = Client(data_provider)
inv = c.get_stations(network=network, starttime=start, endtime=end, level="response")
inv.write("{:}_{:}.xml".format(network, data_provider), format="STATIONXML")

station_list = [sta.code for net in inv for sta in net]
datelist = pd.date_range(start.datetime, end.datetime, freq="D")
nslc_all = []


for station in station_list:
	inv_sta = inv.select(station=station, channel="*HZ")
	channel = inv_sta[0][0][0].code
	location = inv_sta[0][0][0].location_code
	nslc = "{}.{}.{}.{}".format(network, station, location, channel)
	nslc_all.append(nslc)

	for day in datelist:
		datestr = day.strftime("%Y-%m-%d")
		fn = "/DATA/COVID_NOISE/data_{:}/{:}_{:}.mseed".format(dataset, datestr, nslc)
		fn2 = "data/{:}_{:}.{:}.{:}.{:}.npz".format(datestr, network, station, location, channel)
		if day == datelist[-1]:
                        print("Not getting today's data")
                        continue
		if os.path.isfile(fn):
                        print("Already have mseed")
                        continue
		if os.path.isfile(fn2):
                        print("Already have npz")
                        continue
		try:
			st = c.get_waveforms(
					network, station, location, channel,
					UTCDateTime(day)-1801, UTCDateTime(day)+86400+1801,
					attach_response=True)
		except:
			continue
		print(st)
		st.write(fn)

with open('params_{:}.pkl'.format(dataset), 'wb') as f:  # Python 3: open(..., 'wb')
	    pickle.dump([datelist, station_list, dataset, nslc_all], f)



