#name: GetStations
#language: python
#environment: ulmo
#input: string country
#output: dataframe stations

import ulmo
import pandas

stations = ulmo.ncdc.ghcn_daily.get_stations(country=country, as_dataframe=True)