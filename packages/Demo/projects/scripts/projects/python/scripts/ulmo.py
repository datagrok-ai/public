#name: Ulmo
#language: python
#environment: ulmo
#output: dataframe df

from ulmo.nasa import daymet
ornl_lat, ornl_long = 35.9313167, -84.3104124
df = daymet.get_daymet_singlepixel(
  longitude=ornl_long, latitude=ornl_lat, 
  years=[2012, 2013])