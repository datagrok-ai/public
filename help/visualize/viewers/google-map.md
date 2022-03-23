<!-- TITLE: Google Map viewer -->
<!-- SUBTITLE: -->

# Google Map viewer

Google Map Viewer overlays latitude/longitude data from the corresponding table on top of the Google Map.

By default, the map would pick up columns named "lon" / "lat", or "longitude" / "latitude". Open map settings to specify
two columns that represent longitude and latitude. To control marker settings, click on the hamburger menu and select "
Marker settings".

In case the dataset contains geographical data (such as addresses), but does not contain longitude and latitude
coordinates, you might want to extract coordinates using one of
the [geographical functions](https://public.datagrok.ai/functions?q=%23geo), such as #{cmd(AddressToCoordinates)}.

![Google Map](../../uploads/viewers/google-map.png "Google Map")

![Map big data](google-map-city-perf.gif "Map big data")

## Videos

[![Google Map Video](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=3392s)

See also:

* \#{cmd(AddressToCoordinates)}
* \#{cmd(IpToCoordinates)}
* \#{cmd(AddressToCoordinates)}
* [Viewers](../viewers.md)
* [Table View](../../overview/table-view.md)
* [JS API: Google Map](https://public.datagrok.ai/js/samples/ui/viewers/types/google-map)
