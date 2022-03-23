<!-- TITLE: Parallel coordinates plot -->
<!-- SUBTITLE: -->

# Parallel coordinates plot

Parallel coordinates is a common way of visualizing high-dimensional geometry and analyzing multivariate data. To show a
set of points in an n-dimensional space, a backdrop is drawn consisting of n parallel lines, typically vertical and
equally spaced. A point in n-dimensional space is represented as a polyline with vertices on the parallel axes; the
position of the vertex on the i-th axis corresponds to the i-th coordinate of the point. This visualization is closely
related to time series visualization, except that it is applied to data where the axes do not correspond to points in
time, and therefore do not have a natural order. Therefore, different axis arrangements may be of interest.

To change columns, set "Column Names" via the property panel. To rearrange columns, drag column name into the desired
location.

General:

|                    |                   |
|--------------------|-------------------|
| Right click        | Context menu      |
| Drag column name   | Rearrange columns |
| Drag column filter | Filter data       |

![PC Plot](../../uploads/gifs/pc-plot.gif "PC Plot")

## Videos

[![PC Plot](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=1798s)

See also:

* [Viewers](../viewers.md)
* [Table view](../../overview/table-view.md)
* [JS API: PC plot](https://public.datagrok.ai/js/samples/ui/viewers/types/pc-plot)
