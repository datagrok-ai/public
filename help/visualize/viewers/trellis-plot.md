<!-- TITLE: Trellis plot -->
<!-- SUBTITLE: -->

# Trellis plot

Trellis Charts are useful for finding the structure and patterns in complex data. A Trellis Chart is a layout of smaller
charts in a grid with consistent scales. Each smaller chart represents rows that belong to a corresponding category. The
grid layout looks similar to a garden trellis, hence the name Trellis Chart.

There are two ways to add a trellis plot, visualized below:

* click on the "Trellis Plot" icon in the toolbox, and then customize the inner chart by clicking on the "gear" icon on
  the left
* create a viewer that you want to eventually become an inner chart, customize it the way you like, and then click
  on `Viewer | Use in Trellis`

![Trellis Plot](../../uploads/gifs/trellis-plot.gif "Trellis Plot")

![viewers-as-trellis](../viewers-as-trellis.gif)

Typically, you want the data split by one or two columns. Use combo boxes on top of the control for that. Note that you
can split data by one column per dimension.

To change the inner viewer type, click on the viewer icon in the left top corner. To edit inner viewer's settings, use
the "gear" icon next to it.

Trellis Plot automatically picks up element renderers for rendering categories. For instance, this is how it looks for
chemical structures after performing [R-Group Analysis](../../domains/chem/r-group-analysis.md):

![R-Group Analysis](../../uploads/chem/r-group-analysis.png "R-Group Analysis")

## Videos

[![Trellis Plot](../../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=1560s)

See also:

* [Viewers](../viewers.md)
* [Table view](../../overview/table-view.md)
* [R-Group analysis](../../domains/chem/r-group-analysis.md)
* [JS API: Trellis plot](https://public.datagrok.ai/js/samples/ui/viewers/types/trellis-plot)
