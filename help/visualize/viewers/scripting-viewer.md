<!-- TITLE: Scripting Viewer -->
<!-- SUBTITLE: -->

# Scripting Viewer

Scripting Viewers are viewers implemented in R, Python or Julia. Internally, they
use [Scripting](../../compute/scripting.md) for integration with the Datagrok platform.

While not as interactive as the core Datagrok [viewers](../viewers.md), they allow to easily use 
thousands of visualizations already developed for these languages.

To **add an existing viewer** to a table view, select it from the
 **Add | Scripting Viewers** menu.

To edit the rest of the properties, either click on the "gear" icon on top of the viewer,
or press F4 when the viewer has focus, or right-click and select `Viewer | Propeties`.

## Customize Scripting Viewer script

Since any of Scripting Viewers is [Scripting](../../compute/scripting.md), it is easy to customise 
existing or create your own viewer. To customise script code, open script code by clicking
on the popup menu's **Edit script**. To add new Scripting Viewer to the main menu, add "viewers" 
tag to script header.

## Scripting Viewer code example

Example shows code of simple scatterplot written on R, using "ggplot2" library.

```
#name: Scatter Plot
#language: r
#tags: demo, viewers
#input: dataframe t
#input: column xColumnName
#input: column yColumnName
#input: column colorColumnName
#output: graphics

require(ggplot2)

# Compose input columns into data frame with required names
data <- data.frame(x=t[[xColumnName]], y=t[[yColumnName]], color=t[[colorColumnName]])

# Plots
plotScatter <- ggplot(data, 
  aes(x, y, colour=color), xlab=x, ylab=y) +
  labs(x=xColumnName, y=yColumnName) +
  geom_point()
print(plotScatter)
```

## Demo project

Open #{x.demo:TimeSeriesDecomposition."Time series decomposition"} project as an example of 
Scripting Viewer usage for Time series decomposition. 

[Short video lesson](https://drive.google.com/uc?export=download&id=17D-X_5_wPJGeWd9Oc_ZfU48tOvYwf1Do)

### Videos

<iframe width="560" height="315" src="https://www.youtube.com/embed/jHRpOnhBAz4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

See also: 
  
* [Scripting](../../compute/scripting.md)
* [Viewers](../viewers.md)
