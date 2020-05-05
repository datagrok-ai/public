This [package](https://datagrok.ai/help/develop/develop#packages) stores a collection of sample [scripts](https://datagrok.ai/help/compute/scripting) written using the internal language - [GrokScript](https://datagrok.ai/help/overview/grok-script)

[Grok script](https://datagrok.ai/help/overview/grok-script) language is used to control or automate everything within the [Datagrok](datagrok.ai) platform. 
Use it to transform data, automate workflows, run queries, evaluate numerical expressions, execute commands, record macros, perform statistical computations, execute other scripts, etc.


# Info-panel

Separately, here are samples of using [Grok script](https://datagrok.ai/help/overview/grok-script) to create [Info Panels](https://datagrok.ai/help/discover/info-panels).

[Info Panels](https://datagrok.ai/help/discover/info-panels) provide additional information about the current context (which can be a table, a column, or pretty much any other [object](https://datagrok.ai/help/overview/objects) in [Datagrok](datagrok.ai) platform). 
Info panels are meant to be easily developed by the users of the platform, and shared with other users. You can use all features of the Grok platform, such as scripting, data querying and transformation pipeline, user-defined functions, markup, viewers, predictive models.

[Info Panels](https://datagrok.ai/help/discover/info-panels) is displayed as a separate tab in the [Property Panel](https://datagrok.ai/help/overview/property-panel).
[Property Panel](https://datagrok.ai/help/overview/property-panel) is part of [Datagrok](datagrok.ai) which displays the properties of the current [object](https://datagrok.ai/help/overview/objects) in [Datagrok](datagrok.ai). You can read more about it on our [Wiki](https://datagrok.ai/help/overview/property-panel).

# Example
```
#name: Scatter Plot
#description: Panel that contains an interactive scatter plot
#language: grok
#tags: panel
#input: dataframe table
#condition: table.name == "demog" && table.columns.containsAll(["height", "weight", "age", "sex"])
#output: viewer plot

plot = table.ScatterPlot("height", "weight", "age", "sex")
plot.showRegressionLine = true
```

This example is an [Info Panel](https://datagrok.ai/help/discover/info-panels) written in a [Grok script](https://datagrok.ai/help/overview/grok-script).
Definition of [Info Panel](https://datagrok.ai/help/discover/info-panels) is determined by the special tag ```#panel``` in script header.

User will see this [panel](https://datagrok.ai/help/discover/info-panels), which shows a [scatter plot](https://datagrok.ai/help/visualize/viewers/scatter-plot), after opening a table with the name *"demog"*, in which there will be columns *height*, *weight*, *age* and *sex*.
[Panel](https://datagrok.ai/help/discover/info-panels) display conditions are defined in the [script](https://datagrok.ai/help/overview/grok-script) header.
