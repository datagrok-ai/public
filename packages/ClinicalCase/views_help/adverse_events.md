# Adverse events

View to explore Adverse events domain. Contains couple of charts along with Adverse events table view.

* **All events**

Plots all adverse events registered during the study. Detailed information is available in tooltip. Also when selecting Adverse event on the plot corresponding row becomes current in Adverse event table. So Adverse event of interest can easily be analyzed in details. Selection also works vice versa - when selecting a row in a table corresponding Adverse event is selected on a scatter plot.

Scatter plot can be zoomed in and out to drill down to particular patient or see picture in general.

Color indicates severity of an Adverse event.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/ae_all_events.gif" height="500" width='800'/>

* **Events per week**

Histogram with number of Adverse events occurred per week throughout the study. Can help to analyze overall dynamics of Adverse events occurrence.

Color corresponds to severity of an Adverse event.

* **Barcharts**

There are barcharts plotting events distribution by Type, Body system, Causality and Outcome.

Color corresponds to treatment arm. So you can visually assess proportion of each treatment group.

Barcharts are also interactive. By selecting one of the groups in a barchart (for instance, adverse events related to study drug for patients who took placebo) corresponding rows will be selected in the table and 'All events' scatter plot.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/ae_barcharts.gif" height="500" width='800'/>
