# ClinicalCase

Provides support for dealing with clinical data represented in the
[SDTM](https://www.cdisc.org/standards/foundational/sdtm) format.

Issue tracker: https://github.com/datagrok-ai/public/milestone/1

The project objectives include:

* Working with either databases or file folders 
    * Automatic recognition of SDTM data
* Explaining attribute names
    * Domains: `CV` -> `Cardiovascular System Findings`
    * Attributes: `AEMODIFY` -> `Modified Reported Term`
* Content validation
    * Completeness (columns with required variables are present)
    * Data types match
    * Values comply with controlled terminology
    * Out-of-range values
* Patient selection criteria
    * Global filter with conditions based on the values in any table 
* Clinical data-specific visualizations
    * Timelines
    * Hy's Law
* Pre-defined views
    * Study Summary
      * Dates, cohorts
    * Patient Profile
  * Lab Results View
    * Pre-defined, configurable normal ranges
* Defining groups (treatment / control)
    * Automatic discovery of statistically significant changes between groups   
* Subject profile view (in progress)
    * Swim lanes
* Miscellaneous
    * Automatic calculation of study days, if necessary


* `Timelines` viewer for visualizing a flow of events (see the [community forum discussion](https://community.datagrok.ai/t/visualization-related-updates/521/4))


![](timelines.gif)

Views
===============

### Summary

View contains basic information about the study and subject population.
Here you can find number total number of subjects and sites along with the cumulative enrollment line chart which is useful to assess the enrollment dynamics throughout the study.

'Errors' section contains number of errors revealed by validation process in each domain. Validation of SDTM tables (domains) is performed once at the start of application. By clicking on errors number you will be redirected to Validation view.

Also there is a couple of charts with essential population characteristics – age, race, sex, treatment arm – which can be used to assess their distribution within population. Charts are interactive and linked with each other. For instance, by clicking on 'M' sex the other charts will show the proportion of male subjects.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/summary.gif" height="500" width='800'/>

In case trial is registered on [clinicaltrials.gov](https://clinicaltrials.gov/) property panel on the right will basic study information extracted from database along with the link to the study on [clinicaltrials.gov](https://clinicaltrials.gov/).

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/Summary_property_panel.PNG"/>


### Timelines

Timelines view allows to visualize events in time. X axis is time axis reflecting study days. On Y axis there are subjects.
Events are shown either as a point(in case event lased for one day) or as a line (for cases events were prolonged in time)

Adverse events, investigationsl drug exposure and concomitant medication domains are available for analysis(in case SDTM data contains corresponding tables). Filters can also be applyed for each domain.

By zooming in and out you can drill down to particular patient and event or otherwise see the picture of events in general.

Information about particular event is shown in tooltip on mouse hover.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/timelines.gif" height="500" width='800'/>

Several domains can be shown simultaneously on the graph. For instance, the following screenshot shows severe general and cardiac disorders VS aspirin intake.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/timelines.PNG" height="500" width='800'/>


### Patient profile

Patient profile is useful for analysing events related to particular patient.
You can analyse data from laboratory, adverse events, dug exposue and concomitant medication domains in time and see relations between events. All graphs are linked to the same X axis representing study days and it can be zoomed in and out simultaneously. 

Information about events in available in tooltips on mouse hover. For convenience domains han be collapsed or extended.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/patient_profile_zoom.gif" height="500" width='800'/>

* **Lab values chart**

By clicking on settings button you can choose laboratory values to show on chart. List of available values is extraced from 'lb' domain in provided SDTM data.

Values within normal ranges are colored green, values outside normal ranges are red.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/patient_profile_lab.gif" height="500" width='800'/>

* **Lab values line chart**

You can also choose laboratory values by clicking settings button.

Laboratory line chart provides two ways of calculating values dynamics:
1. Relative changes from baseline
2. Relative values between min and max contained in dataset

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/patient_profile_lab_line.gif" height="500" width='800'/>

### Adverse events

View to explore Advere events domain. Contains couple of charts along with Adverse events table view.

* **All events**

Plots all adverse events registered during the study. Detailed information is available in tooltip. Also when selecting Adverse event on the plot corresponding row becomes current in Adverse event table. So Adverse event of inteest can easily be analyzedin details. Selection also works vice versa - when selecting a row in a table corresponding Adverse event is selected on a scatter plot.

Scatter plot can be zoomed in and out to drill down to particular patient or see picure in general. 

Color indicates severity of an Adverse event.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/ae_all_events.gif" height="500" width='800'/>


* **Events per week**

Historgam with number of Adverse events occurred per week throughout the study. Can help to analyze overall dynamics of Adverse events occurrence.

Color corresponds to severity of an Adverse event.

* **Barcharts**

There are barcharts plotting number of events by Type, Body system, Causality and Outcome.

Color corresponds to treatment arm. So you can visually assess proportion of each treatment group.

Barcharts are also interactive. By selecting one of the groups in a barchart (for instance, adverse events related to study drug for patients who took placebo) corresponding rows will be selected in the table and 'All events' scatter plot.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/ae_barcharts.gif" height="500" width='800'/>


### Laboratory

This view contains several specific charts for analyzing laboratory results.

* [**Hy's law**](https://en.wikipedia.org/wiki/Hy%27s_law)

Scatter plot for analyzing possible risk of a fatal drug-induced liver injury.

Shows peak bilirubin values versus peak ALT/AST values across the study. Reference lines are shown at 3*ULN for ALT and AST and 2*ULN for BILI and ALP. Possible Hy's law is defined as AST or ALT >= 3*ULN with bilirubin >=2*ULN.

Color coresponds to treatment arm.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/hys_law.PNG" height="500" width='800'/>

* **Baseline endpoint**

Scatter plot which shows ratio between laboratory values at some selected baselibe and enpoint timepoints.

Scatter plot is divided to 9 parts each of which is annotated with corresponding ration. For instance 'Normal-High' quadrant corresponds to subject who had laboratory value within normal ranges baseline visit but ended up with increased value at endpoint. Thus it can be useful, for example, to identify groups of subjects who developed increasing of some laboratory values compared to baseline or vice versa who started with values out of range but ended within normal ranges.

Baseline, enpoint visits as well as laboratory value can be selected using dropdown lists above the scatter plot.

Color coresponds to treatment arm.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/bl_ep.PNG" height="500" width='800'/>

* **Laboratoty distribution**

This box plot shows distribution of selected laboratory value among all subjects depending on study day. 
In particular you can analyze median, min and max values, upper and lower quartiles and detect outliers. Additionally you can evaluate difference between distributions on different study days. 

Laboratory value as well as study visit can be selected via dropdown lists above box plots.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/lab_distr.PNG" height="500" width='800'/>

* **Results**

This tab contains laboratory domain table.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/lab_table.PNG" height="500" width='800'/>

### Biomarkers distribution

The view shows distribution of selected biomarker values for all subjects at selected visit splitted by selected parameter. It is useful to evaluate median, upper and lower quartiles and outliers. And since distributions in each boxplot can be grouped by certain parameter it is also possible to evaluate difference in biomarker values distribution between different groups.

By default 4 biomarkers with min p-value at the earliest study visit splitted by treatment arm are shown. Baseline visit, split parameter as well as biomarkers can be further changed via dropdown lists and settings button above the boxplots.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/biomarkers_dist.gif" height="500" width='800'/>








