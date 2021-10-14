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

* Lab values chart

By clicking on settings button you can choose laboratory values to show on chart. List of available values is extraced from 'lb' domain in provided SDTM data.

Values within normal ranges are colored green, values outside normal ranges are red.

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/patient_profile_lab.gif" height="500" width='800'/>

* Lab values line chart

You can also choose laboratory values by clicking settings button.

Laboratory line chart provides two ways of calculating values dynamics:
- relative changes from baseline
- relative values between min and max contained in dataset

<img src="https://github.com/datagrok-ai/public/blob/clinical-case-app/packages/ClinicalCase/img/patient_profile_lab_line.gif" height="500" width='800'/>








