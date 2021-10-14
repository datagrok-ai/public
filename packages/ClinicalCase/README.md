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


