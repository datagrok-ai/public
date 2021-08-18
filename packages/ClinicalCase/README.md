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
