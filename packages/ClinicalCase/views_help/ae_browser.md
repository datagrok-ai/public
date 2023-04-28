# AE browser

This view is useful for exploring adverse events in details including preceding events from domains other then ‘ae’.

The view is basically 'ae' domain table which contains list of all registered adverse events. The table can be filtered.

By selecting row in a table the following information will be displayed on a context panel:

- subject ID (in a tooltip on mouse hover you will see basic demographic characteristics - age, sex, race, treatment arm)
- AE name preceded by AE severity
- Days of study in which AE occurred (in a tooltip on mouse hover you will see real AE dates)
- input with number of days before AE for which you want to analyze events in other domains (by default it's 5)
- list of expandable domain panels which contains rows with events occurred during selected period before the AE (by default adverse event, drug exposure and concomitant domains are selected, but you can add other domains by clicking on `+` button)

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/ae_browser.gif" height="500" width='800'/>
