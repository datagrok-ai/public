### Visits

View contains summary information about all the visits performed in frames of clinical trial protocol (unscheduled visits are NOT included).
Information is represented in two ways:
* **Grid**
Columns represent visits, rows correspond to patients. In each cell you will see several numbers of different colors. Colors correspond to domains. On the ribbon panel you can see the names of selected domains. For example 'lb' domain is colored blue.
The numbers themselves mean the following:
- lb and vs - number of values evaluated at the visit for selected subject
- ex - number of investigational product distributions to selected subject at the visit
- ae - number of new adverse events occurred since last visit
- cm - number of new concomitant medication subject took since last visit
The domains can be selected/removed using `Domains` button on the ribbon panel.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/visits_grid_1.gif" height="500" width='800'/>

When selecting the cell you will see more detailed information about patient visit in a context panel. There is an expandable panel for each domain which contains rows from corresponding domain table.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/visits_grid_2.gif" height="500" width='800'/>

When selecting the column the context panel will show study visit summary such as total number of patients, min and max visit dates and some charts representing statistics for different domains:
- pie chart for investigational product exposure
- bar charts for new adverse events and concomitant medications since last visit
- box plots with values distributions for laboratory and vital signs domains

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/visits_grid_3.gif" height="500" width='800'/>

* **Heat maps**
Intensity of color in a heat map represents number of events. Domain can be selected in a drop down list on a ribbon panel.
Map can be sorted by treatment arm. 
Selecting cell or column works as in grid regimen - corresponding context panel will appear. 

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/visits_heatmap.gif" height="500" width='800'/>