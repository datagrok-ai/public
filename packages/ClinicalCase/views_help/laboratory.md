# Laboratory

This view contains several specific charts for analyzing laboratory results.

* [**Hy's law**](https://en.wikipedia.org/wiki/Hy%27s_law)

Scatter plot for analyzing possible risk of a fatal drug-induced liver injury. To generate scatter plot you should select values for ALT, AST and Bilirubin on a context panel.

Shows peak bilirubin values versus peak ALT/AST values across the study. Reference lines are shown at 3*ULN for ALT and AST and 2*ULN for Bilirubin. Possible Hy's law is defined as AST or ALT >= 3*ULN with bilirubin >=2*ULN.

Color corresponds to treatment arm.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/hys_law.gif" height="500" width='800'/>

* **Baseline endpoint**

Scatter plot which shows ratio between laboratory values at some selected baseline and endpoint timepoints.

Scatter plot is divided to 9 parts each of which is annotated with corresponding ration. For instance 'Normal-High' quadrant corresponds to subject who had laboratory value within normal ranges baseline visit but ended up with increased value at endpoint. Thus it can be useful, for example, to identify groups of subjects who developed increasing of some laboratory values compared to baseline or vice versa who started with values out of range but ended within normal ranges.

Baseline, endpoint visits as well as laboratory value can be selected using dropdown lists on a context panel.

Color corresponds to treatment arm.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/bl_ep.PNG" height="500" width='800'/>

* **Laboratory distribution**

This box plot shows distribution of selected laboratory value among all subjects depending on study day.
In particular you can analyze median, min and max values, upper and lower quartiles and detect outliers. Additionally you can evaluate difference between distributions on different study days.

Laboratory value as well as study visit can be selected via dropdown lists on a context panel.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/lab_distr.PNG" height="500" width='800'/>

* **Results**

This tab contains laboratory domain table.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/lab_table.PNG" height="500" width='800'/>
