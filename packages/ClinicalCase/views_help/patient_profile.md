# Patient profile

Patient profile is useful for analyzing events related to particular patient.
You can analyze data from laboratory, adverse events, dug exposure and concomitant medication domains in time and see relations between events. All graphs are linked to the same X axis representing study days and it can be zoomed in and out simultaneously.

Information about events is available in tooltips on mouse hover. For convenience graphs can be collapsed or extended.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/patient_profile_zoom.gif" height="500" width='800'/>

* **Lab values chart**

By clicking on settings button you can choose laboratory values to show on chart. List of available values is extracted from 'lb' domain in provided SDTM data.

Values within normal ranges are colored green, values outside normal ranges are red.

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/patient_profile_lab.gif" height="500" width='800'/>

* **Lab values line chart**

You can also choose laboratory values by clicking settings button.

Laboratory line chart provides the following of calculating values dynamics:

1. Relative changes from baseline
2. Relative values between min and max contained in dataset
3. Relative changes between normalized normal ranges

<img src="https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ClinicalCase/img/patient_profile_lab_line.gif" height="500" width='800'/>
