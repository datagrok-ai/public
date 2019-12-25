<!-- TITLE: Tests: Dialogs -->
<!-- SUBTITLE: -->

# Tests: Dialogs 

Most of the work in platform based on opening of different dialogs. In order to, for example, upload a project, 
apply model or create new column, user needs to open the corresponding dialog.

Dialogs have [input fields](../tests/dialogs-input-fields-test.md) that can be mandatory and optional and have different types of data.

When testing platform dialogs, you should pay attention to following:

* Moving dialog - dialogs can be moved and attached to different parts of the working screen
* Focus - if dialog is in focus, then it has a black border
* Help - in all the dialog headers there is a help button (?). 
If dialog is in focus, then help is shown automatically
* Closing without execution - in all the dialog headers there is a close button (×).
If the dialog is in focus, then you can use ```ESC``` 
* History - history button under which  parameters with which dialogs were execute earlier
* Execute dialog - to execute  dialog there is an "OK" button on it. Same, it can be executed after pressing ```Enter```
If for some reason dialog can not be executed, then "OK" button is gray, ```Enter``` does not work, and the user is shown tooltip with the reason

