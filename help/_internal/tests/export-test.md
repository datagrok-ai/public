<!-- TITLE: Tests: Export -->
<!-- SUBTITLE: -->

# Tests: Export

## Testing scenario

Platform provides feature for saving data in various formats (*csv*, *json*, *xlsx*). Also possible to save entire
project in *d42* format or in ZIP-archive that will contain *csv* files of all project tables and *PNG* files

1. Open *demog.csv* file

1. Export *demog* table as *csv* file. ( **File** | **Save** | **Table as CSV**)

* File downloaded to local storage
* Downloaded file is called *"demog.csv"* (table name)

1. Open downloaded file *"demog.csv"*

* File opened in platform fully (data in cells correspond to the original table)
* All column names, tables correspond to original table
* Column types correspond to original table

1. Close table from exported file (**Workspace | Scratchpad | Tables | Close**)

1. Open **File** | **Save** | **Table as CSV (options)** for export settings to *csv* format (you need to export
   original table *"demog"*).

* *"Save as CSV"* dialog is open

1. Change the values ​​of all boolean fields to **FALSE**.(*Include header*, *Selected columns only*
   , *Filtered rows only*, *Selected rows only*)

1. Export *"demog"* table to CSV files with different delimiters choosing appropriate in **
   Delimiters** field. (*Comma*, *Tab*, *Semicolon*, *Pipe*)

* Files downloaded to local storage

1. Open downloaded files in any suitable text editor

* Delimiters in *CSV* files match the delimiters selected during export

1. Repeat **3-4** steps for each open file

1. Open *"Save as CSV"* dialog (**File** | **Save** | **Table as CSV (options)**) for export settings to *csv* format (
   you need export original table *"demog"*).

1. Export table with different values of **New line** field (*/r/n*, */n*). All Boolean fields must be **FALSE**

1. Open downloaded files in any suitable boolean text editor

* Delimiters in *CSV* files match the new line symbol selected during export

1. Repeat **3-4** steps for each open file

1. Apply filter on value *"M"* to **"SEX"** column

1. Select all rows *​​"Asian"* in **"RACE"** column and select **"SEX"**, **"RACE"** and **"
   Height"** columns (for original *demog* table)

1. Open *"Save as CSV"* dialog (**File** | **Save** | **Table as CSV (options)**) for export settings to *csv* format (
   you need export original table *"demog"*).

1. Type *"test"* into **Missing Value** field

1. Set values ​​for boolean fields: **"Include header**" = FALSE,  **"Selected columns only"** = TRUE, **"Filtered rows
   only"** = TRUE, **"Selected rows only"** = TRUE

1. Click on ```"Set as default"``` button

1. Export table in CSV with given parameters

1. Open downloaded file in any suitable text editor

* There are three columns in downloaded *CSV* file. (**"SEX"**, **"RACE"** and **"Height"**)
* The file has only rows that correspond to values *​​"Asian"* in **"RACE"** column and *"M"*
  in **"SEX"** column
* Column headers are missing in file
* Null values ​​in **Height** column replaced by *"test"*

1. Repeat **3-4** steps

1. Export original table as *csv* file. ( **File** | **Save** | **Table as CSV**)

* File downloaded to local storage
* Downloaded file is called *"demog.csv"* (table name)
* Downloaded data correspond to settings in step #18

1. Export *demog* table as *JSON* file. ( **File** | **Save** | **Table as JSON**)

1. Repeat **3-4** steps

1. Export *demog* table as *MS Excel* file. ( **File** | **Save** | **As Excel Document**)

1. Repeat **3-4** steps for each open file

1. Save *demog* table as presentation file. ( **File** | **Save** | **As Presentation**)

1. Open saved *demog.pptx* file

* On presentation slide *"demog"* table is displayed

1. Open one more table in platform. (for example *"Products"*)

1. Make sure that when exporting to different file formats (*CSV*, *JSON*, *Excel*), current table is saved.

1. Save tables as presentation file. (**File** | **Save** | **As Presentation**)

1. Open saved file

* On presentation slides all tables is displayed
* File is called *"presentation.pptx"*

1. Save open tables as project (**File** | **Save** | **Tables as Project**)

* Downloaded file is named *"project.d42"*.

1. Repeat **3-4** steps

1. Save open tables as *ZIP* (**File** | **Save** | **CSV and PNG as ZIP**)

* Downloaded file is named *"project.zip"*

1. Open downloaded *"project.zip"*

* Inside archive there are *CSV* and *PNG* files for all tables.
* *PNG* files are saved for each created viewers on table layout.

1. Add new column to *"demog"* table (corresponding icon on toolbar). For new column use math formulas and
   constants ([Add new column](../../transform/add-new-column.md))

1. Check correctness of export with new column by repeating steps **2-35**.

1. Add new viewer to table layout ([Viewers](../../visualize/viewers.md))

1. Check correctness of export with added viewer by repeating steps **2-35**.

1. Apply filter for table ([Filters](../../visualize/viewers/filters.md))

1. Check correctness of export with filtering by repeating steps **2-35**.

1. Change date format in **STARTED** column

1. Check correctness of export with different date formats by repeating steps **2-35**.

1. Rename table and its columns

1. Check correctness of export with changing names by repeating steps **2-35**.
