<!-- TITLE: Use Cases: Selecting franchise locations -->
<!-- SUBTITLE: -->

# Use cases: Selecting franchise locations

Owner: Vasiliy

Goal: to show workflow how to create project; get data from external datasource and preprocess it; train and use model;
analyse basic components of input data; create and share report.

Features: creating project, accessing data, creating query, post-processing, Jupyter, scripting, predictive modeling

1. Register new User (Guest | Sign Up, enter email and password)
2. Create new project (Ctrl+Q, "Project View")
3. Get data from external source:

* Open "Connect to Data" view (File | Connect to Data...)
* Create new connection (Context menu on PostgresSQL or Postgres | Add connection
    * Name: Starbucks
    * server: reddata.org
    * db: starbucks
    * login: TBD
    * password: TBD
* Create new Data Query (Context menu on "Starbucks" Data Connection | Add query)
    * Name: Get Starbucks US
    * Query: select * from starbucks_us
    * To get data click: "RUN" of (Context menu on "Get Starbucks US" Data Query | Run)

4. Select only NY state rows (State | Column category button "≡" | Select "NY") and open as separate table
   (Property panel | Extract rows (or Keep rows))
5. Create derived column: "${street_address}, ${city}, ${state}".

* Open Columns View (View | Columns)
* Open "Add new column" dialog (on Toolbar)
* Drag-and-Drop corresponding columns and separate with commas
* OK

6. Convert new Column to coordinates

* Rename Column "${street_address}, ${city}, ${state}" to "Address"
    * Open context menu on column header
    * Click "Rename"
    * Rename to "Address"
* Open context menu on column header
* Click "Address to coordinates"

7. Coordinates to statistics:

* Click on table view header
* Execute: Property panel | Algorithms | Coordinates to statistics

8. Multivariate analysis:

* Open PLS dialog (Tools | Data Science | Multivariate Analysis (PLS)...)
* Select table "train" in Table selector
* Select features: around 10 from statistics including "population_density", "
  us_population_low_income", "us_housing_units_one_person"
* Select outcome row "Sales"

9. t-Test of selected feature using Scripting:

* Open Console (Tools | Console), needed to see result of equation
* Create new script (Tools | Scripting | New Script...)
* Click "Erase script" button "⬒"
* Paste following code:

  ```
  #name: t-test
  #description: Welch's t-test
  #language: r
  #input: dataframe data [input data table]
  #input: column x {type:numerical} [x axis column name]
  #input: column y {type:numerical} [y axis column name]
  #output: double pvalue {action:show} [p-value of t-statistics]
  require(stats)
  ttest = t.test(data[[x]], data[[y]])
  pValue = ttest$p.value
  ```

* Run script (Click "▶" button)
* Select "Sales" as "x" any of features as "y"
* Click "OK"
* Repeat for all features
* Use "History" button to get previous values and change only "y"

10. Create model:

* Open "Predictive Modelling" dialog (Tools | Predictive modeling | Train)
* Select "Method": "Distributed Random Forest"
* Select "Features": "population_density", "us_population_low_income", "us_housing_units_one_person"
    + some more
* Select "Predict": "Sales"
* Click "Train" button

11. Fill "research" table

* Open Google Map Viewer (Navigation Panel | Google Map or Add | Geo | Google Map)
* Select point of interest using double click
* Open selected coordinates as table (Context menu of viewer "≡" | Save selection as table)

12. Apply model on "research" table

* Convert Coordinates to statistics
* Open "Predictive Modeling Browser" (Tools | Predictive modeling | Browse Models)
* Search for trained model (applicable model will be at the top of list, or use search by defined previously name)

13. Open Google maps view to visualize results on map

* Click on "research" table
* Open "Google Map" view (Click on Navigation Panel "Google Map" or Add | Geo | Google Map)

14. Jupyter Notebooks report

* Open Jupyter Notebooks editor (Tools | Jupyter Notebooks)
* Create report using Script language

15. Share created Project

* Open "Project View" (Ctrl+Q, "Project View")
* Open context menu on project that should be shared
* Click "Share"

Sales formula:
log(${population_density} + 1) * (${us_population_eighteen_to_twenty_four_years_old} + ${us_population_low_income} +
${us_housing_units_one_person} + ${us_population_bachelors_degree}) / 4 * 100000
