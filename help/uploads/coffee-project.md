# Coffee Place Project

## [Introduction](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=0s)

Welcome to the Coffee Place Project! Right now, there is sales data in front of your eyes, and based on it, we will be selecting a new spot for the coffee franchise. To repeat these steps on your data, check out [these instructions](https://datagrok.ai/help/access/data-connection), or watch [this video](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=23s). Now, let's get started!

## [Data Exploration](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=60s)

First, take your time to explore the data. Hover your mouse over the table name to see general information.

**Hint:** table `us_locations_with_sales` contains `12` columns and `13,509` rows.

Go over the column names to find out which categories the data falls into and check if the columns have missing values. All of this is at your fingertips—just hover the mouse! For values in the cells, do the same trick again, and you will get the preview of the whole data record.

Now, click on the column names `Street Address`, `City`, `State/Province`, and `Sales` while holding the `Shift` key. This will select the columns. `Extract` them with the corresponding menu command.

Alternatively, as we are only interested in a few columns, we can safely get rid of the rest. Press `Shift` and click on the headers `Brand`, `Store Number`, `Store Name`, `Ownership Type`, `Country`, `Postcode`, `Phone Number`, and `Timezone`. Once you have selected these columns, click on ![the second icon from the right at the top menu](pic-url) or `Shift + Del` to delete them. Selected the wrong column? No worries, use `Ctrl + Shift` to deselect it or undo selection by pressing the `Esc` key.

Next, let's restrict our analysis only to one state, let it be `NY`. To do that, use [the hamburger menu to the right of the column name](pic-url). Find `NY` in the list and click on the name of the state. In the right pane, you will see actions available for `645` records filtered by location. Go for the `Keep rows` option to work with them.

We can create a column containing the full address of the coffee place following these steps:

- choose [the rightmost icon of the menu](pic-url) to add a new column
- type in its name, for example, `Address`
- simply drag the columns into the field preserving the order: `Street Address`, `City`, `State/Province`
- format the record with commas or just paste the following: `${Street Address}, ${City}, ${State/Province}`

Voilà! As we change the formula, we get a preview of the results on the fly, which are refreshed interactively. Save your work and proceed to the next section.

## [Geocoding and Visualization](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=233s)

Now that you have the `Address` column, you might wonder what to do with it. Fortunately, the platform recognizes addresses in the column. Click on its name to see the suggested action `Address to coordinates`. It converts a given address to `latitude` and `longitude`, which is called **geocoding**. After applying this function, two corresponding columns should appear in your table.

Check if you managed to locate all of the addresses. Again, just hover the mouse over the column headers. If you see `null` values in there, find a blank cell, click on it and press `Shift + Enter`. This way, you can select rows with similar values (`nulls`, in this case) and work with them. For analysis purposes, it is better to delete these rows.

Finally, let's visualize the data! Open the `Google Map` viewer in the left pane. On the map, you can see coffeehouse locations in New York State. To start with, let's color-code sales so that a warmer color would represent a larger value, and vice versa. Find `Marker settings` in the hamburger menu of your `Map` viewer and select column `Sales` for marker `Color`. This sort of visual representation sheds light on general trends, e.g., coffee shops in Manhattan did remarkably well. Yet, based on that alone, it is difficult to predict sales. So our next step is to map demographic data to the given coordinates.

## [Mapping U.S. Census Data to Coordinates](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=371s)

Now that the table contains `latitude` and `longitude`, the platform suggests converting `Coordinates To Statistics`. Click on the table name, find this function in the `Algorithms` section on the right, and give it a try. The conversion will take some time, normally, less than a minute. Once the data is retrieved, explore new columns by browsing the `Columns` list under the `Viewers` section on the left. The data includes population density, number of housing units, landscapes and elevation, age and education, etc. Just to check, point your cursor on the table name: you should have `44` columns.

## [Missing Values Imputation](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=506s)

To make use of the diverse data, we can apply a technique called **multivariate analysis**. Go to the top menu and select `ML > Multivariate Analysis (PLS)`. We are going to predict `sales` based on all numerical values, except for `latitude` and `longitude`. To do that, open `Features` and select `All`, then uncheck the first two columns (`latitude`, `longitude`). As the platform warns you, there are missing values in the dataset. First, we should deal with it by either removing the rows or imputing these values. As this section name suggests, we decide on the second option.

Select `Missing Values Imputation` from the end of the system warning. This should open a new window. Then you can open `Impute` and sort the columns by the number of empty values (double-click on `nulls`). Select all columns where the number of nulls is greater than `0` (around `31` columns). Concerning the `Data` used for imputation, select all columns, except for `Street Address`, `City`, `State/Province`, `Address`, `latitude`, `longitude`, and `Sales` (`37` columns in total). Then impute the data using the `Nearest Neighbors` algorithm with a default value of `5`. Simply put, a missing value is computed as the average of the corresponding values of `5` records closest to it. At last, we no longer have missing values in the dataset.

## [Multivariate Analysis](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=679s)

Coming back to the multivariate analysis, repeat the steps we performed earlier:

- In the top menu, select `ML > Multivariate Analysis (PLS)`
- Choose to predict values for the `Sales` column
- Select all `Features` and uncheck `latitude` and `longitude` (`36` columns in total)

You may use the preset number of components equal to `3`. At the bottom of the dialogue, check all the charts you want to render: 

- `Scores` (a scatter plot that shows correlation between observations)
- `Explained Variance` (a bar chart with variable variance explained)
- `Correlation Loadings` (a scatter plot that shows correlation between variables)
- `Predicted vs Reference` (a scatter plot with the regression line comparing predicted vs reference outcomes)
- `Regression Coefficients` (a bar chart with regression coefficients)

Hopefully, now you see a bunch of visualizations. We suggest that you look into these findings to figure out which features correlate with each other and which of them have a greater impact on the outcome. We would especially like to draw your attention to the scatter plot `Predicted vs Reference`, which shows actual predictions for sales. Later, you will see how well a predictive model can cope with the same task. Also, you can learn more about multivariate analysis from [our Wiki page](https://datagrok.ai/help/explore/multivariate-analysis/pls).

Later in this example, we will use the following features:

- `us_housing_units_no_vehicle`
- `us_population_twenty_five_to_sixty_four_years_old`
- `population_density`
- `us_population_bachelors_degree`

After **feature selection**, it's time to train your first predictive model!

## [Predictive Modeling](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=779s)

To train a model, go to the `Models` section in the left pane and click on the `Train` button. Above it, sometimes you can see other models trained on the same data. This comes in handy if you want to compare them with your model's performance. Besides, you can manage the available models right from this panel (should you choose to run a model, edit it, or share).

Now, let's proceed to the `Predictive Modeling` window. You will see a number of fields, feel free to tweak all the parameters you need. Here are some basic steps:

- Select `Sales` as the column you want to predict values for
- Select features you want to include into the analysis, e.g., the ones provided above
- Change the method from `Auto ML` to the one you find applicable, e.g., `Distributed Random Forest`

Notice that the model's name is updated along with your changes. This will help you find your model later. And you can name it differently if you like.

Now, return to the table and go to `ML > Browse Models`. You should see a list of models where you can find the one you have just created. If you want to compare it with some of the existing ones, select the models while holding `Ctrl` and click on the `Compare` command. If you want to apply your model to the data, right-click on its name and select the table.

Additionally, there are two more ways to run a model from an open table:

- In the left pane, browse the `Models` section to find a model, then right-click and select `Run`
- In the `ML` section of the top menu, click on the `Apply model` option

So, assuming you have applied your model in one of the mentioned ways, now your table contains a column called `outcome`. Let's see how well the model predicts:

- Open the `Scatter plot` viewer first
- At the bottom of the viewer, select `Sales` as the variable for the `X` axis
- To the left of the plot, select the predicted variable `outcome` for the `Y` axis
- Lastly, press `R` to plot a regression line

Thus, you have a model that predicts coffee sales based on demographics. We can visualize this on the map now.

## [Predicting Sales on Interactive Map](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=1050s)

Let's go back to the `Google Map` viewer and find `Actions settings` in the viewer's menu. By adjusting them, we can switch the viewer to a special mode and bind our actions to the model. For example, when we double-click, we want to show a prediction. Select this option in the settings and pick your model from the list.

Now, take some time to explore the map. For convenience, select `Full screen` in the viewer settings or press `Alt + F`. Double-click on locations to see if it's worth building a new coffee place there. Find the profitable spots, such as in Times Square, and, perhaps, the underrated ones (by the water).

## [Deploying Models into Apps](https://www.youtube.com/watch?v=tVwpRB8fikQ&t=1141s)
