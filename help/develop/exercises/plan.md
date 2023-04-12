---
title: "Learning plan"
sidebar_position: 1
---

Datagrok provides [programming exercises](exercises.md) to help developers get proficient with the platform. This page
offers a recommended learning plan for completing these exercises in one week.

## Day 1

After Day 1, you will learn how to:

* Import data into the platform
* Create and work with viewers
* Create, edit, and transform a grid
* Create, publish, and share dashboards

The goal for today is to gain hands-on experience with the platform and familiarize yourself with it from a user's
perspective.

1. Register on [public.datagrok.ai](https://public.datagrok.ai) and [dev.datagrok.ai](https://dev.datagrok.ai)
1. Now let's import data into the platform:
   1. [Local files](../../access/import-text.md). Download a
      [dataset](https://raw.githubusercontent.com/datagrok-ai/public/master/packages/ApiTests/files/datasets/demog.csv)
      and import it into the platform using different approaches: drag-and-drop, double-click the "Data" icon on the
      sidebar, select "Data > Open local file", press `Ctrl+O`, and so on.
   1. [File shares](../../access/file-shares). Add the dataset to your *Home* directory.
   1. [Databases](../../access/db-exploration).
      1. Watch a [DB exploration video](https://www.youtube.com/watch?v=YJmSvh3_uCM).
      1. Follow the [Data connectors](https://dev.datagrok.ai/apps/tutorials/DataAccess/DataConnectors) tutorial to
         connect to a Postgres DB and run a query.
   1. [Webservices](../../access/open-api). Browse the [Webservices Manager](https://public.datagrok.ai/webservices)
      (Data > Webservices).
1. Viewers:
   1. Learn how to work with visualizations on the platform by completing these tutorials:
      1. [Scatter plot](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/ScatterPlot)
      1. [Viewers](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Viewers)
      1. [Filters](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Filters)
      1. [Embedded viewers](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/EmbeddedViewers)
   1. Watch an [Interactive data visualization video](https://www.youtube.com/watch?v=67LzPsdNrEc)
   1. Read [documentation](../../visualize/viewers/viewers.md)
   1. Add different viewers using the previously imported dataset
1. Grid:
   1. Read [documentation](../../visualize/viewers/grid)
   1. Watch videos:
      1. [Aggregation](https://www.youtube.com/watch?v=1EI1w2HECrM)
      1. [Using formulas in calculated columns](https://www.youtube.com/watch?v=-yTTaS_WOU4)
      1. [Joining tables](https://www.youtube.com/watch?v=dlbK2Zo-eng)
   1. Create a grid. Edit values, add/remove rows/columns, aggregate, join.
   <!-- TODO: Convert to an exercise -->
1. Dashboards:
   1. Complete the [Dashboards](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Dashboards) tutorial
   1. Create a dashboard and share it with your mentor

## Day 2

After Day 2, you will:

* Learn about functions and scripting
* Write a script and test it

The goal for today is to gain experience with the platform by learning about its more advanced features.

1. Functions
   1. Watch a [First-class functions video](https://www.youtube.com/watch?v=p7_qOU_IzLM&t=724s)
   1. Read documentation:
      1. [Functions](../../datagrok/functions/functions.md)
      1. [Parameter annotation](../../datagrok/functions/func-params-annotation.md)
      1. [Parameter enhancement](../../datagrok/functions/func-params-enhancement.md)
      1. See the list of supported [constants](../../transform/functions/constants.md) and
         [operators](../../transform/functions/operators.md). Explore the groups of functions in the table of contents
         (statistical, math, text, conversion, etc.).
   1. Browse the [Functions](https://public.datagrok.ai/functions) view and filter functions by tags (`#stats`, `#math`,
      `#text`, `#conversion`, etc.). Run simple functions, e.g., `Abs`, `Add` or `And`, from the context panel and the
      console.
1. Scripting
   1. Complete the [Scripting](https://dev.datagrok.ai/apps/tutorials/MachineLearning/Scripting) tutorial
   1. Read documentation:
      1. [Scripting for non-developers](../../compute/scripting-for-non-developers.mdx)
      1. [Scripting](../../compute/scripting.md)
   1. Create a JavaScript script and share it with your mentor:
      1. Cylinder Volume Calculator takes parameters `radius` and `height` and returns `volume` rounded to 2 decimal
         places.
      1. Follow the [instructions](../how-to/add-package-tests.md#testing-functions) to add tests to the script (use
         function `Round10` and constant `PI` in your tests)
      <!-- TODO: Convert to an exercise -->

## Day 3

After Day 3, you will:

* Set up development environment
* Create your first project

The goal for today is to set up the environment and learn how to create simple projects.

1. Check out these links:

   1. Wiki: <https://datagrok.ai/help>
   1. JS API reference: <https://datagrok.ai/js-api/>
   1. API samples: <https://public.datagrok.ai/js>

   Developers use these resources on a regular basis, so consider adding them to your bookmarks.
1. Read our [Contributor’s guide](https://github.com/datagrok-ai/public/blob/master/CONTRIB.md) and [Git
   policy](../advanced/git-policy.mdx)
1. Clone our [public repository](https://github.com/datagrok-ai/public):

   ```sh
   git clone https://github.com/datagrok-ai/public.git
   ```

1. Create a branch for exercises:

   ```sh
   git branch <first letter of your name><your surname>/exercises
   ```

   Example: `jdoe/exercises`
1. Follow instructions from [Setting up development environment](../set-up-environment.md)
1. Create your first project following instructions given in [Exercises](exercises.md#setting-up-the-environment)
1. Add package unit tests:
   1. Read about [Package tests](../how-to/test-packages.md)
   1. Follow the [instructions](../how-to/add-package-tests.md#adding-unit-tests) to add sample tests
   1. Publish your package to the platform and run tests using different approaches described in the above instructions

At the end of the day share the package you created with your mentor and push your branch to GitHub.

## Day 4
