---
title: "Learning plan"
---

Datagrok provides [programming exercises](exercises.md) to help developers get proficient with the platform. This page
offers a recommended learning plan for completing these exercises.

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
   1. Local files. Download [dataset] and import it into platform using different approaches (drag-and-drop, File menu)
   1. [File shares](https://datagrok.ai/help/access/file-shares).
   1. [Databases](https://datagrok.ai/help/access/db-exploration). Watch [DB Exploration
      video](https://www.youtube.com/watch?v=YJmSvh3_uCM). Connect to some DB on public, run query.
   1. [Webservices](https://datagrok.ai/help/access/open-api)
1. Viewers:
   1. Tutorials:
      1. [Scatter plot](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/ScatterPlot)
      1. [Viewers](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Viewers)
      1. [Filters](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Filters)
      1. [Embedded viewers](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/EmbeddedViewers)
   1. [Interactive data visualization video](https://www.youtube.com/watch?v=67LzPsdNrEc)
   1. [Scripting viewers video](https://www.youtube.com/watch?v=jHRpOnhBAz4)
   1. Add different viewers using the previously imported dataset
1. Grid:
   1. Read [documentation](https://datagrok.ai/help/visualize/viewers/grid)
   1. Watch videos:
      1. [Aggregation](https://www.youtube.com/watch?v=1EI1w2HECrM)
      1. [Using formulas in calculated columns](https://www.youtube.com/watch?v=-yTTaS_WOU4)
      1. [Joining tables](https://www.youtube.com/watch?v=dlbK2Zo-eng)
   1. Create a grid. Edit values, add/remove rows/columns, aggregate, join.
1. Dashboards:
   1. Complete the [Dashboards](https://dev.datagrok.ai/apps/tutorials/ExploratoryDataAnalysis/Dashboards) tutorial
   1. Create a dashboard and share it with your mentor

## Day 2

After Day 2, you will:

* Set up development environment
* Create your first project

The goal for today is to set up the environment and learn how to create simple projects.

1. Check out these links:
   1. Wiki: <https://datagrok.ai/help>
   1. JS API reference: <https://datagrok.ai/js-api/>
   1. API samples: <https://public.datagrok.ai/js> Developers use these resources on a regular basis, so consider adding
   them to your bookmarks.
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

## Day 3
