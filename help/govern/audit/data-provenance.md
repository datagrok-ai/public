---
title: "Data provenance"
---

Data provenance lets you trace every dependency behind a result: the queries that
retrieved the raw data, the transformations applied, the scripts executed, the
predictive models used, and the datasets those models trained on.

Data provenance enables data scientists to reason about results, especially when they do not work in isolation. It also
makes debugging the pipelines a lot easier, since you can go back and see how a particular change impacted the result.

## Self-documented data flows

Since Datagrok already has a lot of information related to the data flows, in many cases there is no need to create and
maintain documentation for the data flows. Datagrok offers several out-of-the-box jobs that build visual representations
of the data flows. They answer the following questions:

* Who are end users of a particular database / connection / query?
* Which data sources does a particular group of users use?
* Which algorithms and visualizations are used on the data retrieved from the particular query?
