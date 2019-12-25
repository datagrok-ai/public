<!-- TITLE: Data Provenance -->
<!-- SUBTITLE: -->

## Data Provenance

Data provenance is the ability to fully understand everything that the result depends upon. This
includes queries that were used to retrieve the initial raw data, transformations applied to the
data, scripts that were executed, predictive models that were used, datasets on which
these models were trained, etc. 

Data provenance enables data scientists to reason about results, especially when they do 
not work in isolation. It also makes debugging the pipelines a lot easier, since you can go back
and see how a particular change impacted the result.

## Self-Documented Data Flows

Since Grok already has a lot of information related to the data flows, in many cases there is no need 
to create and maintain documentation for the data flows. Grok offers several out-of-the-box 
jobs that build visual representations of the data flows. They answer the following questions:
* Who are end users of a particular database / connection / query?
* Which data sources does a particular group of users use?
* Which algorithms and visualizations are used on the data retrieved from the particular query? 

See also:
* [Data Pipelines](../entities/data-pipeline.md)
