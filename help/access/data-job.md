<!-- TITLE: Data job -->
<!-- SUBTITLE: -->

# Data job

Data job defines all actions that are needed to produce a dashboard.

Each data job consists of the following parts that are executed to produce a dataset:

* Data queries (each with the corresponding connection)
    * Cross-database querying
* Transformations applied to the data
    * Standard, can be defined via UI (such as "Remove rows where $AGE is null")
    * Out-of-the box domain-specific analyses
        * Statistical process control
        * Credit card transactions: fraud detection
        * NLP: Sentiment analysis
        * DSP: Spectral analysis
        * ... and hundreds of other algorithms
    * Dart scripts
    * R / Python / Julia scripts
    * External programs
* Visual layouts applied to the tables
* Definition of the summary page

To browse all data jobs, open **Admin | Data Jobs**. To browse data job executions, open
**Admin | Data Job Runs**.

## Access control

Data jobs are first-class entities in the Datagrok platform, and as such are subjects to the standard checks and
routines performed against them whenever they are used in the specific context. Some of the commonly used privileges
are:

* `create`
* `edit`
* `delete`
* `query`

Those privileges can be given to individuals or to groups (which can be defined via dynamic filters)
. For more information on the access privilege model, refer to the Datagrok - Access Privileges page.

## Filtering

You can use these fields to filter jobs using [smart search](../overview/smart-search.md):

| Field       | Description                                                            |
|-------------|------------------------------------------------------------------------|
| ID          |                                                                        |
| name        |                                                                        |
| runs        | list of [DataActionRun](../overview/functions/function-call.md) object |
| queries     | list of [DataJob](data-job.md) object                                  |
| createdOn   |                                                                        |
| updatedOn   |                                                                        |
| author      | [User](../govern/user.md) object                                       |
| starredBy   | [User](../govern/user.md) object                                       |
| commentedBy | [User](../govern/user.md) object                                       |
| usedBy      | [User](../govern/user.md) object                                       |

See also:

* [Data pipeline](data-pipeline.md)
* [Data connection](data-connection.md)
* [Data query](data-query.md)
* [Function call](../overview/functions/function-call.md)
