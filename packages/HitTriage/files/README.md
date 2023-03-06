# HitTriage

HitTriage helps chemists assess the quality of hits and decide which compounds should make it to the next stage. It does
it in a managed, reproducible manner, with the triage template consisting of separate steps.

Technically, HitTriage is a [package](https://datagrok.ai/help/develop/develop#packages)
for the [Datagrok](https://datagrok.ai) platform that contains the the `HitTriage` application. There are multiple ways
to start it:

* Open the `Apps` pane and double-click on the `HitTriage` app
* Open it via the direct link: [https://public.datagrok.ai/apps/HitTriage](https://public.datagrok.ai/apps/HitTriage)

There are two outcomes for each project:

1. Compounds that proceed to the next steps
2. Rules used to filter out compounds

## Project Template

Project template defines the structure associated with the following:

* Data Ingestion
  * Data query + pre-processing
* Data Expansion
  * Calculated Properties
  * Data from External Sources
* Property Filters
* Outcome Report Format
* Actions and Triggers (notifications, etc)

The project template is structured and persisted in a way that allows capturing of decision and "replaying" the whole
pipeline on demand, should any step change and new source data appear.

## Data Access

Datagrok's built-in [data query](../../help/access/data-query.md) framework is used as a foundation for data access.
Out-of-the box, it provides the following capabilities that are relevant to this project:

1. Ability to access any data source (db, files, web service) efficiently
    1. Managed connection
2. Parameterization
3. Reproducibility
4. Audit and data lineage

## Data Expansion

During the data expansion step, the system brings in additional information necessary to make a decision on the
compound. Depending on the nature of this information, it would be either calculated, or retrieved from external
sources.

For calculated data, Datagrok's functions framework is used. It allows to define operations (such as computing RDKit
descriptors) that would be executed on the set of molecules.

For retrieving additional information from external sources, the system uses
[parameterized queries](../../help/access/data-query.md).

## Data Storage

The data is stored in the Postgres database. RDKit cartridge is used for efficient searches. Templates are stored in the
same database.

## User Interactions

Users would use the built-in Datagrok facilities for filtering that allow interactive filtering for both compounds and
calculated properties. The UI would look like that:

![hit-triage-filtering](images/hit-triage-filtering.png)

Users will be able to filter out compounds manually, with an optional explanation.

Once the filters are setup, a user clicks the "Submit" button, at which point both hits and decisions get stored in the
database, a report is generated, and a notification is sent to the project team.

## Capturing Results

There are two kinds of results: hits and decisions. Both are stored in the database.

## Advanced Analytics

The following HitTriage-related advanced analytics functionality is available out-of-the-box:

* Structure rendering (RDKit or OpenChemLib)
* Substructure search (either server-side or client-side)
* Interactive data exploration
* Clustering
* Property Calculators
* Chemical Space Exploration
* Dimensionality Reduction (t-SNE, SPE)
* Multivariate Analysis

In addition, the solution takes advantage of the built-in data augmentation system. Whenever a user clicks on a
structure, the relevant information (assays, projects, dose-response curves, etc)
gets shown in the context panel on the right.

## Visualization and Reporting Outcomes

* Built-in [viewers](../../help/visualize/viewers.md)
* Export to Excel / PDF
* Share outcomes as URLs

## Security and Privileges

Datagrok's built-in [privileges system](https://datagrok.ai/help/govern/security)
is used for authentication and authorization within HitTriage. With some app-specific customizations, it allows to

* Define groups and associate them with privileges for objects (HitTriage projects, templates, etc)
  and corresponding actions (create, delete, submit, etc)

## Notifications and Statuses

In the context of HitTriage, the built-in notification system is used for the following:

* Requests for action (review, sign-off, etc)
* Requests for privileges and group memberships (routed to group admins)

Each notification

## External APIs, integration and extension points

While the system fit for the purpose of hit triage, it is also designed in a way that allows to easily customize it in
order to integrate with different systems and processes within enterprises. Here are some of the integration and
extension points:

* Postgres DB: could be directly accessed for integration purposes
* Data Ingestion - ability to define company-specific data ingestion routines
* Data Expansion - ability to define custom calculations (including Python/RDKit, R, etc)
* Datagrok REST API - programmatically access exposed functions
* Datagrok events: ability to customize app behavior (potentially with other packages)
