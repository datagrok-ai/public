# HitTriage

HitTriage helps chemists assess the quality of hits and
decide which compounds should make it to the next stage.
It does it in a managed, reproducible manner, with the triage
template consisting of separate steps. 

There are two outcomes for each project:
1. Compounds that proceed to the next steps
2. Rules used to filter out compounds

## Project Template

* Data Ingestion
* Data Expansion
* Data Storage

## Data Access

Datagrok's built-in [data query](../../help/access/data-query.md) framework is 
used as a foundation for data access. Out-of-the box, it provides the following
capabilities that are relevant to this project:

1. Ability to access any data source (db, files, web service) efficiently
   1. Managed connection 
2. Parameterization
3. Reproducibility
4. Audit and data lineage

## Data Expansion

During the data expansion step, the system brings in additional information 
necessary to make a decision on the compound. Depending on the nature
of this information, it would be either calculated, or retrieved from external
sources.

For calculated data, Datagrok's functions framework is used. It allows to 
define operations (such as computing RDKit descriptors) that would be executed
on the set of molecules.

For retrieving additional information from external sources, the system uses
[parameterized queries](../../help/access/data-query.md).

## Data Storage

The data is stored in the Postgres database. 
RDKit cartridge is used for efficient searches.
Templates are stored in the same database.

## User Interactions

## Capturing Results

There are two kinds of results: molecules annotated with the 

## Visualization and Reporting Outcomes

## Security and Privileges

Datagrok's built-in [privileges system](https://datagrok.ai/help/govern/security)
is used for authentication and authorization within HitTriage. With some 
app-specific customizations, it allows to

* Define groups and associate them with privileges for objects (HitTriage projects, templates, etc)
   and corresponding actions (create, delete, submit, etc)
* 

## Notifications and Statuses

In the context of HitTriage, the built-in notification system is used for the following:

* Requests for action (review, sign-off, etc)
* Requests for privileges and group memberships (routed to group admins)

Each notification

## External APIs, integration and extension points

While the system fit for the purpose of hit triage, it is also designed in a way
that allows to easily customize it in order to integrate with different systems
and processes within enterprises. Here are some of the integration and 
extension points:

* Postgres DB: could be directly accessed for integration purposes
* Data Ingestion - ability to define company-specific data ingestion routines
* Data Expansion - ability to define custom calculations (including Python/RDKit, R, etc)
* Datagrok REST API - programmatically access exposed functions
* Datagrok events: ability to customize app behavior (potentially with other packages)