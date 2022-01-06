<!-- TITLE: Architecture Details -->
<!-- SUBTITLE: -->

## Architecture

![](network-diagram1.png)

Datagrok installation consists of two virtual machines:

* [Datagrok Virtual Machine](#datagrok-virtual-machine)
* [Compute Virtual Machine](#compute-virtual-machine)

Also, it needs [database](#database) and [persistent file storage](#storage). Both of the virtual machines includes
several docker containers and can be deployed in
[AWS EC2](deploy-amazon-ec2.md), [AWS ESC](deploy-amazon-ecs.md)
or [regular host machine](deploy-regular.md).

All docker containers are based on Ubuntu 20.04 and use the latest software available.

## Datagrok virtual machine

This machine is the heart of the platform and is required for all activities.

Datagrok Virtual machine components are:

* [Datlas](#datlas)
* [ElasticSearch](#elasticsearch)
* [Credentials Management Service](../../govern/security.md#credentials). Can be installed as a separate service in
  separate container with a separate database.
* [Grok Connect](../../access/data-connection.md). Separate container with Java-based data connectors to 20+ databases.
* Web application
* Nginx server

### Datlas

Grok server (also referred to as "Datlas") is a Dart stand-alone application that exposes REST API.

Ports:

* `8082` for default mode
* `8443` for SSL mode
* `8083` for CLI

Datagrok supports several deployment schemas, which can be configured using `GROK_MODE`
and `GROK_PARAMETERS`
environment variables.

`GROK_MODE` possible values:

* start - starts the application without database and storage deployment
* deploy - Datlas will perform the full deploy
* auto - Datlas will check the existing database and storage and perform deployment only if needed

`GROK_PARAMETERS` is a JSON-formatted string with options:

| Option              | Required      | Default | Description                                                       |
|---------------------|---------------|---------|-------------------------------------------------------------------|
| dbServer            | **
Required**  |         | Postgres database server.                                         |
| dbPort              | Optional      | 5432    | Postgres database server port                                     |
| db                  | **
Required**  |         | Datagrok database name                                            |
| dbLogin             | **
Required**  |         | Username to connect to database                                   |
| dbPassword          | **
Required**  |         | Password to connect to database                                   |
| dbSsl               | Optional      | false   | If set to true, TLS connection will be used to connect to database|
| dbAdminLogin        | Optional      |         | Postgres admin username to create user and database for Datagrok  |
| dbAdminPassword     | Optional      |         | Postgres admin password to create user and database for Datagrok  |
| googleStorageCert   | Optional      |         | Access certificate to Google Cloud Storage. If set, GCS will be used for persistent data storage| 
| amazonStorageRegion | Optional      |         | S3 region                                                         |
| amazonStorageBucket | Optional      |         | S3 bucket name                                                    |
| amazonStorageId     | Optional      |         | S3 credential ID, Datagrok will resolve EC2 role if empty         |
| amazonStorageKey    | Optional      |         | S3 credential secret key, Datagrok will resolve EC2 role if empty |
| adminPassword       | Optional      |         | Datagrok admin user password which will be created on first start |
| debug               | Optional      | false   | Extended logging and saving stack traces                          |
| useSSL              | Optional      | false   | If set to true, Datlas serves TLS connections                     |
| certPath            | Optional      |         | Path to the TLS certificate                                       |
| certKeyPath         | Optional      |         | Path to the TLS certificate key                                   |
| certKeyPwd          | Optional      |         | Password to the TLS certificate key                               |

### ElasticSearch

[ElasticSearch](https://www.elastic.co/) provides Datagrok's full-text search. It searches in Wiki, forums, and
datasets.

ElasticSearch is configured with default settings. The Datlas application is the only one that uses it, so ElasticSearch
does not need to be exposed outside the Datagrok docker container.

## Compute virtual machine

Compute Virtual Machine is used for performing on-server computations.

Compute Virtual machine components are:

* [Jupyter Kernel Gateway](#jupyter-kernel-gateway)
* [Jupyter Notebook](#jupyter-notebook)
* [H2O](#h2o)
* [GrokCompute](grok-compute.md)

All components work in separate docker containers.

Also, it needs [Grok Helper](#grok-helper) in every container and
[load balancer](#load-balancer) in front of the components to manage the routing and balancing.

See also: [Compute VM](compute-vm.md)

### Load Balancer

A load balancer provides a single entry point to the Compute Virtual Machine services. Also, Load Balancer allows us to
scale containers by the needs of the project. To speed up code execution Compute Virtual Machine uses local cache files.
To avoid loss of the context, Load Balancer uses sticky sessions, which route all of the user's requests to a specific
server for the duration of the session.

Compute Virtual Machine supports any Load Balancer, including AWS Application Load Balancer. To run ComputeVM locally
[cvm_nginx](https://hub.docker.com/r/datagrok/cvm_nginx) docker image can be used.

Load Balancer routes traffic to containers by location

* `/grok_connect/` -> Grok Connect
* `/notebook/` -> Jupyter Notebook
* `/notebook/helper/` -> Jupyter Notebook Grok Helper
* `/jupyter/` -> Jupyter Kernel Gateway
* `/jupyter/helper/` -> Jupyter Kernel Gateway Grok Helper
* `/` - returns 204 No Content

### Grok Helper

Grok Helper is a small Python-based HTTP daemon that provides REST API endpoints for the Datagrok application. It is
located in every container except Grok Compute, which already has the functionality.

Grok Helper exposes API for the following features:

* Jupyter Notebooks converter (HTML, PDF, etc.)
* Utilities
  - Cache entities managing
  - Python environments managing

### Jupyter Kernel Gateway

Jupyter Kernel Gateway provides the scripting feature for the platform. It includes standard libraries for development.

Available languages are: Python, R, JS, Octave, Julia.

Ports:

* `8889` - Jupyter Kernel Gateway
* `5005` - Grok Helper

### Jupyter Notebook

[Jupyter Notebook](../../compute/jupyter-notebook.md) works in the docker container along with other applications. It
includes standard libraries for development. Grok Helper in Jupyter Notebook container converts notebooks to HTML or PDF
by request.

Available languages are: Python. Optional: R, JS, Octave, Julia.

Ports:

* `8888` - Jupyter Notebook
* `5005` - Grok Helper

### H2O

H2O is a Java fully open-source application, which supports the most widely used statistical and machine learning
algorithms.

H2O is accessed by Datlas directly, bypassing the load balancer. The ending user can use H2O through the Datagrok
platform.

Ports:

* `54321` - H2O
* `5005` - Grok Helper

## Storage

Most of the metadata associated with users, datasets, algorithms, predictive models, etc. is kept in
the [relational database](#database).

The actual tabular data (which can be uploaded or created via other means) is stored externally, in the highly
optimized [d42 format](architecture-design.md#in-memory-database), on file storage chosen by the company. Datagrok
supports the following storages:

* Local File System
* Network shares
* S3
* Google Cloud

By default, Datagrok works with Amazon S3 storage. Also, Datagrok can use a local file system if you run a docker
container without AWS or another cloud infrastructure. See [Run docker image on a regular machine](deploy-regular.md)

## Database

Metadata associated with users, datasets, algorithms, predictive models, etc., are kept in a Postgres database. Having
the data stored in a relational database brings well-defined relations and data consistency, enables us to work
efficiently with complex queries, and adds transactional support.

Datagrok uses PostgreSQL as a database; however, if necessary, we can switch to a scalable solution like Aurora or
CocroachDB without changing much code.

Datagrok Virtual Machine can use single any PostgreSQL instance, including Amazon RDS out-of-the-box. For security,
Datagrok supports a TLS connection to connect to the database.

The schema has the following tables:

| Table                | Comments                                | 
|------------------    |-----------------------------------------|
| chats                | Chats                                   | 
| chats_reads          | Topics read by users                    | 
| chats_watch          | Topics watched by users                 | 
| comment_votes        | Upvoted comments                        | 
| comments             | Chat comments                           | 
| connections          | Data connections                        | 
| db_queries           | History of the database queries         | 
| db_ups               | History of database upgrades            | 
| email_history        | Sent emails                             | 
| entities             | System entities (base table)            | 
| entities_chats       | Entity-related chats                    | 
| entities_types       | Descriptions of entity classes          | 
| entities_types_permissions | Sets of permissions specific to entity classes | 
| event_parameter_values     | Parameter values              | 
| event_parameters     | Event parameters                        | 
| event_types          | Event types                             | 
| events               | Events                                  | 
| favorites            | Entities marked as favorites            | 
| files                | Files indexed by the platform           | 
| groups               | User groups                             | 
| groups_relations     | Nested groups                           | 
| groups_requests      | Requests to join groups                 | 
| jobs                 | Jobs                                    | 
| keys                 | Key for encrypted data connections passwords| 
| meta_params          | Dynamic entity parameters               | 
| models               | Predictive models                       | 
| notebooks            | Jupyter notebooks                       | 
| notebooks_tables     | Notebooks applied to tables             | 
| notification_preferences | User notification preferences   | 
| notification_types   | Notification types              | 
| notifications        | Notifications                           | 
| permissions          | Permissions, that set to objects        | 
| pm_input_columns     | Predictive model inputs                 | 
| pm_output_columns    | Predictive model outputs                | 
| project_layouts      | View layouts in projects                | 
| project_relations    | Projects and entities relations         | 
| projects             | Projects                                | 
| queries              | Database queries                        | 
| scripts              | Scripts                                 | 
| table_columns        | Columns in a table                      | 
| table_queries        | Structured table queries                | 
| tables               | Table infos (note that the data resides externally) | 
| tags                 | Entity tags                             | 
| users                | Users                                   | 
| users_sessions       | User sessions                           | 
| view_layouts         | View layouts                            | 
| view_layouts_columns | Columns referenced by layouts (used for layout suggestions) | 

See also:

* [Versioning policy](../versioning-policy.md)
