<!-- TITLE: Enterprise Evaluation FAQ -->
<!-- SUBTITLE: -->

* [Architecture](architecture.md)
    * [Data flows]()
    * [Deployment](architecture.md#deployment)

* Security
    * [Security overview](#security)
    * [Authentication](../../govern/authentication.md) and [authorization](../../govern/authorization.md)
    * [Encryption at rest](#encryption-at-rest) and [encryption in transit](#encryption-in-transit)

* Enterprise Readiness
    * [Logging and monitoring](#logging-and-monitoring) (using standard AWS tools)
    * [Backup and restore](#backup-and-restore)
    * [Disaster recovery (HA/DR)](#disaster-recovery)
    * [Infrastructure as code](#infrastructure-as-code) (ability to deploy using standard DevOps tools)

* Interoperability
    * Calling web services from client and server with proper auth and integrate that data with Datagrok api
    * Interacting with Datagrok 
        * [server api](#server-api)
        * [client api](../grok-api.md)
    * Connecting to common data sources
        * [relational databases](https://youtu.be/YJmSvh3_uCM)
        * [local files](https://datagrok.ai/img/slides/access-file-formats.mp4) 
        * datastore files
    * [Embedding a Datagrok visualization into a custom web application](https://datagrok.ai/embed_test.html)
    * [Embedding a custom visualization into Datagrok](../../viewers/markup.md) 

* Developer experience
    * [Debug environment when developing customizations](https://youtu.be/PDcXLMsu6UM)  
    * [Devops process including deployment of packages, their dependencies and versioning](../dev.md)  
    * [Concurrent work by team of developers](../dev.md#development)  

* Scalability and Performance
    * Maximum dataset sizes
    * Stability under concurrent user load
    * Scaling and stability under load
[ComputeVM Scaling](../../compute/compute-vm.md#scalability)

* Extensibility
    * [Creating custom visualizations](https://github.com/datagrok-ai/public/tree/master/packages/Sequence)
    * [Creating custom server side components](https://github.com/datagrok-ai/public/tree/master/packages/Pedometer)
    * [Creating custom scripts](https://datagrok.ai/help/compute/scripting) and utilizing them in other components
    * [Ability to reskin Datagrok to appear as a fit-for-purpose web application](https://public.datagrok.ai/apps/spgi) 
    * [Ability to build custom application including data entry, workflow, data model, state management, persistence, etc](https://github.com/datagrok-ai/public/tree/master/packages) 

* Frontend
    * Capacity of holding data in client with differently sized data sets related to project scenarios, sizes detailed in data sources section below
        * High-throughput screening
        * Virtual screening (minimal requirement)
        * Large-chemical spaces
    * Ability to cross-connect multiple large data tables 
    * Visualization of data sets in different plots 
        * High-throughput screening 
        * Virtual screening
    * 2D | 3D structure rendering and interaction: 
        * Ability to customize structure rendering, and interact with it
        * Atom-based annotations
        * 2D-3D connectivity
    * Interactivity
        * [Live data masking](https://youtu.be/67LzPsdNrEc)
        * [Filter by selection](https://youtu.be/67LzPsdNrEc)
        * Different input methods (2D drawing etc.)
        * API calls
    * Ability to plugin non-native Datagrok pieces (e.g. react containers) and interact with Datagrok frontend (eg react containers) 
    * Scalability of frontend scripting functionality

## Encryption at Rest

For AWS deployment, we rely on Amazon's built-in encryption for 
[RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html) 
and 
[S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/dev/bucket-encryption.html).
  
## Server API

Datagrok client side uses HTTP rest API to interact with server side. Authentication token must be passed to access all features.
[Proof of concept video](https://www.youtube.com/watch?v=TjApCwd_3hw)
  
## Encryption in Transit

All client-server communications use the [HTTPS](https://en.wikipedia.org/wiki/HTTPS) protocol, 
which means it is secure and encrypted.

## Logging and Monitoring

## Backup and Restore

Amazon has scheduled backup for RDS and S3, but you can backup and restore RDS database as usual Postgres database.

## Disaster Recovery

Datagrok supports Docker installation, Amazon cluster will immediately restart failed instance 
https://www.youtube.com/watch?v=oFs9RShkHT8

## Infrastructure as Code