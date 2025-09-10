---
title: FAQ
sidebar_position: 2
---
## Enterprise evaluation FAQ

### Architecture

**Q:** How is Datagrok architected?  
**A:** Datagrok is built as a modular, service-oriented platform with client-server separation.  
Learn more about [Datagrok's architecture](../../develop/under-the-hood/architecture.md#goals).

**Q:** How is Datagrok deployed?  
**A:** Datagrok supports Docker, Kubernetes, on-premises, and cloud deployments (AWS, GCP, Azure).  
Learn more about [Deployment](../../deploy/deploy.md).

### Requirements

**Q:** What are the browser requirements?
**A:** Datagrok is compatible with modern browsers.  Popular choices are Chrome, Edge, or Safari.

### Security

**Q:** How does Datagrok handle authentication and authorization?  
**A:** Datagrok uses role-based access control and integrates with enterprise identity providers such as LDAP, SSO, and OAuth.  
Learn more about [Access control](../../govern/access-control/access-control.md).

**Q:** Is data encrypted at rest?  
**A:** Yes, Datagrok relies on Amazon's built-in encryption for [RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html) and [S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucket-encryption.html).

**Q:** Is data encrypted in transit?  
**A:** Yes, all client-server communications use [HTTPS](https://en.wikipedia.org/wiki/HTTPS), which means it is secure and encrypted.  

### Enterprise readiness

**Q:** How is logging and monitoring handled?  
**A:** Datagrok uses [AWS CloudWatch](https://aws.amazon.com/cloudwatch/) to collect and analyze logs and metrics.  

**Q:** How are backups and restores managed?  
**A:** AWS provides automated RDS and S3 backups; RDS can also be restored as a standard PostgreSQL database.  

**Q:** What is the disaster recovery strategy?  
**A:** High availability is supported via Docker and AWS clusters that automatically restart failed instances.  
See [Disaster recovery demo](https://www.youtube.com/watch?v=oFs9RShkHT8) for more information.

<!--Q: Can infrastructure be managed as code?  
**A:** Yes, using Docker-compose, Jenkins, CloudFormation, or Terraform.  
Datagrok Docker containers are built using Jenkins, all software are upgraded and patched on every build.
Docker-compose manifest is used to describe and deploy Datagrok applications.
Also, there are multiple advanced options to deploy application:
CloudFormation template to deploy to AWS ECS
Terraform scripts to deploy to AWS ECS-->

### Interoperability

**Q:** Can Datagrok call external web services?  
**A:** Yes, both client and server can call authenticated web services and integrate responses with the Datagrok API.  
See [Developer guide](../../develop/develop.md) for details.

**Q:** How do I interact with Datagrok from the server?  
**A:** Use Datagrok’s REST server API with authentication tokens.  
See [Proof of concept video](https://www.youtube.com/watch?v=TjApCwd_3hw) for more information.

**Q:** How do I interact with Datagrok from the client?  
**A:** Use the JavaScript API for custom applications and extensions.  
Learn more about  [JavaScript API](../../develop/packages/js-api.md).

**Q:** What data sources can Datagrok connect to?  
**A:** Datagrok can connect to a comprehensive range of data sources across multiple categories: [database connections](../../access/databases/connectors/connectors.md), [file storage systems](../../access/files/files.md), [web services and APIs](../../access/open-api.md). Datagrok supports 50+ [file formats](../../access/files/supported-formats.md), including domain-specific like SDF, FASTA, and others.

**Q:** Can I embed custom visualizations into Datagrok?  
**A:** Yes, Datagrok supports custom viewer integrations.  
Learn more about [Custom visualizations](../../visualize/viewers/markup.md).

### Developer experience

**Q:** How do I develop and debug Datagrok customizations?   
**A:** Follow the [step-by-step guide](../../develop/how-to/packages/create-package.md) for creating your own package. If your package is based on a non-standard template, you may need to configure it differently for [debugging](../../develop/advanced/debugging.md). See also [Debug demo](https://youtu.be/PDcXLMsu6UM) for details.

**Q:** How does the DevOps process work?  
**A:** The process includes package deployment, dependency management, and versioning.  
See [Developer guide](../../develop/develop.md) for details.

Q: Can multiple developers work concurrently?  
**A:** Yes, with versioning, branching, and package management.  
Learn more about [Concurrent development](../../develop/develop.md#development).

### Scalability and performance

**Q:** What are the maximum dataset sizes?  
**A:** Datagrok is designed to handle millions of data points interactively for visualization and exploration. See [Why Datagrok?](../../datagrok/datagrok.md#why-datagrok) for details.  

**Q:** How stable is the platform under heavy load?  
**A:** Datagrok is designed for horizontal scaling, ensuring stability with concurrent users.  
See [Scaling and stability](../../develop/under-the-hood/infrastructure.md#scalability) for more information.

### Extensibility

Q: Can I create custom visualizations?  
**A:** Yes, using Datagrok packages.  
See [BiostructureViewer example](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer).

Q: Can I build server-side components?  
**A:** Yes, custom back-end logic can be added.  
See [Admetica example](https://github.com/datagrok-ai/public/tree/master/packages/Admetica).

Q: Can I add scripts and reuse them in components?  
**A:** Yes, Datagrok supports multiple scripting languages.  
Learn more about [Scripting](https://datagrok.ai/help/compute/scripting).

Q: Can I reskin Datagrok?  
**A:** Yes, you can build tailored UI applications.  
See [Example app](https://public.datagrok.ai/apps/HitTriage/HitTriage?browse=apps).

Q: Can I build full custom applications?  
**A:** Yes, including workflows, data models, state management, and persistence.  
See [Custom application packages](https://github.com/datagrok-ai/public/tree/master/packages).

### Frontend

**Q:** What capabilities does Datagrok provide in the frontend for handling data, visualizations, interactivity, and integration?     
**A:** Datagrok's frontend supports:
* Capacity of holding data in client with differently sized data sets related to project scenarios, sizes detailed in
    data sources section below
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
  * Ability to plugin non-native Datagrok pieces (e.g. react containers) and interact with Datagrok frontend (eg react
    containers)
  * Scalability of frontend scripting functionality

## Is my data private?

**Q:** What happens to my data when I open a local file in Datagrok?

**A:** When you open a local file in Datagrok (like dragging and dropping a file to your browser), you can analyze it without saving. This data stays in your browser's memory and isn't sent to the server unless you run resource-intensive server-side computations. Your data is gone when you close the browser tab. To save your work, you need to upload it to the server. Note that uploading data does not make it accessible to others. Your data stays private and visible to you only until you explicitly share it. Learn how to [save](../concepts/project/project.md#saving-entities-to-projects) and [share](../navigation/basic-tasks/basic-tasks.md#share) data.

**Q:** What data or telemetry is sent back to Datagrok? What egress ports/protocols are used?

**A:** Datagrok does not send anything without user permissions. If user reports error and explicitly checks "Send report back to Datagrok", email is sent to feedback@datagrok.ai. Datagrok can download images from Docker Hub or packages from NPM.

## Data Connections

**Q:** Can I connect Datagrok to BigQuery or other data stores?

**A:** Yes. Datagrok connects to [BigQuery](../../access/databases/connectors/bigquery.md) and most other [popular databases](../../access/databases/connectors/connectors.md) out of the box. In addition, any machine-readable data source can be easily integrated.

<!--## Deployment

**Q:** Do you have a guide for installing and deploying Datagrok on a GCP Kubernetes cluster?

**A:** The fastest way to deploy Datagrok for evaluation is by using Docker Compose on a virtual machine. This setup takes just a minute or two.  

[Learn more](https://datagrok.ai/help/develop/admin/docker-compose)-->

## Assay Data Processing

**Q:** How can I ingest raw assay data into Datagrok?  

**A:** Raw data can be ingested through multiple channels, including drag-and-drop, [file shares](../../access/files/files.md#file-sharing-and-access-control), 
[database connections](../../access/databases/databases.md#connecting-to-database), [OpenAPI](../../access/open-api.md), and integrations with services like 
Benchling or Revvity Signals. Datagrok is developing a Plates application
to support predefined plate templates, batch ingestion, and
integrated analysis in one place.

**Q:** How can I perform curve fitting, normalization, QC calculations from raw assay data in Datagrok?

**A:** The [Curves](https://github.com/datagrok-ai/public/blob/master/packages/Curves/README.md)
plugin provides a complete workflow for converting raw assay data into fitted curves with QC
calculations. It transforms well-level assay data into fitted dose–response curves, supporting
functions such as 4PL and sigmoid, and can be easily extended to additional analyses, including
Km/Vmax, dose ratio, melt curves (DSF Tm), and qPCR, by defining the corresponding functions
and parameters.

For modeling dynamic systems, see the [Diff-Studio](../../compute/diff-studio.md) plugin that
allows users to solve sets of differential equations interactively through an intuitive UI and 
a declarative approach.  

**Q:** How can I visualizing old and new data side by side?  

**A:** Multiple [viewer](../../visualize/viewers/viewers.md) and grids can be combined in a single [dashboard](../../datagrok/concepts/project/dashboard.md) to display different datasets side by side. 
Each viewer (scatterplots, charts, etc.) has configurable settings, allowing plots to source data from different tables. There are also multiple ways to [join](../../transform/join-tables.md), [link](../../transform/link-tables.md), [pivot](../../transform/aggregate-rows.md), and apply other [transformations](../../transform/transform.md) 
to data within the intuitive UI. 

**Q:** How can I link raw data, metadata, and analysis results to compounds or sequences so that this information can be easily recalled and shared across Datagrok dashboards?  

**A:** The [Sticky meta](../../govern/catalog/sticky-meta.md) feature allows tagging compounds, sequences, or other entities 
with meta-information, such as calculated properties, analysis results, or comments. Tagged information can be recalled anywhere in Datagrok 
when the same entity appears, supporting collaborative workflows and effective communication between team members.  