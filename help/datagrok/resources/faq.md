---
title: FAQ
sidebar_position: 2
---
##  Drug design and discovery

### Assay data processing

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

## Deploy

### Enterprise readiness

**Q:** How is logging and monitoring handled?  
**A:** Datagrok uses [AWS CloudWatch](https://aws.amazon.com/cloudwatch/) to collect and analyze logs and metrics.  

**Q:** How are backups and restores managed?  
**A:** AWS provides automated RDS and S3 backups; RDS can also be restored as a standard PostgreSQL database.  

**Q:** What is the disaster recovery strategy?  
**A:** High availability is supported via Docker and AWS clusters that automatically restart failed instances.  
See [Disaster recovery demo](https://www.youtube.com/watch?v=oFs9RShkHT8) for more information.

## Access

**Q:** What data sources can Datagrok connect to?  
**A:** Datagrok can connect to a comprehensive range of data sources across multiple categories: [database connections](../../access/databases/connectors/connectors.md), [file storage systems](../../access/files/files.md), [web services and APIs](../../access/open-api.md). Datagrok supports 50+ [file formats](../../access/files/supported-formats.md), including domain-specific like SDF, FASTA, and others.

**Q:** Can I connect Datagrok to BigQuery or other data stores?  
**A:** Yes. Datagrok connects to [BigQuery](../../access/databases/connectors/bigquery.md) and most other [popular databases](../../access/databases/connectors/connectors.md) out of the box. In addition, any machine-readable data source can be easily integrated.

## Govern

**Q:** Does Datagrok follow secure development standards and industry best practices?  
**A:** Yes, Datagrok adheres to secure development and industry best practices across infrastructure, development, and enterprise security:

- Security-first infrastructure – designed with security from the ground up. Features include secure credentials management, flexible authentication (OAuth, SSO, Active Directory), and role-based access control.  
  [Learn more](../../deploy/GCP/deploy-gcp-gke-terraform.md#security) | [Authentication & authorization](../../develop/how-to/apps/build-an-app.md#authentication) | [Role-based access](../../compute/compute.md#privileges-and-visibility)

- Quality assurance – Automated check for vulnerabilities (Snyk, Grype) and multi-layered testing including unit, integration, UI, and performance tests.  
  [Learn more](../../develop/qa/quality-assurance.md#continuous-integration-and-deployment-system) | [Testing](../../develop/qa/quality-assurance.md#automated-testing)

- Secure development lifecycle – CI/CD pipelines enforce automated build, test, and security checks; semantic versioning; and secure credential management.  
  [Learn more](../../develop/develop.md#continuous-integration) | [Version control](../..//develop/dev-process/versioning-policy) | [Credential management](../../develop/how-to/packages/manage-credentials.md)

### Authentication and access control

**Q:** How does Datagrok handle authentication and authorization?  
**A:** Datagrok uses role-based access control and integrates with enterprise identity providers such as LDAP, SSO, and OAuth.  
Learn more about [Access control](../../govern/access-control/access-control.md)

**Q:** Is data encrypted at rest?  
**A:** Yes, Datagrok relies on Amazon's built-in encryption for [RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html) and [S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucket-encryption.html).

**Q:** Is data encrypted in transit?  
**A:** Yes, all client-server communications use [HTTPS](https://en.wikipedia.org/wiki/HTTPS), which means it is secure and encrypted. 

**Q:** How is MFA enforced for GitHub accounts?  
**A:** GitHub requires all code contributors to enable two-factor authentication (2FA) as of March 2023. Developers comply by using device-tied passcodes.  

**Q:** Who has access to Bitbucket, and what are the authentication requirements?  
**A:** Only core developers have access. MFA is being enabled to strengthen authentication. Bitbucket also provides built-in security scanning features.  

### Vulnerability and patch management

**Q:** How does Datagrok manage vulnerabilities in the application and cloud infrastructure?  
**A:** We monitor CISA alerts and run daily [Snyk](https://snyk.io/) scans on container builds. Vulnerabilities are triaged, remediated via infrastructure-as-code pipelines, and verified through CI/CD testing. Customers are notified and supported in version upgrades as needed.  
Learn more about [Infrastructure](../../develop/under-the-hood/infrastructure.md)

**Q: How are OS and server patching handled?**  
**A:** Jenkins and development servers are rebuilt quarterly with a fresh OS. Documentation for builds is maintained in internal repositories.  
Learn more about [Deployment](../../deploy/deploy.md)

### Endpoint security

**Q:** What controls mitigate risks on BYOD devices?  
**A:** We have a policy that all devices must have active malware protection, updated signatures, and drive encryption enabled.  

**Q:** How is data exfiltration via removable media prevented?  
**A:** We have a policy that sensitive data must not be stored on removable media. Passwords or credentials must not be transmitted unencrypted.  

### Endpoint Detection and Response (EDR)

**Q:** What is Datagrok’s current EDR approach?  
**A:** No dedicated EDR is deployed on associates’ devices. However, all devices must have up-to-date and active malware protection. Associates are required to report any suspicious activity via the alert channel. Incident response team triages and executes the incident management process 

**Q:** How are security logs collected and monitored?  
**A:** AWS resources use centralized logging. Currently, there is no proactive log review, but reporting of failed logins with alert thresholds is being implemented.  

### Is my data private?

**Q:** What happens to my data when I open a local file in Datagrok?  
**A:** When you open a local file in Datagrok (like dragging and dropping a file to your browser), you can analyze it without saving. This data stays in your browser's memory and isn't sent to the server unless you run resource-intensive server-side computations. Your data is gone when you close the browser tab. To save your work, you need to upload it to the server. Note that uploading data does not make it accessible to others. Your data stays private and visible to you only until you explicitly share it. Learn how to [save](../concepts/project/project.md#saving-entities-to-projects) and [share](../navigation/basic-tasks/basic-tasks.md#share) data.

**Q:** What data or telemetry is sent back to Datagrok? What egress ports/protocols are used?  
**A:** Datagrok does not send anything without user permissions. If user reports error and explicitly checks "Send report back to Datagrok", email is sent to feedback@datagrok.ai. Datagrok can download images from Docker Hub or packages from NPM.

## Visualize

**Q:** Can I create custom visualizations?  
**A:** Yes, Datagrok offers multiple options for creating and embedding custom visualizations:

- Custom viewers via JavaScript API – [build custom viewers](../../develop/how-to/viewers/develop-custom-viewer.md) using the
[JavaScript API](../../develop/packages/js-api.md)
- Scripting viewers – use R, Python, or Julia to embed visualizations via
[scripting](../../compute/scripting/scripting.mdx).
- Third-party libraries – integrate frameworks like ECharts, D3.js, Circos, or three.js 
(see [Charts package](https://github.com/datagrok-ai/public/tree/master/packages/Charts))
- Custom file and folder viewers – extend Datagrok with  
file viewers for specific formats 
(e.g., PDB files via [NGL viewer](../../visualize/viewers/ngl.md) through the [Biostructure Viewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#biostructure-viewer) package)
- [Custom cell renderers](../../develop/how-to/grid/custom-cell-renderers.md) - create custom visualization for cells in data [grid/table](../../visualize/viewers/grid.md)

**Q:** How can I visualizing old and new data side by side?  
**A:** Multiple [viewer](../../visualize/viewers/viewers.md) and grids can be combined in a single [dashboard](../../datagrok/concepts/project/dashboard.md) to display different datasets side by side. 
Each viewer (scatterplots, charts, etc.) has configurable settings, allowing plots to source data from different tables. There are also multiple ways to [join](../../transform/join-tables.md), [link](../../transform/link-tables.md), [pivot](../../transform/aggregate-rows.md), and apply other [transformations](../../transform/transform.md) 
to data within the intuitive UI. 

## Collaborate

**Q:** How can I link raw data, metadata, and analysis results to compounds or sequences so that this information can be easily recalled and shared across Datagrok dashboards?  
**A:** The [Sticky meta](../../govern/catalog/sticky-meta.md) feature allows tagging compounds, sequences, or other entities 
with meta-information, such as calculated properties, analysis results, or comments. Tagged information can be recalled anywhere in Datagrok 
when the same entity appears, supporting collaborative workflows and effective communication between team members.  

## Develop

### Architecture

**Q:** How is Datagrok architected?  
**A:** Datagrok is built as a modular, service-oriented platform with client-server separation.  
Learn more about [Datagrok's architecture](../../develop/under-the-hood/architecture.md#goals)

**Q:** How is Datagrok deployed?  
**A:** Datagrok supports Docker, Kubernetes, on-premises, and cloud deployments (AWS, GCP, Azure).  
Learn more about [Deployment](../../deploy/deploy.md)

### Requirements

**Q:** What are the browser requirements?  
**A:** Datagrok is compatible with modern browsers.  Popular choices are Chrome, Edge, or Safari.

### Interoperability

**Q:** Can Datagrok call external web services?  
**A:** Yes, both client and server can call authenticated web services and integrate responses with the Datagrok API.  
See [Developer guide](../../develop/develop.md) for details

**Q:** How do I interact with Datagrok from the server?  
**A:** Use Datagrok’s REST server API with authentication tokens.  
See [Proof of concept video](https://www.youtube.com/watch?v=TjApCwd_3hw) for more information

**Q:** How do I interact with Datagrok from the client?  
**A:** Use the JavaScript API for custom applications and extensions.  
Learn more about  [JavaScript API](../../develop/packages/js-api.md)

### Developer experience

**Q:** How do I develop and debug Datagrok customizations?   
**A:** Follow the [step-by-step guide](../../develop/how-to/packages/create-package.md) for creating your own package. If your package is based on a non-standard template, you may need to configure it differently for [debugging](../../develop/advanced/debugging.md).  
See also [Debug demo](https://youtu.be/PDcXLMsu6UM) for details

**Q:** How does the DevOps process work?  
**A:** The process includes package deployment, dependency management, and versioning.  
See [Developer guide](../../develop/develop.md) for details

**Q:** Can multiple developers work concurrently?  
**A:** Yes, with versioning, branching, and package management.  
Learn more about [Concurrent development](../../develop/develop.md#development)

### Scalability and performance

**Q:** What are the maximum dataset sizes?  
**A:** Datagrok is designed to handle millions of data points interactively for visualization and exploration. See [Why Datagrok?](../../datagrok/datagrok.md#why-datagrok) for details 

**Q:** How stable is the platform under heavy load?  
**A:** Datagrok is designed for horizontal scaling, ensuring stability with concurrent users.  
See [Scaling and stability](../../develop/under-the-hood/infrastructure.md#scalability) for more information

### Extensibility

**Q:** Can I create custom visualizations?  
**A:** Yes, Datagrok offers multiple options for creating and embedding custom visualizations.  
See [Custom visualizations](#visualize) for details

**Q:** Can I build server-side components?  
**A:** Yes, custom back-end logic can be added.  
See [Admetica example](https://github.com/datagrok-ai/public/tree/master/packages/Admetica)

**Q:** Can I add scripts and reuse them in components?  
**A:** Yes, Datagrok supports multiple scripting languages.  
Learn more about [Scripting](../../compute/scripting/scripting.mdx)

**Q:** Can I reskin Datagrok?  
**A:** Yes, you can build tailored UI applications.  
See [Example app](https://public.datagrok.ai/apps/HitTriage/HitTriage?browse=apps)

**Q:** Can I build full custom applications?  
**A:** Yes, including workflows, data models, state management, and persistence.  
See [Custom application packages](https://github.com/datagrok-ai/public/tree/master/packages)

### Frontend

**Q:** What frontend capabilities does Datagrok provide for scientific data analysis?     
**A:** Datagrok delivers a high-performance, interactive frontend for handling, visualizing, and analyzing large-scale scientific datasets. Key capabilities include:

- Data capacity and client-side handling – load and explore millions of rows directly in the browser using Datagrok’s in-memory data engine. Supports high-throughput and virtual screening scenarios with GPU-accelerated visualizations.  
  See also: [Example: 2.7M ChEMBL molecules](../../datagrok/datagrok.md#why-datagrok) | [WebGPU acceleration](../../visualize/viewers/scatter-plot#webgpu-acceleration)
- Cross-connecting multiple large data tables – perform sophisticated data manipulation including joins, filters, and aggregation. Linked tables enable interactive cross-referencing across datasets.  
  Learn more about [Transformations](../../transform/transform.md)
- Visualization capabilities – Over 50 interactive [viewers](../../visualize/viewers/viewers.md) for synchronized dashboards, including [chemically-aware viewers](../../datagrok/solutions/domains/chem/chemically-aware-viewers.md)  
- 2D/3D structure rendering and interaction.  
See [Cheminformatics](../../datagrok/solutions/domains/chem/chem.md) for details
- Interactivity features – live data masking, filter by selection, synchronized updates across viewers, multiple input methods including 2D sketching.  
See [Cheminformatics](../../datagrok/solutions/domains/chem/chem.md) for details
- API integration and extensibility. See [JavaScript API](../../develop/packages/js-api.md) for details
- Frontend scripting scalability. See [Compute](../../compute/compute.md) for details

<!--Q: Can infrastructure be managed as code?  
**A:** Yes, using Docker-compose, Jenkins, CloudFormation, or Terraform.  
Datagrok Docker containers are built using Jenkins, all software are upgraded and patched on every build.
Docker-compose manifest is used to describe and deploy Datagrok applications.
Also, there are multiple advanced options to deploy application:
CloudFormation template to deploy to AWS ECS
Terraform scripts to deploy to AWS ECS-->

<!--## Deployment

**Q:** Do you have a guide for installing and deploying Datagrok on a GCP Kubernetes cluster?   
**A:** The fastest way to deploy Datagrok for evaluation is by using Docker Compose on a virtual machine. This setup takes just a minute or two.    
[Learn more](https://datagrok.ai/help/develop/admin/docker-compose)-->



