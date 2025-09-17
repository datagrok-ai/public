---
title: FAQ
sidebar_position: 2
format: mdx
---

## Quick reference

* [Supported databases](../../access/databases/connectors/connectors.md)
* [Supported file shares](../../access/files/shares/shares.md)
* [Supported file and data formats](../../access/files/supported-formats.md)
* [Supported visualizations](../../visualize/viewers/viewers.md)

##  Drug design and discovery

##### <b>Q: Can Datagrok handle multiple complex data types in one interface?</b>

Datagrok allows access to many different and complex data types in a single user
interface, supporting use cases like
[DMTA](https://datagrok.ai/solutions/dmta-cycle) and other integrated chemical
and assay data workflows. See [Example: Complex data types in one user interface](https://public.datagrok.ai/p/skalkin.medchem_1/Med_Chem)  

##### <b>Q: Is there a sample dashboard for chemists?</b>

* [This dashboard](https://public.datagrok.ai/p/oahadzhaniandatagrokai.spgi100_2/spgi-100) has multiple visualizations, including chemical spreadsheet, multiple scatterplot views, R-groups, and molecule cards
* [This DMTA spreadsheet](https://public.datagrok.ai/p/skalkin.medchem_1/Med_Chem) shows multiple modalities in one view
* The [Demo App](https://public.datagrok.ai/apps/Tutorials/Demo/Cheminformatics) showcases common analyses like MMP or chemical space, and how they look in Datagrok 

##### <b>Q: Is there a sample dashboard for biologists?</b>

The [Demo App](https://public.datagrok.ai/apps/Tutorials/Demo/Bioinformatics) showcases common analyses like sequence space or Peptide SAR, and how they look in Datagrok.

### Assays, plates, curves

##### <b>Q: Does Datagrok support compound-centric and assay-centric data views?</b> 

Yes. Datagrok supports both:  

* Compound view (compound-centric): multiple endpoints per compound  
[Chemically aware viewers](../../datagrok/solutions/domains/chem/chemically-aware-viewers.md) | 
[Example: rich chemical dataset](https://public.datagrok.ai/p/skalkin.medchem_1/Med_Chem)
* Table view (assay-centric): multiple compounds in a table with configurable endpoints  
[Table view](../../datagrok/navigation/views/table-view.md) | [Forms](../../visualize/viewers/forms.md) 

Users can configure coloring, highlighting, and switch between compounds as rows or columns, and filter, search, and sort the data. [Common actions](../../datagrok/navigation/views/table-view.md#common-actions) | [Substructure search](../../datagrok/solutions/domains/chem/chem.md#substructure-search--filtering)

##### <b>Q: How can I ingest raw assay data into Datagrok?</b>

Raw data can be ingested through multiple channels, including drag-and-drop, [file shares](../../access/files/files.md#file-sharing-and-access-control), 
[database connections](../../access/databases/databases.md#connecting-to-database), [OpenAPI](../../access/open-api.md), and integrations with services like 
Benchling or Revvity Signals. Datagrok is developing a Plates application
to support predefined plate templates, batch ingestion, and
integrated analysis in one place.

##### <b>Q: How can I perform curve fitting, normalization, QC calculations from raw assay data in Datagrok?</b>

The [Curves](https://github.com/datagrok-ai/public/blob/master/packages/Curves/README.md)
plugin provides a complete workflow for converting raw assay data into fitted curves with QC
calculations. It transforms well-level assay data into fitted dose–response curves, supporting
functions such as 4PL and sigmoid, and can be easily extended to additional analyses, including
Km/Vmax, dose ratio, melt curves (DSF Tm), and qPCR, by defining the corresponding functions
and parameters.

For modeling dynamic systems, see the [Diff-Studio](../../compute/diff-studio.md) plugin that
allows users to solve sets of differential equations interactively through an intuitive UI and 
a declarative approach. 

##### <b>Q: Can Datagrok visualize dose-response curves for multiple compounds?</b> 

Datagrok supports visualizing DR curves for multiple compounds simultaneously and allows plotting curves using custom equations and parameters on the fly (e.g., 4-parameter curve fit over a 10-point dose-response).  
[Multicurve viewer](https://github.com/datagrok-ai/public/blob/master/packages/Curves/README.md#multi-curve-viewer)  | [Use case for assay plates and DRC](https://datagrok.ai/solutions/assay-plates)

## Deploy

##### <b>Q: How is Datagrok deployed?</b>

Datagrok supports Docker, Kubernetes, on-premises, and cloud deployments (AWS, GCP, Azure).  
Learn more about [Deployment](../../deploy/deploy.md)

##### <b>Q: What are the end-user computer and browser requirements?</b>

Datagrok runs completely in the browser. There is no
need to install any additional software. Datagrok is compatible with modern
browsers, including Chrome, Edge, and Safari.

##### <b>Q: How is logging and monitoring handled?</b>

Datagrok uses [AWS CloudWatch](https://aws.amazon.com/cloudwatch/) to collect and analyze logs and metrics.  

##### <b>Q: How are backups and restores managed?</b>

AWS provides automated RDS and S3 backups; RDS can also be restored as a
standard PostgreSQL database.  

##### <b>Q: What is the disaster recovery strategy?</b>

High availability is supported via Docker and AWS clusters that automatically
restart failed instances.  
See [Disaster recovery demo](https://www.youtube.com/watch?v=oFs9RShkHT8) for more information.

<!--## Deployment

#### **Q:** Do you have a guide for installing and deploying Datagrok on a GCP Kubernetes cluster?   
**A:** The fastest way to deploy Datagrok for evaluation is by using Docker Compose on a virtual machine. This setup takes just a minute or two.    
[Learn more](https://datagrok.ai/help/develop/admin/docker-compose)-->

## Access

##### <b>Q: What data sources can Datagrok connect to?</b>

Datagrok connects to any data source, including [databases](../../access/databases/connectors/connectors.md), [file storage systems](../../access/files/files.md), [web services and APIs](../../access/open-api.md). The platform supports [50+ file formats](../../access/files/supported-formats.md), including domain-specific like SDF, FASTA, and others.

##### <b>Q: Can I connect Datagrok to BigQuery or other data stores?</b>

Yes. Datagrok connects to [BigQuery](../../access/databases/connectors/bigquery.md) and most other [popular databases](../../access/databases/connectors/connectors.md) out of the box. In addition, any machine-readable data source can be easily integrated.

##### <b>Q: Does Datagrok support intuitive data query?</b>

Datagrok allows intuitive data queries using scientifically logical and searchable lists. Users can search across all data, including:  
- Logical conditions (e.g., `is not null`)  
- Combined logic with `AND` / `OR`  
- Date ranges and comparisons (`<`, `>`)  
- Text substring and fuzzy matching  

See [Data filtering and search](../../visualize/viewers/filters.md#search) for detailed guidance

## Govern

##### <b>Q: Do you follow secure development standards and industry best practices?</b>

Yes, we adhere to secure development and industry best practices across infrastructure, development, and enterprise security:

* Security-first infrastructure – designed with security from the ground up.
  Features include secure credentials management, flexible authentication
  (OAuth, SSO, Active Directory), and role-based access control.  
  [Learn more](../../deploy/GCP/deploy-gcp-gke-terraform.md) | [Authentication & authorization](../../develop/how-to/apps/build-an-app.md#authentication) | [Role-based access](../../compute/compute.md#privileges-and-visibility)

* Quality assurance – automated check for vulnerabilities (Snyk, Grype) and
  multi-layered testing including unit, integration, UI, and performance tests.  
  [Learn more](../../develop/qa/quality-assurance.md#continuous-integration-and-deployment-system) | [Testing](../../develop/qa/quality-assurance.md#automated-testing)

* Secure development lifecycle – CI/CD pipelines enforce automated build, test,
  and security checks; semantic versioning; and secure credential management.  
  [Learn more](../../develop/develop.md#continuous-integration) | [Version control](../../develop/dev-process/versioning-policy) | [Credential management](../../develop/how-to/packages/manage-credentials.md)

### Governance & access 

#### Authentication & Authorization 

##### <b>Q: How does Datagrok handle authentication and authorization?</b>

Datagrok uses role-based access control and integrates with enterprise identity
providers such as LDAP, SSO, and OAuth. 
[Learn more about access control](../../govern/access-control/access-control.md)

##### <b>Q: Who has access to Bitbucket, and what are the authentication requirements?</b>

Only datagrok core developers have access. MFA is being enabled to strengthen authentication. Bitbucket also provides built-in security scanning features.  

##### <b>Q: How is MFA enforced for GitHub accounts?</b>

GitHub requires all code contributors to enable two-factor authentication (2FA) as of March 2023. Developers comply by using device-tied passcodes.  

#### Data privacy & usage 

##### <b>Q: What happens to my data when I open a local file in Datagrok?</b>

When you open a local file in Datagrok (like dragging and dropping a file to your browser), you can analyze it without saving. This data stays in your browser's memory and isn't sent to the server unless you run resource-intensive server-side computations. Your data is gone when you close the browser tab. To save your work, you need to upload it to the server. Note that uploading data does not make it accessible to others. Your data stays private and visible to you only until you explicitly share it. Learn how to [save](../concepts/project/project.md#saving-entities-to-projects) and [share](../navigation/basic-tasks/basic-tasks.md#share) data.

##### <b>Q: What data or telemetry is sent back to Datagrok? What egress ports/protocols are used?</b>

* No automatic telemetry: Datagrok does not send any data to Datagrok servers by default.
* Optional error report/feedback: Users can optionally send feedback or error reports to Datagrok by selecting the "Send report to Datagrok" checkbox in the corresponding dialog.
* Package/image pulls: Datagrok can download images from Docker Hub or packages from NPM.
* Your connectors: Any external calls come from plugins you install (e.g., APIs, databases) over their standard ports/protocols.


### Data security 

#### Encryption 

##### <b>Q: Is data encrypted at rest?</b>

Yes, Datagrok relies on Amazon's built-in encryption for [RDS](https://docs.aws.amazon.com/AmazonRDS/latest/UserGuide/Overview.Encryption.html) and [S3 buckets](https://docs.aws.amazon.com/AmazonS3/latest/userguide/bucket-encryption.html).

##### <b>Q: Is data encrypted in transit?</b>

Yes, all client-server communications use [HTTPS](https://en.wikipedia.org/wiki/HTTPS), which means it is secure and encrypted. 

#### Vulnerability & patch management 

##### <b>Q: How does you manage vulnerabilities in the application and cloud infrastructure?</b>

We monitor CISA alerts and run daily [Snyk](https://snyk.io/) scans on container
builds. Vulnerabilities are triaged, remediated via infrastructure-as-code
pipelines, and verified through CI/CD testing. Customers are notified and
supported in version upgrades as needed. Learn more about
[Infrastructure](../../develop/under-the-hood/infrastructure.md)

##### <b>Q: How are OS and server patching handled?</b>

A: Jenkins and development servers are rebuilt quarterly with a fresh OS.
Documentation for builds is maintained in internal repositories. Learn more about [Deployment](../../deploy/deploy.md)

### Endpoint & device security 

#### Device controls 

##### <b>Q: What controls mitigate risks on BYOD devices?</b>

We have a policy that all devices must have active malware protection, updated
signatures, and drive encryption enabled.  

##### <b>Q: How is data exfiltration via removable media prevented?</b>

We have a policy that sensitive data must not be stored on removable media.
Passwords or credentials must not be transmitted unencrypted.  

#### EDR 

##### <b>Q: What is Datagrok’s current EDR approach?</b>

No dedicated EDR is deployed on associates’ devices. However, all devices must
have up-to-date and active malware protection. Associates are required to report
any suspicious activity via the alert channel. Incident response team triages
and executes the incident management process.  

##### <b>Q: How are security logs collected and monitored?</b>

A: AWS resources use centralized logging. Currently, there is no proactive log
review, but reporting of failed logins with alert thresholds is being
implemented.  

## Transform

##### <b>Q: Can I add a calculated column from a user interface?</b>

Yes. The [Add new column](../../transform/add-new-column.md) feature is very powerful. You can: 

* write your own expressions manually
* use 500+ available [functions](../../datagrok/concepts/functions/functions.md), or
* add your own functions using custom packages and scripts. 

For example, from a SMILES column you can generate molecular properties (e.g., MW, cLogP) or ADMET predictions.

## Visualize

##### <b>Q: What are the maximum dataset sizes?</b>

Datagrok handles millions of data points interactively for visualization and exploration. See [Why Datagrok?](../../datagrok/datagrok.md#why-datagrok) for details. 

##### <b>Q: What visualization options does Datagrok provide for data analysis and chemical structures?</b>  

Datagrok supports rich visualization for both chemical and general data, including:  

* Common chart types: [scatterplot](../../visualize/viewers/scatter-plot.md), [bar chart](../../visualize/viewers/bar-chart.md),  [pie chart](../../visualize/viewers/pie-chart.md), [box plot](../../visualize/viewers/box-plot.md)—[viewers](../../visualize/viewers/viewers.md)
* [Chemically aware viewers](../../datagrok/solutions/domains/chem/chemically-aware-viewers.md) and [Forms](../../visualize/viewers/forms.md)  
* Control charting for assay consistency: [line chart](../../visualize/viewers/line-chart.md), [statistical process control](../../visualize/viewers/line-chart.md#statistical-process-control), scatterplot with [formula lines](../../visualize/viewers/scatter-plot.md#formula-lines)  
* Statistical analysis support: [statistics viewer](../../visualize/viewers/statistics.md), [correlation plot](../../visualize/viewers/correlation-plot.md), scatter plot with [regression lines](../../visualize/viewers/scatter-plot.md#regression-lines), [statistical hypothesis testing](../navigation/views/table-view.md#statistical-hypothesis-testing)
* Advanced visual features: [coloring](../../visualize/viewers/grid.md#color-code-columns), [shaping and sizing](../../visualize/viewers/grid.md#resizing-columns), [formatting](../../visualize/viewers/grid.md#format-cells), [labeling](../../visualize/viewers/scatter-plot.md#labels), [trellising](../../visualize/viewers/trellis-plot.md)
* Table view integration for visualization, including [charts in cells](https://github.com/datagrok-ai/public/blob/master/help/visualize/viewers/charts-in-cells.md) 

See also: [Viewer gallery](https://github.com/datagrok-ai/public/blob/master/help/visualize/viewers/viewer-gallery.md), [Table view](../navigation/views/table-view.md), [Grid](../../visualize/viewers/grid.md).

##### <b>Q: Can I create custom visualizations?</b>

Yes. There're several options:

* Custom viewers via JavaScript API – 
[build custom viewers](../../develop/how-to/viewers/develop-custom-viewer.md) using the
[JavaScript API](../../develop/packages/js-api.md)
* Scripting viewers – use R, Python, or Julia to embed visualizations via
[scripting](../../compute/scripting/scripting.mdx).
* Third-party libraries – integrate frameworks like ECharts, D3.js, Circos, or three.js 
(see [Charts package](https://github.com/datagrok-ai/public/tree/master/packages/Charts))
* Custom file and folder viewers – extend Datagrok with 
file viewers for specific formats 
(e.g., PDB files via [NGL viewer](../../visualize/viewers/ngl.md) through the [Biostructure Viewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer#biostructure-viewer) package)
* [Custom cell renderers](../../develop/how-to/grid/custom-cell-renderers.md) - create custom visualization for cells in data [grid/table](../../visualize/viewers/grid.md).

##### <b>Q: How can I visualizing old and new data side by side?</b>

Multiple [viewer](../../visualize/viewers/viewers.md) and grids can be combined in a single [dashboard](../../datagrok/concepts/project/dashboard.md) to display different datasets side by side. 
Each viewer (scatterplots, charts, etc.) has configurable settings, allowing plots to source data from different tables. There are also multiple ways to [link tables](../../transform/link-tables.md), [join](../../transform/join-tables.md),  [pivot](../../transform/aggregate-rows.md), and apply other [transformations](../../transform/transform.md) 
to data within the intuitive UI. 

##### <b>Q: Can Datagrok tables be exported for presentations or analysis?</b>

Yes. Tables can be exported in multiple formats, including images, CSV/TSV, and SDF.  
Learn more about [Downloading](../navigation/basic-tasks/basic-tasks.md#download)

##### <b>Q: Can I merge multiple columns into a single "form cell" that displays label-value pairs (e.g., cLogP: 2.09) for these columns as a formatted list within that cell?</b>

Yes. Use [smart forms](../../visualize/viewers/grid.md#summary-columns).

<details>
<summary>See visual</summary>

![smart forms](../../deploy/releases/img/smart-forms.gif)

</details>

##### <b>Q: Can I pivot a table so each object (e.g., molecule) becomes a column and values (e.g., molecular properties) appear as rows?</b>

Yes. Use [Forms](../../visualize/viewers/forms.md) viewer.

<details>
<summary>See visual</summary>

![Forms viewer](../../visualize/viewers/img/forms.gif)

</details>

## Collaborate

##### <b>Q: Can users save and share custom dashboards, analyses, and visualizations in Datagrok?</b>  

Yes. Users can create projects that include data, analyses, and dashboards, then share them with others. Dashboards are dynamic and customizable.   
[Creating and managing projects](../concepts/project/project.md#creating-and-managing-projects) | [Creating dynamic dashboards](../../access/databases/databases.md#creating-dynamic-dashboards-for-query-results)  

##### <b>Q: How can I link raw data, metadata, and analysis results to compounds or sequences so that this information can be easily recalled and shared across Datagrok dashboards?</b>

The [Sticky meta](../../govern/catalog/sticky-meta.md) feature allows tagging
compounds, sequences, or other entities with meta-information, such as
calculated properties, analysis results, or comments. Tagged information can be
recalled anywhere in Datagrok when the same entity appears, supporting
collaborative workflows and effective communication between team members.  

## Develop

### Architecture

##### <b>Q: How is Datagrok architected?</b>

Datagrok is built as a modular, service-oriented platform with client-server separation.  
Learn more about [Datagrok's architecture](../../develop/under-the-hood/architecture.md#goals)

### Interoperability

##### <b>Q: Can Datagrok call external web services?</b>

Yes, both client and server can call authenticated web services and integrate
responses with the Datagrok API. See [Developer guide](../../develop/develop.md) for details

##### <b>Q: How do I interact with Datagrok from the server?</b>

Use Datagrok’s REST server API with authentication tokens.  
See [Proof of concept video](https://www.youtube.com/watch?v=TjApCwd_3hw) for more information

##### <b>Q: How do I interact with Datagrok from the client?</b>

Use the JavaScript API for custom applications and extensions.  
Learn more about  [JavaScript API](../../develop/packages/js-api.md)

### Developer experience

##### <b>Q: How do I develop and debug Datagrok customizations?</b>

Follow the [step-by-step guide](../../develop/how-to/packages/create-package.md) for creating your own package. If your package is based on a non-standard template, you may need to configure it differently for [debugging](../../develop/advanced/debugging.md). See also a [debug demo](https://youtu.be/PDcXLMsu6UM) for details.

##### <b>Q: How does the DevOps process work?</b>

The process includes package deployment, dependency management, and versioning.  
See [Developer guide](../../develop/develop.md) for details

##### <b>Q: Can multiple developers work concurrently?</b>

Yes, with versioning, branching, and package management.  
Learn more about [Concurrent development](../../develop/develop.md#development)

### Scalability and performance

##### <b>Q: How stable is the platform under heavy load?</b>

Datagrok is designed for horizontal scaling, ensuring stability with concurrent users.  
See [Scaling and stability](../../develop/under-the-hood/infrastructure.md#scalability) for more information

### Extensibility

##### <b>Q: Can I create custom visualizations?</b>

Yes, Datagrok offers multiple options for creating and embedding custom visualizations.  
See:
* [Custom viewers](../../develop/how-to/viewers/develop-custom-viewer.md)
* [Custom cell renderers](../../develop/how-to/grid/custom-cell-renderers.md)
* [Custom file viewers](../../develop/how-to/files/create-custom-file-viewers.md)

##### <b>Q: Can I build server-side components?</b>

Yes, custom back-end logic can be added.  
See [Admetica example](https://github.com/datagrok-ai/public/tree/master/packages/Admetica)

##### <b>Q: Can I add scripts and reuse them in components?</b>

Yes, Datagrok supports multiple scripting languages.  
Learn more about [Scripting](../../compute/scripting/scripting.mdx)

##### <b>Q: Can I reskin Datagrok?</b>

Yes, you can build tailored UI applications.  
See [Example app](https://public.datagrok.ai/apps/HitTriage/HitTriage?browse=apps)

##### <b>Q: Can I build full custom applications?</b>

Yes, including workflows, data models, state management, and persistence.  
See [Custom application packages](https://github.com/datagrok-ai/public/tree/master/packages)

### Frontend

##### <b>Q: What frontend capabilities does Datagrok provide for scientific data analysis?</b>

Datagrok delivers a high-performance, interactive frontend for handling, visualizing, and analyzing large-scale scientific datasets. Key capabilities include:

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



