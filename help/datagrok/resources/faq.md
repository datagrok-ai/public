---
title: FAQ
sidebar_position: 2
---

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