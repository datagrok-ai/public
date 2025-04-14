---
title: "Plugins"
sidebar_position: 1.2 
---

## General

|Area <div style={{ width:130 }}></div>  |Plugin <div style={{ width:140 }}></div>|  Tag <div style={{ width:115 }}></div> |Description  | Release|
|------|-------|------------|---------------|-------|
| Platform |[Power Pack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack)| Recommended| Commonly used platform enhancements| Stable|
| Resources | [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials)|Recommended| App for learning Datagrok with interactive tutorials, Demo app| Stable |

## Access

|Area <div style={{ width:130 }}></div> |Plugin<div style={{ width:140 }}></div>|  Tag <div style={{ width:115 }}></div> |Description <div style={{ width:270 }}></div> | Release|
|------|-------|------------|----------------------|---|
| Data format |[SQLite](https://github.com/datagrok-ai/public/tree/master/packages/SQLite)|Optional|Support for SQLite |Stable|
| Data format |[Arrow](https://github.com/datagrok-ai/public/tree/master/packages/Arrow)|Optional |Support for parquet and feather |Stable|
| Data format |[File Editors](https://github.com/datagrok-ai/public/tree/master/packages/FileEditors)|Optional |Support for PDF, DOCX |Stable|
| Misc |[Alation](https://github.com/datagrok-ai/public/tree/master/packages/Alation)|Optional| Integration with [Alation](https://www.alation.com/) (commercial 3rd party license)|Stable|
| Jira |[Jira Connect](https://github.com/datagrok-ai/public/tree/master/packages/JiraConnect)|Optional| Integration with Jira |Stable|

## Govern

|Area <div style={{ width:130 }}></div> |Plugin<div style={{ width:140 }}></div>|  Tag <div style={{ width:90 }}></div> |Description  | Release|
|------|-------|---------------|-------------------|---|
| Usage analysis | [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis) | Optional|App that monitors platform and plugins usage, event log ([wiki](../govern/audit/usage-analysis.md)) | Beta |


## Visualize

|Area <div style={{ width:130 }}></div>  |Plugin<div style={{ width:140 }}></div>|  Tag <div style={{ width:110 }}></div> |Description | Release  |
|------|-------|---------------|-----|--------------|
|Spreadsheet |[Power Grid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid)| Recommended| [Grid](../visualize/viewers/grid.md) viewer enhancements | Stable|
|Viewers, misc |[Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts)|Optional| Visualizations based on external libraries ([ECharts](https://echarts.apache.org/en/index.html), [D3.js](https://d3js.org/), [Circos](https://github.com/nicgirault/circosJS), and [three.js](https://threejs.org/))| Stable |
| Trees |[Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram)|Fit-for-purpose| [Dendrogram](../visualize/viewers/dendrogram.md) viewer | Stable |
| Trees |[PhyloTree Viewer](https://github.com/datagrok-ai/public/tree/master/packages/PhyloTreeViewer)| Fit-for-purpose | [PhylocanvasGL](../visualize/viewers/img/phylocanvas-gl-viewer.gif) viewer|  Stable |
|GIS|[GIS](https://github.com/datagrok-ai/public/tree/master/packages/GIS)|Fit-for-purpose|GIS (geographic information system) functionality, like maps and geocoding| Stable |
|NMR data |[NMRium](https://github.com/datagrok-ai/public/tree/master/packages/nmrium)|Fit-for-purpose| Visualization of NMR diagram | Stable |
|Spectral data |[Spectra Viewer](https://github.com/datagrok-ai/chem-spectra-viewer/tree/main)|Fit-for-purpose| Preview of JDX format for various types of spectra, including IR, UV, Mass, HPLC, absorption and others | Stable |
<!--
|Viewers |[Forms]||???|Beta|
|Misc |[ChaRPy]|Fit-for-purpose|Adds two commands, "To Python script" and "To R script," to Datagrok viewers. These commands generate Python or R code for the selected viewer and execute the script to show the corresponding plot.|Beta|
-->

## Compute

|Area <div style={{ width:125 }}></div> |Plugin| Tag <div style={{ width:110 }}></div> |Description  | Release|
|------|-------|------------|----------------------|---|
|General |[Compute](https://github.com/datagrok-ai/public/tree/master/packages/Compute)|Required|Provides analytical and UI blocks for scientific computing|Stable|
|ODE Solver |[Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/ODEs)|Fit-for-purpose|An app for solving ordinary differential equations (ODE)|Stable|

## Scripting

|Area <div style={{ width:125 }}></div> |Plugin| Tag <div style={{ width:110 }}></div> |Description  | Release|
|------|-------|------------|----------------------|---|
|Scripting |[Pyodide](https://github.com/datagrok-ai/public/tree/master/packages/Pyodide)|Optional|Enables running Python scripts in the browser|Stable|
<!----

|Modeling |[Bioreactors](https://github.com/datagrok-ai/public/tree/master/packages/Bioreactors)|Fit-for-purpose|Simulation of the Controlled Fab-Arm Exchange mechanism|Alpha|
|Modeling |SimPKPD|Fit-for-purpose|App for PKPD simulations|Labs|

---->

## Learn

|Area <div style={{ width:125 }}></div> |Plugin <div style={{ width:140 }}></div>| Tag <div style={{ width:110 }}></div> |Description  | Release|
|------|-------|------------|----------------------|---|
|ML Toolkit |[EDA](https://github.com/datagrok-ai/public/tree/master/packages/EDA)|Required|ML toolkit: dimensionality reduction, multivariate analysis, supervised ML, ANOVA, etc.|Stable|
|Jupyter Notebooks |[Notebooks](https://github.com/datagrok-ai/public/tree/master/packages/Notebooks)|Optional|Integration with [JupyterLab Notebooks](https://jupyter.org/)|Stable|
|MLflow models |[MLflow](https://github.com/datagrok-ai/public/tree/master/packages/MLFlow)| Fit-for-purpose| Integrates MLflow models into the Datagrok platform|Beta|


## Develop

|Area <div style={{ width:125 }}></div> |Plugin <div style={{ width:140 }}></div>| Tag <div style={{ width:110 }}></div> |Description  | Release|
|------|-------|------------|----------------------|---|
|DevTools |[Dev Tools](https://github.com/datagrok-ai/public/tree/master/packages/DevTools)|Recommended|Developer tools (TestManager, DevPanel, etc.)|Stable|
|Resources |[API Samples](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples)|Recommended|Examples of Grok API|Stable|

## Solutions

### Chem

|Area <div style={{ width:125 }}></div> |Plugin| Tag <div style={{ width:110 }}></div> |Description  | Release|
|------|-------|-------------------|---------------|---|
| General |[Chem](https://github.com/datagrok-ai/public/tree/master/packages/Chem) | Required| Cheminformatics support ([wiki](solutions/domains/chem/chem.md)).<br/>Comes with OpenChemLib sketcher | Stable|
| Sketchers |[Marvin](https://github.com/datagrok-ai/labs/tree/master/packages/Marvin) | Optional| Integration with Marvin JS (commercial 3rd party license) | Stable|
| Sketchers |[Ketcher Sketcher](https://github.com/datagrok-ai/public/tree/master/packages/KetcherSketcher) | Optional| Integration with [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html) (Apache license, version 2.0) | Stable|
| Sketchers |[Chem Draw Sketcher](https://github.com/datagrok-ai/labs/tree/master/packages/ChemDraw) | Optional| Integration with ChemDraw (commercial 3rd party license) | Stable|
| Database search |[ChEMBL](https://github.com/datagrok-ai/public/tree/master/packages/Chembl) | Recommended| ChEMBL database for on-premise deployment | Stable|
| Database search |[DrugBank](https://github.com/datagrok-ai/public/tree/master/packages/DrugBank)   | Optional| Information on 11,300 drugs from [DrugBank](https://go.drugbank.com/)| Stable|
| Database search |[ChEMBL API](https://github.com/datagrok-ai/public/tree/master/packages/ChemblAPI) | Optional| Webservice integration: ChEMBL | Stable|
| Database search |[PubChem](https://github.com/datagrok-ai/public/tree/master/packages/PubChemApi)   | Optional| Webservice integration: [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Stable|
| Database search |[Chemspace](https://github.com/datagrok-ai/public/tree/master/packages/Chemspace) | Optional| Webservice integration: [Chemspace](https://chem-space.com/)  |  Stable|
| Database search |[SureChEMBL](https://github.com/datagrok-ai/public/tree/master/packages/SureChEMBL)| Optional| Performs searches through a locally deployed [SureChEMBL](https://www.surechembl.org) database | Beta|
| Virtual screening | [Docking](https://github.com/datagrok-ai/public/tree/master/packages/Docking)| Optional |Let's you batch screen libraries against AutoDock-prepared targets with interactive visualization | Stable |
| Virtual screening<br/>Hit to lead |[Hit Triage](https://github.com/datagrok-ai/public/tree/master/packages/HitTriage)| Fit-for-purpose| Apps for virtual screening (Hit Triage) and hit design (Hit Design) | Stable|
| Synthetic planning | [Retrosynthesis](https://github.com/datagrok-ai/public/tree/master/packages/Retrosynthesis)| Optional |Creates retrosynthesis paths for the selected molecule | Beta |
| Registration system integration | [CDD Vault Link](https://github.com/datagrok-ai/public/tree/master/packages/CddVaultLink)| Optional |Provides integration with [CDD Vault](https://www.collaborativedrug.com/cdd-informatics-platform) registration system | Beta |
<!--
| Misc |[EnamineStore](https://github.com/datagrok-ai/public/tree/master/packages/EnamineStore) | Optional| Webservice integration: [Enamine](https://enaminestore.com/search), a service for online shopping for chemical building blocks |  Stable|


|Misc |[Chemspace](https://github.com/datagrok-ai/public/tree/master/packages/Chemspace)|Misc|Integration with the Chemspace, a service for online shopping for the chemical building blocks|Alpha|
-->

### Bio

|Area  |Plugin|Tag <div style={{ width:110 }}></div> | Description | Release|
|------|-------|-----------------|-----------------|---|
| General |[Bio](https://github.com/datagrok-ai/public/tree/master/packages/Bio) | Required| Bioinformatics support ([wiki](solutions/domains/bio/bio.md))|Stable|
|Data format|[HELM](https://github.com/datagrok-ai/public/tree/master/packages/HELM)|Recommended |Support for [HELM notation](https://pistoiaalliance.atlassian.net/wiki/spaces/HELM/overview), HELM editor |Stable|
| Visualization |[Biostructure Viewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer) | Recommended| Visualization of biological structures ([wiki](../visualize/viewers/biostructure.md)) |Stable|
| SAR  | [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides)| Fit-for-purpose| App for sequence-activity relationship analysis for peptides ([wiki](solutions/domains/bio/peptides-sar.md))|Stable|
| Curve fitting |[Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves)| Fit-for-purpose| Support for fitted curves (like dose-response curves), including in-grid rendering, storing charts in cells, interactivity, and automatic fitting | Stable |
| Oligonucleotides  |[Sequence Translator](https://github.com/datagrok-ai/public/tree/master/packages/SequenceTranslator) | Fit-for-purpose| App that converts oligonucleotides into various formats  |Stable|
| Oligonucleotides  |[Oligo Batch Calculator](https://github.com/datagrok-ai/public/tree/master/packages/OligoBatchCalculator) |Fit-for-purpose| App that calculates oligonucleotide properties |Beta|
| ADMET | [Admetica](https://github.com/datagrok-ai/public/tree/master/packages/Admetica) | Optional | App that lets you evaluate ADMET properties | Stable |
| Docking | [Docking](https://github.com/datagrok-ai/public/tree/master/packages/Docking) | Optional | Integration with [Autodock GPU](https://catalog.ngc.nvidia.com/orgs/hpc/containers/autodock) that let's you run docking and analyze the results in Datagrok| Stable |
|Protein structure prediction, docking | [BioNeMo](https://github.com/datagrok-ai/public/tree/master/packages/BioNeMo) | Fit-for-purpose | Integrates advanced models for protein structure prediction and molecular docking| Beta |

### NLP

|Area <div style={{ width:110 }}></div> |Plugin <div style={{ width:140 }}></div>|  Tag <div style={{ width:110 }}></div> |Description   | Release|
|------|-------|---------------|-------------------|---|
|NLP |[NLP](https://github.com/datagrok-ai/public/tree/master/packages/NLP)|Optional|Integration with AWS Translate, a neural machine translation service. Extends Datagrok with info panels for text files|Stable|

### Digital health, clinical

|Area <div style={{ width:110 }}></div> |Plugin <div style={{ width:140 }}></div>|  Tag <div style={{ width:110 }}></div> |Description   | Release|
|------|-------|---------------|-------------------|---|
| Clinical data | [Clinical Case](https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase) | Fit-for-purpose|App for analyzing clinical data in SDTM format| Demo |
