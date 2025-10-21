---
title: "Plugins"
sidebar_position: 1.2 
---

## General

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|---------------|-------|
|[Power Pack](https://github.com/datagrok-ai/public/tree/master/packages/PowerPack)| Platform | Recommended| Commonly used platform enhancements| Stable|
| [Tutorials](https://github.com/datagrok-ai/public/tree/master/packages/Tutorials)| Resources |Recommended| App for learning Datagrok with interactive tutorials, Demo app| Stable |

## Access

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|----------------------|---|
|[SQLite](https://github.com/datagrok-ai/public/tree/master/packages/SQLite)| Data format |Optional|Support for SQLite |Stable|
|[Arrow](https://github.com/datagrok-ai/public/tree/master/packages/Arrow)| Data format |Optional |Support for parquet and feather |Stable|
|[File Editors](https://github.com/datagrok-ai/public/tree/master/packages/FileEditors)| Data format |Optional |Support for PDF, DOCX |Stable|
|[Alation](https://github.com/datagrok-ai/public/tree/master/packages/Alation)| Misc |Optional| Integration with [Alation](https://www.alation.com/) (commercial 3rd party license)|Stable|
|[Jira Connect](https://github.com/datagrok-ai/public/tree/master/packages/JiraConnect)| Jira |Optional| Integration with Jira |Stable|

## Govern

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|---------------|-------------------|---|
| [Usage Analysis](https://github.com/datagrok-ai/public/tree/master/packages/UsageAnalysis) | Usage analysis | Optional|App that monitors platform and plugins usage, event log ([wiki](../govern/audit/usage-analysis.md)) | Beta |


## Visualize

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|---------------|-----|--------------|
|[Power Grid](https://github.com/datagrok-ai/public/tree/master/packages/PowerGrid)| Spreadsheet | Recommended| [Grid](../visualize/viewers/grid.md) viewer enhancements | Stable|
|[Charts](https://github.com/datagrok-ai/public/tree/master/packages/Charts)| Viewers, misc |Optional| Visualizations based on external libraries ([ECharts](https://echarts.apache.org/en/index.html), [D3.js](https://d3js.org/), [Circos](https://github.com/nicgirault/circosJS), and [three.js](https://threejs.org/))| Stable |
|[Dendrogram](https://github.com/datagrok-ai/public/tree/master/packages/Dendrogram)| Trees |Fit-for-purpose| [Dendrogram](../visualize/viewers/dendrogram.md) viewer | Stable |
|[PhyloTree Viewer](https://github.com/datagrok-ai/public/tree/master/packages/PhyloTreeViewer)| Trees | Fit-for-purpose | [PhylocanvasGL](../visualize/viewers/img/phylocanvas-gl-viewer.gif) viewer|  Stable |
|[GIS](https://github.com/datagrok-ai/public/tree/master/packages/GIS)| GIS |Fit-for-purpose|GIS (geographic information system) functionality, like maps and geocoding| Stable |
|[NMRium](https://github.com/datagrok-ai/public/tree/master/packages/nmrium)| NMR data |Fit-for-purpose| Visualization of NMR diagram | Stable |
|[Spectra Viewer](https://github.com/datagrok-ai/chem-spectra-viewer/tree/main)| Spectral data |Fit-for-purpose| Preview of JDX format for various types of spectra, including IR, UV, Mass, HPLC, absorption and others | Stable |
|[Excalidraw](https://github.com/datagrok-ai/public/tree/master/packages/Excalidraw)| Diagrams |Fit-for-purpose| Enables working with diagrams and drawings in the .excalidraw format| Stable |
<!--
|Viewers |[Forms]||???|Beta|
|Misc |[ChaRPy]|Fit-for-purpose|Adds two commands, "To Python script" and "To R script," to Datagrok viewers. These commands generate Python or R code for the selected viewer and execute the script to show the corresponding plot.|Beta|
-->

## Compute

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|----------------------|---|
|[Compute](https://github.com/datagrok-ai/public/tree/master/packages/Compute)| General |Required|Provides analytical and UI blocks for scientific computing|Stable|
|[Diff Studio](https://github.com/datagrok-ai/public/tree/master/packages/DiffStudio)| ODE Solver |Fit-for-purpose|An app for solving ordinary differential equations (ODE)|Stable|

## Scripting

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|----------------------|---|
|[Pyodide](https://github.com/datagrok-ai/public/tree/master/packages/Pyodide)| Scripting |Optional|Enables running Python scripts in the browser|Stable|

## Learn

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|----------------------|---|
|[EDA](https://github.com/datagrok-ai/public/tree/master/packages/EDA)| ML Toolkit |Required|ML toolkit: dimensionality reduction, multivariate analysis, supervised ML, ANOVA, etc.|Stable|
|[Notebooks](https://github.com/datagrok-ai/public/tree/master/packages/Notebooks)| Jupyter Notebooks |Optional|Integration with [JupyterLab Notebooks](https://jupyter.org/)|Stable|
|[MLflow](https://github.com/datagrok-ai/public/tree/master/packages/MLFlow)| MLflow models | Fit-for-purpose| Integrates MLflow models into the Datagrok platform|Beta|


## Develop

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|------------|----------------------|---|
|[Dev Tools](https://github.com/datagrok-ai/public/tree/master/packages/DevTools)| DevTools |Recommended|Developer tools (TestManager, DevPanel, etc.)|Stable|
|[API Samples](https://github.com/datagrok-ai/public/tree/master/packages/ApiSamples)| Resources |Recommended|Examples of Grok API|Stable|

## Solutions

### Chem

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|-------------------|---------------|---|
|[Chem](https://github.com/datagrok-ai/public/tree/master/packages/Chem) | General | Required| Cheminformatics support ([wiki](solutions/domains/chem/chem.md)).<br/>Comes with OpenChemLib sketcher | Stable|
|[Marvin](https://github.com/datagrok-ai/labs/tree/master/packages/Marvin) | Sketchers | Optional| Integration with Marvin JS (commercial 3rd party license) | Stable|
|[Ketcher Sketcher](https://github.com/datagrok-ai/public/tree/master/packages/KetcherSketcher) | Sketchers | Optional| Integration with [Ketcher](https://lifescience.opensource.epam.com/ketcher/index.html) (Apache license, version 2.0) | Stable|
|[Chem Draw Sketcher](https://github.com/datagrok-ai/labs/tree/master/packages/ChemDraw) | Sketchers | Optional| Integration with ChemDraw (commercial 3rd party license) | Stable|
|[ChEMBL](https://github.com/datagrok-ai/public/tree/master/packages/Chembl) | Database search | Recommended| ChEMBL database for on-premise deployment | Stable|
|[DrugBank](https://github.com/datagrok-ai/public/tree/master/packages/DrugBank) | Database search | Optional| Information on 11,300 drugs from [DrugBank](https://go.drugbank.com/)| Stable|
|[ChEMBL API](https://github.com/datagrok-ai/public/tree/master/packages/ChemblAPI) | Database search | Optional| Webservice integration: ChEMBL | Stable|
|[PubChem](https://github.com/datagrok-ai/public/tree/master/packages/PubChemApi) | Database search | Optional| Webservice integration: [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Stable|
|[Chemspace](https://github.com/datagrok-ai/public/tree/master/packages/Chemspace) | Database search | Optional| Webservice integration: [Chemspace](https://chem-space.com/)  |  Stable|
|[SureChEMBL](https://github.com/datagrok-ai/public/tree/master/packages/SureChembl)| Database search | Optional| Performs searches through a locally deployed [SureChEMBL](https://www.surechembl.org) database | Beta|
| [Docking](https://github.com/datagrok-ai/public/tree/master/packages/Docking)| Virtual screening | Optional |Let's you batch screen libraries against AutoDock-prepared targets with interactive visualization | Stable |
|[Hit Triage](https://github.com/datagrok-ai/public/tree/master/packages/HitTriage)| Virtual screening<br/>Hit to lead | Fit-for-purpose| Apps for virtual screening (Hit Triage) and hit design (Hit Design) | Stable|
| [Retrosynthesis](https://github.com/datagrok-ai/public/tree/master/packages/Retrosynthesis)| Synthetic planning | Optional |Creates retrosynthesis paths for the selected molecule, built on top of [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder) | Beta |
| [CDD Vault Link](https://github.com/datagrok-ai/public/tree/master/packages/CddVaultLink)| Registration system integration | Optional |Provides integration with [CDD Vault](https://www.collaborativedrug.com/cdd-informatics-platform) registration system | Beta |
| [Reinvent4](https://github.com/datagrok-ai/public/tree/master/packages/Reinvent4)| Generative molecular design | Fit-for-purpose |Provides integration with the [Reinvent4](https://github.com/MolecularAI/REINVENT4) molecular design tool | Beta |

### Bio

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|-----------------|-----------------|---|
|[Bio](https://github.com/datagrok-ai/public/tree/master/packages/Bio) | General | Required| Bioinformatics support ([wiki](solutions/domains/bio/bio.md))|Stable|
|[HELM](https://github.com/datagrok-ai/public/tree/master/packages/Helm)| Data format |Recommended |Support for [HELM notation](https://pistoiaalliance.atlassian.net/wiki/spaces/HELM/overview), HELM editor |Stable|
|[Biostructure Viewer](https://github.com/datagrok-ai/public/tree/master/packages/BiostructureViewer) | Visualization | Recommended| Visualization of biological structures ([wiki](../visualize/viewers/biostructure.md)) |Stable|
| [Peptides](https://github.com/datagrok-ai/public/tree/master/packages/Peptides)| SAR  | Fit-for-purpose| App for sequence-activity relationship analysis for peptides ([wiki](solutions/domains/bio/peptides-sar.md))|Stable|
|[Curves](https://github.com/datagrok-ai/public/tree/master/packages/Curves)| Curve fitting | Fit-for-purpose| Support for fitted curves (like dose-response curves), including in-grid rendering, storing charts in cells, interactivity, and automatic fitting | Stable |
|[Sequence Translator](https://github.com/datagrok-ai/public/tree/master/packages/SequenceTranslator) | Oligonucleotides  | Fit-for-purpose| App that converts oligonucleotides into various formats  |Stable|
|[Oligo Batch Calculator](https://github.com/datagrok-ai/public/tree/master/packages/OligoBatchCalculator) | Oligonucleotides  | Fit-for-purpose| App that calculates oligonucleotide properties |Beta|
| [Admetica](https://github.com/datagrok-ai/public/tree/master/packages/Admetica) | ADMET | Optional | App that lets you evaluate ADMET properties | Stable |
| [Docking](https://github.com/datagrok-ai/public/tree/master/packages/Docking) | Docking | Optional | Integration with [Autodock GPU](https://catalog.ngc.nvidia.com/orgs/hpc/containers/autodock) that let's you run docking and analyze the results in Datagrok| Stable |
| [BioNeMo](https://github.com/datagrok-ai/public/tree/master/packages/BioNeMo) | Protein structure prediction, docking | Fit-for-purpose | Integrates advanced models for protein structure prediction and molecular docking| Beta |

### NLP

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|---------------|-------------------|---|
|[NLP](https://github.com/datagrok-ai/public/tree/master/packages/NLP)| NLP |Optional|Integration with AWS Translate, a neural machine translation service. Extends Datagrok with info panels for text files|Stable|

### Digital health, clinical

|Plugin <div style={{ width:140 }}></div> |Area <div style={{ width:130 }}></div>|  Tag <div style={{ width:110 }}></div> |Description <div style={{ width:315 }}></div> | Release|
|------|-------|---------------|-------------------|---|
| [Clinical Case](https://github.com/datagrok-ai/public/tree/master/packages/ClinicalCase) | Clinical data | Fit-for-purpose|App for analyzing clinical data in SDTM format| Demo |
