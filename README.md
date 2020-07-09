# Datagrok package repository

This is a public repository for [packages](https://datagrok.ai/help/develop/develop#packages) 
available for the [Datagrok™](https://datagrok.ai) platform, 
curated by the Datagrok team. It contains different scientific methods integrated with the platform, 
as well as custom viewers, application, and functions.

Organizations that deploy Datagrok on their premises can access public packages. In addition to that,
they might establish their own private repositories (most likely on their premises on the virtual private
cloud) that would contain their proprietary applications built on top of the platform. 

# Ideas for the contributions

If you want to get familiar with the platform, here are some ideas. Pick whatever interests you,
and reach out to Andrew (askalkin@datagrok.ai) or post on our [community forum](https://community.datagrok.ai/). 

* Visualizations
  * Sankey diagram
  * Timeline chart
  * Gantt chart
  * Sunburst chart
  * Radar chart
  * Port visjs-based [network diagram](https://datagrok.ai/help/visualize/viewers/network-diagram) from Dart to JavaScript
  * WebGL-based rendering of the 2D scatter plot to work with 10M+ points
  * GIS: leaflet-m
  * [Event drops](https://github.com/marmelab/EventDrops)
* Scientific methods
  * Statistical hypothesis testing
  * Bayesian statistics
  * Computer vision
  * NLP 
    * integration with the modern models (BERT, etc)
    * pre-defined tasks (sentiment extraction, etc)
* Apps  
  * Examples: Real-time sentiment analysis of tweets 
  * Examples: Amazon store clone 
* File editors and viewers
* File metadata extractors (see Apache Kafka)
* WASM-based support for digital signal processing
* Domain-specific algorithms  
* Connectors to web services
* Bioinformatics
* Telecom
* Fintech

# Package Development Guide

Here are some of the best practices for developing and publishing Datagrok plugins: 

1. [Presentation](#presentation)
    * [Repository](#repo)
    * Name
    * Description
    * Links
    * Metadata
2. [Code](#code)
    * File Structure
    * Plugin API
3. Publishing  

## Presentation

### Repository

## Code



### Plugin API

See also:
* [Datagrok overview](https://datagrok.ai/)
* [JavaScript development](https://datagrok.ai/help/develop/develop)