<!-- TITLE: Datagrok for life sciences -->
<!-- SUBTITLE: -->

# Datagrok for life sciences

Datagrok is a great tool for working with tabular data of any origin. It will retrieve, analyze, visualize and transform
traditional business data, such as sales, inventory, or marketing data as well as any other system (but of course
faster). However, it really shines in dealing with scientific data.

"Scientific data" is a very broad term, but there are certain attributes that quite often apply to it. Let's take a look
at them and see how Datagrok shows itself.

## Volume

Due to the nature of being measured by instruments, the volume of the data can greatly exceed the typical "business"
data.

While no single computer is capable of working with datasets of unlimited size, Datagrok has been engineered to
efficiently work with as big datasets as possible. It offers two modes:

* **Local data**—the whole dataset is loaded into the browser memory. Thanks to our proprietary
  [in-memory database](../develop/advanced/performance.md#in-memory-database), the platform can efficiently work with
  billions of rows or millions of columns. This is the preferred way of working with data since it offers interactive
  visualizations and instant access to all underlying data points. Most of the datasets fall into that category.
* **Remote data**—sometimes, it is impractical or impossible to download the whole dataset, but you still need to work
  with it. In this case, Datagrok offers several solutions that will transparently translate the necessary work to be
  performed on the server-side. Here are some of them:
    * [DB Exploration](../access/db-exploration.md) to explore database schemas
    * [DB Visual Query](../access/db-visual-query.md) to aggregate and pivot on the database side

## Complexity

Scientific data has a complex, multidimensional structure. Fortunately, the platform is well-equipped to handle data of
diverse complexity:

* Unstructured or semi-structured data
* Graph data
* Instrument data
* Incomplete data

## Data discovery

Given the volume, variety, and complexity of the scientific data available to companies, being able to efficiently
discover relevant datasets is essential.

## Rendering

Unlike business data, more often than not raw values in scientific datasets can not be interpreted by simply looking at
them. Most common reasons for that are:

1. A value is encoded. Examples: molecules encoded as SMILES
2. A value can be interpreted, but only makes sense as part of the bigger picture. Examples: a particular pixel of an
   image
3. A value makes sense only together with other attributes. Examples: (LAT, LNG) points
4. A value can be interpreted only in the context of the continuous stream of data. Examples: ECG values

In order to present users with a meaningful interpretation of a dataset, the platform tries to infer
the [true meaning of the raw data](../discover/semantic-types.md) and then offers a number of visualizations and actions
that can be applied either to the whole dataset or an individual point.

Examples:

* [rendering molecules]
* [lat/lng]
* [address]

## Metadata

Due to the complexity of scientific data, having the proper metadata associated with it is pivotal for computers (as
well as people!) to understand it. Datagrok provides a comprehensive
[metadata management framework](../discover/metadata.md).

We also support [FAIR Principles](../discover/fair.md) where applicable.

## Formats

Scientific data often comes in formats that are not understood by a typical business application. In addition to
commonly used data formats, such as csv, txt, xsls, xml, and json, Datagrok
[supports](../access/importing-data.md#supported-file-types) the following ones:

| Extension     | Description          |
|---------------|----------------------|
| edf           | European Data Format |
| sas7bdat      | SAS                  |
| rds           | R RDS                |
| rda           | R RDA                |
| h5            | H5                   |
| nc            | NetCDF               |
| mat           | MATLAB MAT           |
| ipynb         | Jupyter Notebook     |

## Reproducibility

* [Jupyter notebooks](../compute/jupyter-notebook.md)

## Integration

* LIMS systems

## References

* [FAIR principles](https://www.go-fair.org/fair-principles/)
* [Software Engineering for Machine Learning: A Case Study](https://www.microsoft.com/en-us/research/publication/software-engineering-for-machine-learning-a-case-study/)
* [Insight into complex scientific data using a graph data store](https://medium.com/blackfynn/insight-into-complex-scientific-data-using-a-graph-data-store-f2b540684c84)
