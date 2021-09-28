# Compute

Provides support for scientific computations.

Scientific model is an arbitrary [function](https://datagrok.ai/help/overview/functions/function) with typed inputs and outputs.

* Inputs and outputs can be scalars, vectors, matrices
* Typically written in R, Python or Matlab/Octave and run on the server
* Other languages are possible, including compiled to WebAssembly and run in the browser
* Can be both simple and fast, or require considerable compute power to run, demanding
  automatic scaling
* Evaluated in multiple contexts of scientific process, such as variability analysis
  of output values in respect to input values variation
* Natural to compose, so that an output of one function becomes an input to the other
* The data source to the function input can be a CSV file, a table from a database query,
  a file from a file share or a call to a Web REST service
* The function itself and historical results of its runs are subject to sharing with
  collaborators and referencing in logs for decisions traceability

`Compute` package allows converting arbitrary functions in R, Python and any other language
into full-fledged scientific models, providing for easy to use and highly automated evaluation
and computation environment. In addition, Datagrok platform supports a layout markup, so that
the model function becomes a GUI-rich application with no manual coding.

## Roadmap

1. Stacked functions

2. Persistent, sharable historical runs

3. Sensitivity analysis

4. Modeling input parameters

5. Outlier detection and annotation

6. Model Repository and discoverability

7. Scalable and asyncronous computations

8. Export and reporting

9. Functions REST endpoints

10. Data annotation

11. Test data for functions

12. Functions versioning

13. Scaling on demand

Most of the above features are implemented in the package, some of them are part of the platform.