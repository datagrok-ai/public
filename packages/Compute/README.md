# Compute

Provides support for scientific computations.

Scientific model is an arbitrary [function](https://datagrok.ai/help/overview/functions/function)
with typed inputs and outputs:

* Inputs and outputs can be scalars, vectors, matrices
* Typically written in R, Python or Matlab/Octave and run on the server
* Other languages are possible, including compiled to WebAssembly and run in the browser
* Can be both simple and fast, or require considerable compute power to run, demanding
  automatic scaling
* Evaluated in multiple contexts of experiment design, such as variability analysis
  of output values with respect to input values
* Natural to compose, so that an output of one function becomes an input to the other
* A data source to the function input can be a CSV file, a table from a database query,
  a file from a file share or a return from a call to a Web REST service
* The function itself and historical results of its runs are subject to sharing with
  collaborators and referencing in logs for decisions traceability

`Compute` package:

* Allows using arbitrary functions in R, Python and any other language as full-fledged
  scientific models
* Enables a full lifecycle of mathematical models: create, tune, share, use, validate
* Provides for easy to use and highly automated evaluation and computation environment
* Enables contextual process for Design of Experiment (sensitivity analysis and more)

In addition, Datagrok platform supports a UI layout markup, so that the model function becomes
a GUI-rich application with no manual coding.

Project planning board: [link](https://github.com/datagrok-ai/public/projects/8).

## Roadmap

1. **UI markup**

    Annotate function inputs and outputs to produce highly interactive, visually reach GUI:
    
    * arrange inputs and outputs in blocks and tabs
    * add captions, units of measure and other information to inputs and outputs
    * automatically produce additional plots with Datagrok
    [viewers](https://datagrok.ai/help/visualize/viewers)

2. **Input providers**

    Produce inputs to functions in-place as outputs of other functions (aka input providers),
    including:
    
    * queries to databases
    * dialog-based functions (outlier detection, data annotation)
    * queries to OpenAPI and REST endpoints
    * other computing functions with or without GUI

    These may include UI parts as well. The input provider is specified as part of the Universal UI
    markup.

3. **Persistent, sharable historical runs**

    It is already possible to provide a link to a function (with specified input parameters in
    the URI), which will open a function view and run it.
    
    Once a certain version of a specific function is run with specific inputs, the result should
    be stored in the immutable database log along with the inputs. Later it will be used to verify
    the grounds for decisions made from these calculations.

4. **Sensitivity analysis**

    * Sample inputs:
      * by specified number of samples
      * by a specified distribution or within a range
      * for a specified set of scalar inputs and/or columns of the matrix input
    * Produce variability analysis for outputs based on the sampled inputs
    * Visualize the results of analysis with Datagrok viewers

5. **Modeling input parameters**

    Solve an inversion problem: identify input conditions leading to specified output constraints.

6. **Outlier detection and annotation**

    * Automatic outliers detection
    * Manual outliers markup and annotation
    * Used as an input provider in other functions

7. **Model Repository and discoverability**

8. **Scaling on demand**

9. **Asynchronous computation**

10. **Export and reporting**

11. **REST endpoints**

12. **Data annotation**

13. **Test data for functions**

14. **Functions versioning**

15. **Audit**

Most of the features are implemented in this package, some of them are part of the platform.