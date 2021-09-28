# Compute

Provides support for scientific computations.

Scientific model is an arbitrary function with typed inputs and outputs.

* Usually written in R or Python, it can be both simple and fast, or require considerable time
  to run, which may demand automatic scaling
* Input types can be scalars, vectors and matrices
* Such function is subject to evaluation in multiple contexts, essential for scientific studies,
  such as variability analysis of the output values in respect to input values variation
* It is natural to compose such functions, so that an output of one becomes an input to the other
* The data source for the function input can be a CSV file, a table from a database query, or
  a file from a file share or an Internet address
* Both the function itself and historical results of its runs are subject to sharing with
  collaborators and referencing in logs for traceability of decisions based on them

`Compute` package allows converting arbitrary scripts in R, Python and any other language,
into full-fledged scientific models, providing for easy to use evaluation and computation
environment. In addition, Datagrok platform supports a script markup, so that the script
becomes a GUI-rich application with no manual coding.