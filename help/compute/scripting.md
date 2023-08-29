---
title: "Scripting in Datagrok"
---

This article describes Datagrok scripting.
This is an extremely powerful concept allowing you to drastically improve
your data analysis experience with Datagrok.

Scripting combines fast
interactive visualizations and other features of the
Datagrok platform with thousands of statistical packages and
visualizations available in
[R](https://www.r-project.org/about.html), [Python](https://www.python.org),
[Octave](https://octave.org/), [Julia](https://julialang.org), or
[JavaScript](https://www.javascript.com).

## How scripting in Datagrok works?

Briefly, you only need to do the follows:

* Create a new script in the
  [Datagrok script editor](scripting-for-non-developers.mdx#working-with-datagrok-script-editor)
  or publish it as a part of a [package](../develop/packages/extensions).
* Write the code in any supported languages: Python, R, Julia, etc.
* Add **header**: set of annotation comments specifying the script language, tags,
  input and output variables.

That's all. These simple actions turn your script to a **Datagrok function**.
Without any additional effort you have received a wide range of possibilities:

* When you execute the script "as is",
  Datagrok automatically creates UI to specify input parameters and captures output data.
  The [RichFunctionView](scripting-advanced.mdx#using-richfunctionview-input-control)
  editors allows you to create complex UI with parameters tabs and groups
  just by adding annotation comments.
* Seamlessly transfer data between Datagrok and your script.
  For example, retrieve data table via Datagrok connections, process it by your custom script,
  return data table back and visualize it with any of the Datagrok viewers.
* Access computational history. Datagrok saves all runs in a history,
  so you can easily recall every run of the script.
* Sharing and collaborative work. You can easily share your scripts with your colleagues,
  specify access groups, share link to specific script run.
  See the [sharing documentaion](../collaborate/sharing) section for details.
* Use your script as input for other scripts.
  Probably, you want to process your data in Python and plot graphs in R? Easily.
  Just write two scripts and annotate parameters. Datagrok will care about all data transfers for you.
* Integrate your script in the Datagrok platform.
  For example, you can create a script plotting some graph for a molecule.
  By using [semantic types](../catalog/semantic-types)
  Datagrok recognizes the meaning of the data and automatically applies the script
  when you browse macromolecule details.

## Where to find how it works

* Review the [Scripting for non-professional developers](scripting-for-non-developers.mdx)
  for basic concepts of Datagrok scripting.
* For advanced features, review the [Advanced scription](scripting-advanced.mdx) section.
* For deep dive into Datagrok scripting, using JS api, controlling and developing your own components,
  review the [developer documentation](../develop/develop.md)
