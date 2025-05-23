---
title: "Jupyter Notebook"
---

The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain live
code, equations, visualizations, and narrative text. Uses include data cleaning and transformation, numerical
simulation, statistical modeling, data visualization, machine learning, and much more.

Datagrok lets you create, edit, import, link, and apply Notebooks into tables.

![Jupyter Notebook](../uploads/gifs/jupyter-notebooks.gif "Jupyter Notebook")

## Create a notebook

There are several ways to create a new notebook:

1. From the **Top Menu**, select **ML** > **Notebooks** > **New Notebook..**  
2. Go to **Browse** > **Platform** > **Notebooks** and click **NEW NOTEBOOK** from the notebook toolbar.  
3. Run `CmdNewNotebook()` from the 
   [Console](../datagrok/navigation/panels/panels.md#console).  
4. Select a table and click **Actions** > **Open in Notebooks** in the 
   [Context Panel](../datagrok/navigation/panels/panels.md#context-panel).
5. Right-click the table's header and select **Table** > **Open in Notebook**.

If you want to link a notebook to a table, use one of the last two methods.

## Notebooks browser

In the [Notebooks Browser](https://public.datagrok.ai/notebooks), you can navigate over all available notebooks. If a notebook can be applied to any of the opened tables, an **Apply to** option will be available in its context menu.

## Apply existing notebooks into tables

To use one of the notebooks available in the platform, open a table and search for applicable notebooks in
the [Notebooks Browser](https://public.datagrok.ai/notebooks). Select a notebook and use its context menu to apply it to
your data.

## Return a table into the platform

Use simple function called "grok". Example:

```python
    grok(< table >)
```

where table in a Pandas dataframe.

## Supported language

* Python (3.8)

<!--Optionally: R, Julia and JavaScript.-->

## Environments

Each script can be run in an isolated environment, with predefined packages configuration, same as for
[Script Environments](scripting/advanced-scripting/under-the-hood.mdx#environment-isolation). The environment can be specified in the notebook properties.

## Importing notebooks

To import a notebook, just drag-and-drop the corresponding `.ipynb` file to the platform or use
[import](../access/files/files.md).

## Videos

[![Notebooks](../uploads/youtube/visualizations2.png "Open on Youtube")](https://www.youtube.com/watch?v=7MBXWzdC0-I&t=3880s)

See also:

* [https://jupyter.org](https://jupyter.org)
* [A gallery of interesting Jupyter notebooks](https://github.com/jupyter/jupyter/wiki/A-gallery-of-interesting-Jupyter-Notebooks#statistics-machine-learning-and-data-science)
