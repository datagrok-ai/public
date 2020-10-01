<!-- TITLE: Jupyter Notebook -->
<!-- SUBTITLE: -->

# Jupyter Notebook

The Jupyter Notebook is an open-source web application that allows you to create and 
share documents that contain live code, equations, visualizations and narrative text. 
Uses include: data cleaning and transformation, numerical simulation, statistical modeling, 
data visualization, machine learning, and much more.

Datagrok allows to create, edit, import, link and apply Notebooks into tables.

![Jupyter Notebook](../uploads/gifs/jupyter-notebooks.gif "Jupyter Notebook")

## Create a Notebook

Creating a new notebook is very easy. There are several ways to do this:

  1. Click `Functions | Notebooks | New Notebook`
  2. Run #{x.CmdNewNotebook} from the [Console](https://datagrok.ai/help/overview/console)
  3. Select one or more tables and click `Actions | Open in Notebooks` on the [Property Panel](../overview/property-panel.md)

If you want to link a notebook to one or more tables, please use the third method.

## Notebooks Browser

In the [Notebooks Browser](https://public.datagrok.ai/notebooks), you can navigate over all 
available notebooks. If any of notebooks can be applicable to any of opened tables, it will be marked with 
"Applicable to" note. Also filtering can be used to filter only applicable notebooks.

## Apply Existing Notebooks into Tables

To use one of the notebooks available in the platform, open a table and search for applicable notebooks in the [Notebooks Browser](https://public.datagrok.ai/notebooks). Select a notebook and use its context menu to apply it to your data.

## Return a Table into the Platform

Use simple function called "grok". 
Example:
```python
    grok(<table>)
```
where table in a Pandas dataframe.
  
## Supported Languages

* Python (3.7)

Optionally: R, Julia and JavaScript.

## Environments

Each script can be run in an isolated environment, with predefined packages configuration, same as for 
[Script Environments](../compute/scripting.md#Environments). The environment can be specified in the notebook properties. 

# Importing Notebooks

To import a notebook, just drag-and-drop the corresponding `.ipynb` file to the platform or use 
[import](../access/importing-data.md).

See also:

  * [http://jupyter.org](http://jupyter.org)
  * [A gallery of interesting Jupyter notebooks](https://github.com/jupyter/jupyter/wiki/A-gallery-of-interesting-Jupyter-Notebooks#statistics-machine-learning-and-data-science)
