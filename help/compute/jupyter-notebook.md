<!-- TITLE: Jupyter Notebook -->
<!-- SUBTITLE: -->

# Jupyter Notebook

The Jupyter Notebook is an open-source web application that allows you to create and 
share documents that contain live code, equations, visualizations and narrative text. 
Uses include: data cleaning and transformation, numerical simulation, statistical modeling, 
data visualization, machine learning, and much more.

Datagrok allows to create, edit, import, link and apply Notebooks into tables.

![Jupyter Notebook](../uploads/gifs/jupyter-notebooks.gif "Jupyter Notebook")

## Create Notebook

Creating a new Notebook is very easy. There some ways to do this:
1. Click "Tools | Notebooks | New Notebook"
2. Run #{x.CmdNewNotebook} from console
3. Select one or more tables and click on [Property Panel](../features/property-panel.md) "Actions | Open in Notebooks"
Note: to link the Notebook with table(s) you should use way "3".

## Notebooks Browser

To open Notebooks Browser click "Tools | Notebooks | Browse Notebooks". There are you can navigate over all 
available notebooks. If any of notebooks can be applicable to any of opened tables, it will be marked with 
"Applicable to" note. Also filtering can be used to filter only applicable notebooks.

## Apply Existing Notebook into Table

Open required table and in "Notebooks Browser" search for notebooks for applicable table. 
To apply: select Notebook and using it's context menu apply notebook into required table.
  
## Supported Languages

* Python (3.7)

Optionally: R, Julia and JavaScript.

## Environments

Each script can be run in isolated environment, with predefined packages configuration, same as for 
[Script Environments](../compute/scripting.md#Environments). Environment an be specified in Notebook properties. 

# Importing Notebook

To import Notebook just drag-and-drop the corresponding "ipynb" file to the Platform or use 
[import](../access/importing-data.md).

See also:
* [http://jupyter.org](http://jupyter.org)
* [A gallery of interesting Jupyter notebooks](https://github.com/jupyter/jupyter/wiki/A-gallery-of-interesting-Jupyter-Notebooks#statistics-machine-learning-and-data-science)
