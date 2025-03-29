# PowerPack

Commonly used platform enhancements

Owner: [Andrew Skalkin](https://github.com/skalkin) 
Issues: [Tracker](https://github.com/datagrok-ai/public/projects/2)

## Power widgets

A start page that contains widgets (annotated with the `dashboard` tag) that are dynamically
discovered from the packages available to the current user.

## Power search

Ability to search for anything from the start screen, with the special support for the following:

* Widgets
* Functions
* Applications
* User-defined external apps (via iframe)
* Entities (connections, queries, etc)

### Search templates

The goal is to provide an efficient and extensible way to react to user input
in the global search box. A search template consists of the following:

* regular expressions (parameters would become properties)
* widget identifier
* property mapping

A collection of search templates is curated by the community (either global or
within your organization). You can mix and match template collections by editing
the `searchTemplatePaths` package property.

The following example illustrates the template that binds a ChEMBL molecule identifier
(such as CHEMBL123456) with a web page containing the molecule summary. Here is the
corresponding entry in the chembl.json file:

```
     {
      "id": "chembl-id-representations",
      "name": "Representations",
      "widget": "webWidget",
      "templates": [{
        "template": "(CHEMBL[0-9]+)",
        "url": "https://www.ebi.ac.uk/chembl/embed/#compound_report_card/${1}/representations"
      }]
    }
```

Now, when you enter `CHEMBL123456` in the search box, the following page appears:
`https://www.ebi.ac.uk/chembl/embed/#compound_report_card/CHEMBL123456/representations`

## Power Expressions

`Add New Column` dialog that lets you add derived columns by applying a formula to existing columns.
To activate, click on the `Add New Column` icon on the toolbar. Here are some of the features:

* Support for functions implemented with Python, R, Julia, JavaScript, C++, and others
* Interactive preview of results as you type
* Functions registry
* Auto-suggested functions based on the column type and semantic type
* History, saving and reusing formulas

# Windows Manager

Use icons on the right of the status bar to control visibility of tool windows.