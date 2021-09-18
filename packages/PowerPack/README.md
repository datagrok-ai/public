# PowerPack

Commonly used platform enhancements

Owner: [Andrew Skalkin](https://github.com/skalkin)  
Issues: https://github.com/datagrok-ai/public/projects/2

## Power widgets

A start page that contains widgets (annotated with the `dashboard` tag) that are dynamically discovered from the packages available to the current user.

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

Now, when a user enters `CHEMBL123456` in the search box, he sees the content of the following page:
`https://www.ebi.ac.uk/chembl/embed/#compound_report_card/CHEMBL123456/representations`