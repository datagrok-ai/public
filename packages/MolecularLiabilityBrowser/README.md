# Molecular Liability Browser

MolecularLiabilityBrowser is a [package](https://datagrok.ai/help/develop/develop#packages) for
the [Datagrok](https://datagrok.ai) platform. MLB application allows to select variable regions 
linked with antigen and perform Composition Analysis ([WebLogo](https://datagrok.ai/help/visualize/viewers/web-logo))
for complementarity-determining regions (CDR) defined as 'chothia', 'aroop', 'north', 'martin', 'kabat' 
on antibody's residue position numbering schemes 'imgt', 'aho', 'chothia', 'kabat' (requires preliminary 
calculated data to be loaded into database) separately for 'Light' and 'Heavy' chains. Allows filtering variable 
regions for predicted and observed post-translational modifications (PTM) of residues. Exploring the ngl/3D view for PDB 
structures allows highlighting CDR or Paratope regions/areas or particular residues annotated for PTM.  

### About connections and .gitignore file

We can NOT put /connections/* files to .gitignore, particularly if the connection is referenced in any
of `queries/*.sql` files, you will get an error
`Variable "<con_name>" not found.` trying to deploy your package. But there may be a way to prevent security sensitive
strings in `connections/*.json`. Leave only meta attibutes in `connections/<con_name>.json` file:

```
{
    "name": "MLB",
    "#type": "DataConnection",
    "dataSource": "PostgresDart",
    "description": "Molecular Liability Browser data",
    "tags": ["demo"]
}
```

and remove `parameters` and `credentials`. Deploy your package and connection will be created in Datagrok, but without
values of `Server`, `Db`, `Login` and `Password`
properties. You will be able to specify values for them. After that your application running can use this connection and
perform queries through it, sadly only until first redeployment when all props will be cleared. So
commit `connections/<con_name>.json` to repository without values and set values in your local files only and prevent
commit carefully.
