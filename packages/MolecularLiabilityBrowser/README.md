# Molecular Liability Browser

MolecularLiabilityBrowser is a [package](https://datagrok.ai/help/develop/develop#packages) for
the [Datagrok](https://datagrok.ai) platform.

### About connections and .gitignore file

We can NOT put /connections/* files to .gitignore, particularly if the connection is referenced in any
of `queries/*.sql` files, you will get an error
`Variable "<con_name>" not found.` trying to deploy your package. But there may be a way to prevent security sensitive
strings in `connections/*.json`. Leave only meta attibutes in `connections/<con_name>.json` file:

```
{
    "name": "MLB",
    "#type": "DataConnection",
    "dataSource": "PostgreSQL",
    "description": "Molecular Liability Browser data",
    "tags": ["demo"]
}
```

and remove `parameters` and `credentials`. Deploy your package and connection will be created in Datagrok, but without
values of `Server`, `Db`, `Login` and `Password`
properties. You will be able to specify values for them. After that your aplication running can use this connection and
perform queries through it, saddly only until first redeployment when all props will be cleared. So
commit `connections/<con_name>.json` to repository without values and set values in your local files only and prevent
commit carefully.
