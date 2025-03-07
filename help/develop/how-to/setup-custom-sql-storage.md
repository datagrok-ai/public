---
title: "Setup Custom SQL Storage"
---

Datagrok allows you to create custom data storage for your package. The data is stored server-side in the PostgreSQL database. This makes creating powerful CRUD applications possible directly within Datagrok. The data is stored persistently and can be accessed throughout the platform with [Datagrok Queries](./access-data.md).

## Storage Creation

Custom databases can be defined at the package level. To create the database, add a corresponding directory to the package: `{Package root}/databases/{Name}`. SQL files under this directory will be used to configure the database.

## Querying Custom SQL Storage

The data stored in the created storage can be accessed through package queries. For the connection name, use the database name, which is the same as the name of the created directory.

## Updating the Database

Databases for packages support forward migrations. This means you can add new columns but cannot make updates that are incompatible with the previous configuration.

Migrations are applied in lexicographical filename order, so it is advised to use incrementing filenames, such as `0000_init.sql`, `0001_add_optional_column.sql`.

Each new migration requires publishing the package in `--release` mode to ensure all users receive the updates.

:::warning

Please note that no rollback functionality is provided, so we advise performing each migration carefully.

:::

## Example: Compound Registration App

In this section, we provide an example of setting up the database for the compound registration app. We assume the package name is `CompoundRegistrator`.

1. Create a new folder `CompoundRegistrator/databases/compounds`.
2. In the folder, create a file `0000_init.sql` with the following content:

```sql
create table compounds.list (
    id int primary key,
    smiles varchar(512)
);

GRANT ALL ON TABLE compounds.list to :LOGIN;
```

3. Note that `compounds` is the name of the PostgreSQL database, and login configuration might be needed.
4. Under `CompoundRegistrator/queries`, create queries to register compounds and fetch the list of compounds:

```sql
--name: insertElement
--connection: Compounds
--input: string smiles
insert into compounds.list (id, smiles) values 
(floor(random() * 1000 + 1)::int, @smiles);
--end

--name: getElements
--connection: Compounds
SELECT id, smiles from compounds.list;
--end
```

5. Access the queries through the JS-API in your application! For example, you can create a form that registers a molecule:

```js
const smilesInput = ui.input.molecule('Molecule');
const addButton = ui.button('Add Compound', async () => {
  await grok.data.query('CompoundRegistrator:insertElement', {'smiles': smilesInput.value!});
  grok.shell.info('Compound added successfully!');
});
```