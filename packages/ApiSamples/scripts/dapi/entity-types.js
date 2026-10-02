//name: entity-types
//language: javascript

// Entity types are the named kinds a sticky meta schema attaches to (`/api/entities/types`).
// grok.dapi.entityTypes is an ordinary data source, so it takes the usual query calls.

const types = await grok.dapi.entityTypes.by(100).order('name').list();

const df = DG.DataFrame.create(types.length);
df.name = 'Entity types';
df.columns.addNewString('name').init((i) => types[i].name);
df.columns.addNewString('matching').init((i) => types[i].matching);

grok.shell.addTableView(df);
