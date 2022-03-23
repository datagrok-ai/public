import * as DG from 'datagrok-api/dg';


export type EntityType = DG.Entity | DG.DataFrame | DG.Column | DG.ViewBase;

export const supportedEntityTypes = [
  'FileInfo', 'DataConnection', 'DataQuery', 'User', 'Group', 'DataFrame',
  'Column', 'Package', 'Project', 'Script', 'Func', 'ViewLayout', 'View',
];

export const dfExts = ['csv', 'xlsx'];

export const entExtract = {
  FileInfo: (ent: DG.FileInfo) => dfExts.includes(ent.extension) ? `let df = await grok.data.files.openTable("${ent.fullPath}");` : ``,
  DataConnection: (ent: DG.DataConnection) => `let connection = await grok.dapi.connections.find('${ent.id}');`,
  DataQuery: (ent: DG.DataQuery) => `let query = await grok.dapi.queries.find("${ent.id}");`,
  User: (ent: DG.User) => `let user = await grok.dapi.users.find("${ent.id}");`,
  Group: (ent: DG.Group) => `let group = await grok.dapi.groups.find("${ent.id}");`,
  DataFrame: (ent: DG.DataFrame) => `let df = await grok.data.openTable("${ent.id}");`,
  Column: (ent: DG.Column) => `let df = await grok.data.openTable("${ent.dataFrame.id}");\nlet col = df.col("${ent.name}");`,
  Package: (ent: DG.Package) => `let package = await grok.dapi.packages.find("${ent.id}");`,
  Project: (ent: DG.Project) => `let project = await grok.dapi.projects.find("${ent.id}");`,
  Script: (ent: DG.Script) => `let script = await grok.dapi.scripts.find("${ent.id}");`,
  Func: (ent: DG.Func) => `let func = DG.Func.find({ name: "${ent.name}" })[0];`,
  ViewLayout: (ent: DG.ViewLayout) => `let layout = await grok.dapi.layouts.find("${ent.id}");`,
  View: (ent: DG.View) => `let view = grok.shell.view("${ent.name}");`,
};

export const helpUrls = {
  Column: { wiki: 'https://datagrok.ai/help/overview/table', class: 'https://datagrok.ai/js-api/classes/dg.Column' },
  DataConnection: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/classes/dg.DataConnection' },
  DataFrame: { wiki: 'https://datagrok.ai/help/overview/table', class: 'https://datagrok.ai/js-api/classes/dg.DataFrame' },
  DataQuery: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/classes/dg.DataQuery' },
  FileInfo: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/classes/dg.FileInfo' },
  Script: { wiki: 'https://datagrok.ai/help/compute/scripting', class: 'https://datagrok.ai/js-api/classes/dg.Script' },
  ViewLayout: { wiki: 'https://datagrok.ai/help/develop/how-to/layouts', class: 'https://datagrok.ai/js-api/classes/dg.ViewLayout' },
  View: { wiki: 'https://datagrok.ai/help/develop/ui#views', class: 'https://datagrok.ai/js-api/classes/dg.View' },
};

export const tags = {
  DataFrame: ['construction', 'modification', 'events'],
  Column: ['construction', 'modification', 'access'],
};

export const templates = {
  FileInfo: (ent: DG.FileInfo) =>
`// Read as text
const str = await grok.dapi.files.readAsText("${ent.fullPath}");

// Read as dataframe
const df = await grok.data.files.openTable("${ent.fullPath}");`,
  DataConnection: (ent: DG.DataConnection) =>
`const connection = await grok.dapi.connections.find('${ent.id}');

// Get the connection's parameters
grok.shell.info(JSON.stringify(connection.parameters));

// Find the connection's queries
const queries = await grok.dapi.queries.filter(\`connection.id = "\${connection.id}"\`).list();
grok.shell.info('Found queries: ' + queries.length);`,
  DataQuery: (ent: DG.DataQuery) =>
`const q = await grok.dapi.queries.find("${ent.id}");

// Run a query
const res = await grok.data.query("${ent.nqName}", {});`,
  User: (ent: DG.User) =>
`const user = await grok.dapi.users.include("group.memberships, group.adminMemberships").find("${ent.id}");
const userGroup = user.group;

// Find the groups the user belongs to
grok.shell.info(\`Memberships of \${user.friendlyName}: [\${userGroup.memberships}]\`);
grok.shell.info(\`Admin memberships of \${user.friendlyName}: [\${userGroup.adminMemberships}]\`);

// Manage permissions
const entity = await grok.dapi.layouts.first();
let canEdit = false;

await grok.dapi.permissions.grant(entity, userGroup, canEdit);
console.log(await grok.dapi.permissions.get(entity));
await grok.dapi.permissions.revoke(userGroup, entity);`,
  Group: (ent: DG.Group) =>
`const group = await grok.dapi.groups.find("${ent.id}");

// Manage permissions
const entity = await grok.dapi.layouts.first();
let canEdit = false;

await grok.dapi.permissions.grant(entity, group, canEdit);
console.log(await grok.dapi.permissions.get(entity));
await grok.dapi.permissions.revoke(group, entity);`,
  DataFrame: (ent: DG.DataFrame) =>
`const df = grok.shell.table("${ent.name}");`,
  Column: (ent: DG.Column) =>
`const df = grok.shell.table("${ent.dataFrame.name}");
const col = df.col("${ent.name}");`,
  Package: (ent: DG.Package) =>
`// Find package functions
const funcs = DG.Func.find({ package: "${ent.name}" });

// Call a function
const res1 = await grok.functions.call("${ent.name}:test", {});
const res2 = await funcs[0].apply({});

// Read credentials
const package = await grok.dapi.packages.find("${ent.id}");
const credentialsResponse = await package.getCredentials();
if (credentialsResponse == null) {
  grok.shell.info("Credentials are not set.");
} else {
  grok.shell.info(credentialsResponse.parameters);
}`,
  Project: (ent: DG.Project) =>
`// Find a project by its name and open it in the workspace
const project = await grok.dapi.projects.open("${ent.friendlyName}");

// Read permissions
console.log(await grok.dapi.permissions.get(project));`,
  Script: (ent: DG.Script) =>
`const script = await grok.dapi.scripts.find("${ent.id}");

// Call a script
const res = await grok.functions.call("${ent.nqName}", {});

// Read a script
grok.shell.info(script.script);`,
  Func: (ent: DG.Func) =>
`// Find a function by package, name, tags, or returnType
const f = DG.Func.find({ name: "${ent.name}" })[0];

// Work with parameters
const paramDefaults = { num: 1, string: 'a', bool: true };
const params = f.inputs.reduce((obj, p) => {
  obj[p.name] = p.defaultValue ?? paramDefaults[p.propertyType];
  return obj;
}, {});

// Call a function
const res = await f.apply(params);`,
  ViewLayout: (ent: DG.ViewLayout) =>
`// Apply to the original table
const layout = await grok.dapi.layouts.find("${ent.id}");
const tableId = JSON.parse(layout.viewState).tableId;
const df = await grok.data.openTable(tableId);
const view = grok.shell.addTableView(df);
view.loadLayout(layout);

// Get layouts applicable to the dataframe 
const layouts = await grok.dapi.layouts.getApplicable(df);
grok.shell.info("Layouts found: " + layouts.length);

// Save and serialize
const savedLayout = view.saveLayout();
console.log(savedLayout.toJson());`,
  View: (ent: DG.View) =>
`// Get a view by its name
const view = grok.shell.view("${ent.name}");

// Set the current view
grok.shell.v = view

// Append new elements
view.append(ui.h1('New fascinating content'))`,
};

export const viewerConst = Object.entries(DG.VIEWER)
  .reduce((obj, entry) => {
    const [key, value] = entry;
    obj[value] = key;
    return obj;
  }, {});
