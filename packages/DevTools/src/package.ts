import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './styles.css';

export const _package = new DG.Package();
const dfExts = ['csv', 'xlsx'];

const templates = {
  FileInfo: (ent: DG.FileInfo) =>
`// Read as text
const str = await grok.dapi.files.readAsText("${ent.fullPath}");

// Read as dataframe
const df = await grok.data.files.openTable("${ent.fullPath}");

// Delete
await grok.dapi.files.delete("${ent.fullPath}");
grok.shell.info("${ent.fileName} exists: " + (await grok.dapi.files.exists("${ent.fullPath}")));`,
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
`
DataFrame snippet
`,
  Column: (ent: DG.Column) =>
`
Column snippet
`,
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
console.log(savedLayout.toJson());`
};

const entExtract = {
  FileInfo: (ent: DG.FileInfo) => dfExts.includes(ent.extension) ? `let df = await grok.data.files.openTable("${ent.fullPath}");` : ``,
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
};

const helpUrls = {
  Column: { wiki: 'https://datagrok.ai/help/overview/table', class: 'https://datagrok.ai/js-api/Column' },
  DataConnection: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/DataConnection' },
  DataFrame: { wiki: 'https://datagrok.ai/help/overview/table', class: 'https://datagrok.ai/js-api/DataFrame' },
  DataQuery: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/DataQuery' },
  FileInfo: { wiki: 'https://datagrok.ai/help/develop/how-to/access-data', class: 'https://datagrok.ai/js-api/FileInfo' },
  Script: { wiki: 'https://datagrok.ai/help/develop/scripting', class: 'https://datagrok.ai/js-api/Script' },
  ViewLayout: { wiki: 'https://datagrok.ai/help/develop/how-to/layouts', class: 'https://datagrok.ai/js-api/ViewLayout' },
};

const tags = {
  DataFrame: ['construction', 'modification', 'events'],
  Column: ['creation'],
};

function getTagEditor(type: string): HTMLElement {
  let t = DG.TagEditor.create();
  for (let tag of tags[type]) {
    t.addTag(tag);
  }
  return t.root;
}

function format(s): string {
  s = s.replaceAll('-', ' ');
  return s[0].toUpperCase() + s.slice(1);
}

async function loadSnippets(ent: any): Promise<DG.Script[]> {
  const type = ent.constructor.name;
  let tags = `#demo and #${type}`;
  if (type === 'FileInfo' && dfExts.includes(ent.extension)) {
    tags += 'and #dataframe';
  }
  const snippets = (await grok.dapi.scripts.list({ filter: tags }));
  return snippets.slice(0, 3);
}

//tags: autostart
export function describeCurrentObj(): void {
  grok.events.onAccordionConstructed.subscribe(async (acc) => {
    const ent = acc.context;
    const type = ent.constructor.name;

    if (ent) {
      const snippets = await loadSnippets(ent);
      const template = (type in templates) ? templates[type](ent) : '';
      if (snippets.length === 0 && !template) return;

      let links = helpUrls[type] || [];
      links = Object.keys(links).map(key => ui.link(`${type} ${key}`, links[key], `Open ${key} reference`));

      const snippetNames = snippets.map(s => ui.divText(format(s.friendlyName), { classes: 'd4-link-action' }));
      let editor = ui.textInput('', template);
      (editor.input as HTMLInputElement).style.height = '200px';
      (editor.input as HTMLInputElement).style.overflow = 'hidden';

      snippetNames.forEach((el, idx) => el.addEventListener('click', () => {
        editor.value = snippets[idx].script;
      }));

      const clipboardBtn = ui.button(ui.iconFA('copy'), () => {
        (editor.input as HTMLInputElement).select();
        const copied = document.execCommand('copy');
        if (copied) {
          const copyIcon = clipboardBtn.removeChild(clipboardBtn.firstChild);
          clipboardBtn.appendChild(ui.iconFA('clipboard-check'));
          setTimeout(() => {
            clipboardBtn.removeChild(clipboardBtn.firstChild);
            clipboardBtn.appendChild(copyIcon);
          }, 1000);
        }
      }, 'Copy');
      $(clipboardBtn).addClass('dt-snippet-editor-icon dt-clipboard-icon');

      const editorBtn = ui.button(ui.iconFA('external-link-square'), () => {
        grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, { script: editor.value }));
      }, 'Open in editor');
      $(editorBtn).addClass('dt-snippet-editor-icon dt-editor-icon');

      const resetBtn = ui.button(ui.iconFA('redo'), () => editor.value = template, 'Reset');
      $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

      const topEditorBtn = ui.button(ui.iconFA('edit'), () => {
        grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, { script: entExtract[type](ent) }));
      }, 'Open in editor');
      $(topEditorBtn).addClass('dt-snippet-inline-icon');

      const browserLogBtn = ui.button(ui.iconFA('terminal'), () => {
        console.clear();
        console.log(grok.shell.o);
        grok.shell.info('The object was printed to console. Press F12 to open the developer tools.');
      }, 'Log to console');
      $(browserLogBtn).addClass('dt-snippet-inline-icon');

      let devPane = acc.getPane('Dev');      
      if (!devPane) devPane = acc.addPane('Dev', () => {
        return ui.divV([
          ui.divH([ui.divText(`${type} ${ent.name}:`), topEditorBtn, browserLogBtn], { style: { 'align-items': 'baseline' } }),
          ...((type in tags) ? [getTagEditor(type)] : []),
          ...links,
          ...snippetNames,
          ui.divV([clipboardBtn, editorBtn, resetBtn, editor.root], 'dt-textarea-box'),
        ]);
      });
    }
  });
}
