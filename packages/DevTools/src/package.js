import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import './styles.css';

export const _package = new DG.Package();
const dfExts = ['csv', 'xlsx'];

const templates = {
  FileInfo: (ent) =>
`(async () => {
  // Read as text
  const str = await grok.dapi.files.readAsText("${ent.fullPath}");

  // Read as dataframe
  const df = await grok.data.files.openTable("${ent.fullPath}");

  // Delete
  await grok.dapi.files.delete("${ent.fullPath}");
  grok.shell.info("${ent.fileName} exists: " + (await grok.dapi.files.exists("${ent.fullPath}")));
})();`,
  DataQuery: (ent) =>
`(async () => {
  const q = await grok.dapi.queries.find("${ent.id}");

  // Run a query
  const res = await grok.data.query("${ent.nqName}", {});
})();`,
  User: (ent) =>
`(async () => {
const user = await grok.dapi.users.include("group.memberships, group.adminMemberships").find("${ent.id}");
const userGroup = user.group;

// Find the groups the user belongs to
grok.shell.info(\`Memberships of \${user.friendlyName}: [\${userGroup.memberships}]\`);
grok.shell.info(\`Admin memberships of \${user.friendlyName}: [\${userGroup.adminMemberships}]\`);

// Manage permissions
const entity = await grok.dapi.layouts.first();
let canEdit = false;

await grok.dapi.permissions.grant(entity, userGroup, canEdit);
console.log(await grok.dapi.permissions.get(entity));
await grok.dapi.permissions.revoke(userGroup, entity);
})();`,
  Group: (ent) =>
`(async () => {
  const group = await grok.dapi.groups.find("${ent.id}");

  // Manage permissions
  const entity = await grok.dapi.layouts.first();
  let canEdit = false;

  await grok.dapi.permissions.grant(entity, group, canEdit);
  console.log(await grok.dapi.permissions.get(entity));
  await grok.dapi.permissions.revoke(group, entity);
})();`,
  DataFrame: (ent) =>
`
DataFrame snippet
`,
  Column: (ent) =>
`
Column snippet
`,
  Package: (ent) =>
`(async () => {
  // Find package functions
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
  }
})();`,
  Project: (ent) =>
`(async () => {
  const project = await grok.dapi.projects.find("${ent.id}");
  console.log(await grok.dapi.permissions.get(project));
})();`,
  Script: (ent) =>
`(async () => {
  const script = await grok.dapi.scripts.find("${ent.id}");

  // Call a script
  const res = await grok.functions.call("${ent.nqName}", {});

  // Read a script
  grok.shell.info(script.script);
})();`,
  ViewLayout: (ent) =>
`(async () => {
  // Apply to the original table
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
  console.log(savedLayout.toJson());
})();
`};

const entExtract = {
  FileInfo: (ent) => dfExts.includes(ent.extension) ? `let df = await grok.data.files.openTable("${ent.fullPath}");` : ``,
  DataQuery: (ent) => `let query = await grok.dapi.queries.find("${ent.id}");`,
  User: (ent) => `let user = await grok.dapi.users.find("${ent.id}");`,
  Group: (ent) => `let group = await grok.dapi.groups.find("${ent.id}");`,
  DataFrame: (ent) => `let df = await grok.data.openTable("${ent.id}");`,
  Column: (ent) => `let df = await grok.data.openTable("${ent.dataFrame.id}");\nlet col = df.col("${ent.name}");`,
  Package: (ent) => `let package = await grok.dapi.packages.find("${ent.id}");`,
  Project: (ent) => `let project = await grok.dapi.projects.find("${ent.id}");`,
  Script: (ent) => `let script = await grok.dapi.scripts.find("${ent.id}");`,
  Func: (ent) => `let func = DG.Func.find({ name: "${ent.name}" })[0];`,
  ViewLayout: (ent) => `let layout = await grok.dapi.layouts.find("${ent.id}");`,
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

function getTagEditor(type) {
  let t = DG.TagEditor.create();
  for (let tag of tags[type]) {
    t.addTag(tag);
  }
  return t.root;
}

function format(s) {
  s = s.replaceAll('-', ' ');
  return s[0].toUpperCase() + s.slice(1);
}

async function loadSnippets(ent) {
  const type = ent.constructor.name;
  let tags = `#demo and #${type}`;
  if (type === 'FileInfo' && dfExts.includes(ent.extension)) {
    tags += 'and #dataframe';
  }
  const snippets = (await grok.dapi.scripts.list({ filter: tags }));
  return snippets.slice(0, 3);
}

//tags: autostart
export function describeCurrentObj() {
  grok.events.onAccordionConstructed.subscribe(async (acc) => {
    const ent = acc.context;
    const type = ent.constructor.name;

    if (ent) {
      const snippets = await loadSnippets(ent);
      const template = (type in templates) ? templates[type](ent) : '';
      if (snippets.length === 0 && !template) return;

      let links = helpUrls[type] || [];
      links = Object.keys(links).map(key => ui.link(`${type} ${key}`, links[key]));

      const snippetNames = snippets.map(s => ui.divText(format(s.friendlyName), { classes: 'd4-link-action' }));
      let editor = ui.textInput('', template);
      editor.input.style = 'height: 200px; overflow: hidden;';
      // editor.input.style = 'width: 0; height: 0; visibility: hidden;';
      // editor.root.style.display = 'none';

      snippetNames.forEach((el, idx) => el.addEventListener('click', () => {
        editor.value = snippets[idx].script;
        // editor.root.style.display = 'block';
        // editor.input.style = 'width: 200; height: 300; visibility: visible;';
      }));

      const clipboardBtn = ui.button(ui.iconFA('copy'), () => {
        editor.input.select();
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

      const editorBtn = ui.button(ui.iconFA('external-link'), () => {
        grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, { script: editor.value }));
      }, 'Open in editor');
      $(editorBtn).addClass('dt-snippet-editor-icon dt-editor-icon');

      const resetBtn = ui.button(ui.iconFA('redo'), () => editor.value = template, 'Reset');
      $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

      const topEditorBtn = ui.button(ui.iconFA('edit'), () => {
        grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, { script: entExtract[type](ent) }));
      }, 'Open in editor');
      $(topEditorBtn).addClass('dt-snippet-inline-icon');

      const browserLogBtn = ui.button(ui.iconFA('terminal'), () => (console.clear(), console.log(grok.shell.o)), 'Log to console');
      $(browserLogBtn).addClass('dt-snippet-inline-icon');

      let snippetsPane = acc.getPane('Snippets');      
      if (!snippetsPane) snippetsPane = acc.addPane('Snippets', () => {
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
