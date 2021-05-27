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
  DataQuery: (ent) => `grok.data.query("${ent.nqName}", {}).then(t => grok.shell.info(t.rowCount));`,
  User: (ent) =>
`(async () => {
const user = await grok.dapi.users.find("${ent.id}");
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
  Package: (ent) =>
`(async () => {
  // Call a function
  const res = await grok.functions.call("${ent.name}:test", {});

  // Read credentials
  const package = await grok.dapi.packages.find("${ent.id}");
  const credentialsResponse = await package.getCredentials();
  if (credentialsResponse == null) {
    grok.shell.info("Credentials are not set.");
  } else {
    grok.shell.info(credentialsResponse.parameters);
  }
})();`,
  Script: (ent) => `grok.functions.call("${ent.nqName}", {}).then(res => grok.shell.info(res));`,
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

const helpUrls = {
  DataConnection: ['https://datagrok.ai/help/develop/how-to/access-data', 'https://datagrok.ai/js-api/DataConnection'],
  DataQuery: ['https://datagrok.ai/help/develop/how-to/access-data', 'https://datagrok.ai/js-api/DataQuery'],
  FileInfo: ['https://datagrok.ai/help/develop/how-to/access-data', 'https://datagrok.ai/js-api/FileInfo'],
  Script: ['https://datagrok.ai/help/develop/scripting', 'https://datagrok.ai/js-api/Script'],
  ViewLayout: ['https://datagrok.ai/help/develop/how-to/layouts', 'https://datagrok.ai/js-api/ViewLayout'],
};

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
      links = links.map(link => {
        let anc = document.createElement('a');
        anc.innerHTML = link;
        anc.setAttribute('href', link);
        anc.setAttribute('target', '_blank');
        return anc;
      });

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
      $(clipboardBtn).addClass('clipboard-icon');

      let snippetsPane = acc.getPane('Snippets');      
      if (!snippetsPane) snippetsPane = acc.addPane('Snippets', () => {
        return ui.divV([
          ...links,
          ...snippetNames,
          ui.divV([clipboardBtn, editor.root], 'textarea-box')
        ]);
      });
    }
  });
}
