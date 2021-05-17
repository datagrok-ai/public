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
})();`,
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
      const snippetNames = snippets.map(s => ui.divText(format(s.friendlyName), { classes: 'd4-link-action' }));
      let editor = ui.textInput('', (type in templates) ? templates[type](ent) : '');
      editor.input.style.height = '200px';
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
          ...snippetNames,
          ui.divV([clipboardBtn, editor.root], 'textarea-box')
        ]);
      });
    }
  });
}
