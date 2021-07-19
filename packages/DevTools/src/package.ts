import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {
  dfExts,
  entExtract,
  templates,
  helpUrls,
  tags,
  viewerConst
} from './constants';
import './styles.css';

export const _package = new DG.Package();

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

function getViewerScript(viewer: DG.Viewer): string {
  let options = viewer.getOptions(false)['look'];
  delete options['#type'];
  let script = `grok.shell.v.addViewer(DG.VIEWER.${viewerConst[viewer.type]}, ${JSON.stringify(options, null, 2)});`;
  return `<pre><code>${script}</code></pre>`;
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

  grok.events.onContextMenu.subscribe((args) => {
    let ent = args.args.context;
    let menu = args.args.menu;
    if (ent instanceof DG.Viewer) {
      let toScriptGroup = menu.group('To Script');
      let toJsScript = toScriptGroup.find('To JavaScript');
      if (!toJsScript) toScriptGroup.item('To JavaScript', () => grok.shell.info(getViewerScript(ent)));
    }
  });
}
