import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';
import CodeMirror from 'codemirror';
import {
  dfExts,
  entExtract,
  templates,
  helpUrls,
  tags,
  viewerConst,
  EntityType,
  supportedEntityTypes,
} from './constants';


function getGroupInput(codeMirror: CodeMirror.EditorFromTextArea, type: string): HTMLElement {
  const items = tags[type];
  const inp = ui.input.choice('See snippets:', {value: items.length ? items[0] : null, items: items,
    onValueChanged: async (value) => {
      const snippets = await loadSnippets(type, value);
      const container = $('.dt-dev-pane-container > .dt-snippet-section');
      container.empty();
      container.append(formSnippetSection(codeMirror, snippets));
    }});
  return inp.root;
}

async function loadSnippets(type: string, tag: string | null = null): Promise<DG.Script[]> {
  let tags = `#demo and #${type}`;
  if (tag) tags += `and #${tag}`;
  const snippets = (await grok.dapi.scripts.list({filter: tags}));
  return snippets;
}

function format(s: string): string {
  s = s.replace(/-/g, ' ');
  return s[0].toUpperCase() + s.slice(1);
}

function formSnippetSection(codeMirror: CodeMirror.EditorFromTextArea,
  snippets: DG.Script[], count: number = 3): HTMLDivElement[] {
  const snippetNames = snippets.map((s) => ui.divText(format(s.friendlyName), {classes: 'd4-link-action'}));
  const lineEndRegex = /\r\n|\r|\n/;
  const paramAnnotationRegex = /^\/\/(name|tags|language|help-url):.+/;
  snippetNames.forEach((el, idx) => el.addEventListener('click', () => {
    let s = '';
    const lines = snippets[idx].script.split(lineEndRegex);
    for (const line of lines) {
      if (paramAnnotationRegex.test(line))
        continue;
      s += line + '\n';
    }
    codeMirror.setValue(s);
  }));
  if (snippetNames.length > count) {
    const rest = snippetNames.splice(count);
    const ellipsis = ui.iconFA('ellipsis-h', () => {
      $(ellipsis.parentElement.parentElement).append(rest);
      ellipsis.remove();
    }, 'Show more examples');
    snippetNames.push(ui.div(ellipsis));
  }
  return snippetNames;
}

function getViewerScript(viewer: DG.Viewer): string {
  const options = viewer.getOptions(false)['look'];
  delete options['#type'];
  const type = viewerConst[viewer.type];
  const script = `grok.shell.v.addViewer(${type ? `DG.VIEWER.${type}` :
    `'${viewer.type}'`}, ${JSON.stringify(options, null, 2)});`;
  return `<pre><code>${script}</code></pre>`;
}

export function addToJSContextCommand({args: {menu, context}}) {
  const toScriptGroup = menu.group('To Script');
  const toJsScript = toScriptGroup.find('To JavaScript');
  if (!toJsScript) toScriptGroup.item('To JavaScript', () => grok.shell.info(getViewerScript(<DG.Viewer>context)));
}

export async function _renderDevPanel(ent: EntityType, minifiedClassNameMap: {}): Promise<DG.Widget> {
  if (ent == null)
    return DG.Widget.fromRoot(ui.divText('Entity does not exist.', {style: {color: 'var(--failure)'}}));


  let type = ent.constructor.name;
  if (!supportedEntityTypes.includes(type) && type in minifiedClassNameMap)
    type = minifiedClassNameMap[type].find((c) => ent instanceof eval(`DG.${c}`)) ?? type;

  const snippets = await loadSnippets(type,
    (ent instanceof DG.FileInfo && dfExts.includes(ent.extension)) ? 'dataframe' :
      (ent instanceof DG.DataFrame || ent instanceof DG.Column) ? tags[type][0] :
        null);
  const template = (type in templates) ? templates[type](ent) : '';

  if (snippets.length === 0 && !template)
    return DG.Widget.fromRoot(ui.divText(`Unsupported entity: ${type}.`, {style: {color: 'var(--failure)'}}));

  let links = helpUrls[type] || [];
  links = Object.keys(links).map((key) => key === 'additional' ?
    Object.keys(links[key]).map((title) => ui.link(title, links[key][title], 'Open wiki reference')) :
    ui.link(`${type} ${key}`, links[key], `Open ${key} reference`));

  const editor = ui.input.textArea('', {value: template});
  const codeMirror = CodeMirror.fromTextArea(editor.input as HTMLTextAreaElement, {
    mode: 'javascript',
    readOnly: false,
    lineNumbers: false,
    showCursorWhenSelecting: false,
    lineWrapping: true,
  });
  setTimeout(() => codeMirror.refresh(), 0);
  codeMirror.getWrapperElement().classList.add('dt-cd-wrapper');

  const playBtn = ui.button(ui.iconFA('play'), () => {
    eval(`(async () => {\n${codeMirror.getValue()}\n})()`); // TODO: script approval
  }, 'Run');
  $(playBtn).addClass('dt-snippet-editor-icon dt-play-icon');

  const clipboardBtn = ui.button(ui.iconFA('copy'), () => {
    navigator.clipboard.writeText(codeMirror.getValue());
    const copyIcon = clipboardBtn.removeChild(clipboardBtn.firstChild);
    clipboardBtn.appendChild(ui.iconFA('clipboard-check'));
    setTimeout(() => {
      clipboardBtn.removeChild(clipboardBtn.firstChild);
      clipboardBtn.appendChild(copyIcon);
    }, 1000);
  }, 'Copy');
  $(clipboardBtn).addClass('dt-snippet-editor-icon dt-clipboard-icon');

  const editorBtn = ui.button(ui.iconFA('external-link-square'), () => {
    grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, {script: codeMirror.getValue()}));
  }, 'Open in editor');
  $(editorBtn).addClass('dt-snippet-editor-icon dt-editor-icon');

  const resetBtn = ui.button(ui.iconFA('redo'), () => codeMirror.setValue(template), 'Reset');
  $(resetBtn).addClass('dt-snippet-editor-icon dt-reset-icon');

  const topEditorBtn = ui.button(ui.iconFA('edit'), () => {
    grok.shell.addView(DG.View.createByType(DG.View.JS_EDITOR, {script: entExtract[type](ent)}));
  }, 'Open in editor');
  $(topEditorBtn).addClass('dt-snippet-inline-icon');

  const browserLogBtn = ui.button(ui.iconFA('terminal'), () => {
    console.log(grok.shell.o);
    grok.shell.info('The object was printed to console. Press F12 to open the developer tools.');
  }, 'Log to console');
  $(browserLogBtn).addClass('dt-snippet-inline-icon');

  return DG.Widget.fromRoot(ui.divV([
    ui.divH([ui.divText(`${type} ${ent.name}:`), topEditorBtn, browserLogBtn], {style: {alignItems: 'baseline'}}),
    ...links.map((link: HTMLAnchorElement | HTMLAnchorElement[]) =>
      Array.isArray(link) ? ui.div([ui.divText('See also:'), ...link]) : link),
    ...((type in tags) ? [getGroupInput(codeMirror, type)] : []),
    ui.div(formSnippetSection(codeMirror, snippets), 'dt-snippet-section'),
    ui.divV([playBtn, clipboardBtn, editorBtn, resetBtn, editor.root], 'dt-textarea-box'),
  ], 'dt-dev-pane-container'));
}

export function getMinifiedClassNameMap(): { [index: string]: string[] } {
  return supportedEntityTypes.reduce((map, t) => {
    const minClassName = eval(`DG.${t}.name`);
    map[minClassName] = [...(map[minClassName] || []), t];
    return map;
  }, {});
}

export function hasSupportedType(ent: EntityType, minifiedClassNameMap: {}): boolean {
  let type = ent.constructor.name;
  if (!supportedEntityTypes.includes(type) && type in minifiedClassNameMap)
    type = minifiedClassNameMap[type].find((c) => ent instanceof eval(`DG.${c}`)) ?? type;
  return supportedEntityTypes.includes(type);
}
