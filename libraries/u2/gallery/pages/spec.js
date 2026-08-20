import {signal, Scope, Control, SpecContext, parseSpec, renderSpec} from '../../src/index.js';
import {divH, span, button} from '../../src/core/elements.js';
import {registerAll} from '../../src/spec/registrations.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.id = id;
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  document.head.append(l);
}

injectOnce('u2-elements-css', '../../css/elements.css');
injectOnce('u2-inputs-css', '../../css/inputs.css');
injectOnce('u2-choice-css', '../../css/choice.css');
injectOnce('u2-number-css', '../../css/number.css');
injectOnce('u2-form-css', '../../css/form.css');
injectOnce('u2-spec-css', '../../css/spec.css');

registerAll();

const DEMO = `{
  "$schema": "dg-ui/1",
  "root": {
    "tag": "u2-div-v",
    "children": [
      {"tag": "h2", "props": {"text": "Compound"}},
      {"tag": "u2-form", "children": [
        {"tag": "u2-text-input", "props": {"label": "Name", "name": "name"},
         "bind": {"value": "$.name"}},
        {"tag": "u2-number-input", "props": {"label": "Amount", "name": "amount", "min": 0, "step": 5},
         "bind": {"value": "$.amount"}},
        {"tag": "u2-text-input", "props": {"label": "Comment", "value": "not bound, so dump folds it in"}}
      ]},
      {"tag": "u2-button", "props": {"text": "Save", "primary": true}, "on": {"click": "cmd:save"}}
    ]
  }
}`;

const CORPUS = `{
  "$schema": "dg-ui/1",
  "root": {"tag": "u2-form", "children": [
    {"tag": "u2-text-input", "props": {"label": "Name", "name": "name"}, "bind": {"value": "$.name"}},
    {"tag": "u2-number-input", "props": {"label": "Age", "name": "age", "mode": "int", "min": 0},
     "bind": {"value": "$.age"}},
    {"tag": "u2-button", "props": {"text": "Save", "primary": true}, "on": {"click": "cmd:save"}}
  ]}
}`;

const COMMENT = 'Comment';

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function readout(label, source) {
  return divH([span(`${label} = `), span(source)], 'u2-gallery-status');
}

function editor(text, rows) {
  const area = el('textarea');
  area.value = text;
  area.rows = rows;
  area.spellcheck = false;
  area.style.width = '100%';
  area.style.boxSizing = 'border-box';
  area.style.fontFamily = 'var(--dg-font-family-mono)';
  area.style.fontSize = 'var(--dg-font-size-small)';
  return area;
}

/** The live value of the first node carrying `props.label`, as it came back out of dump(). */
function dumpedValue(node, label) {
  if (node.props && node.props.label === label)
    return node.props.value;
  for (const child of node.children ?? []) {
    const found = dumpedValue(child, label);
    if (found !== undefined)
      return found;
  }
  return undefined;
}

export async function render(main) {
  const page = new Control();
  main.append(el('h1', null, 'Spec (dg-ui/1)'));
  const intro = el('p');
  intro.innerHTML = 'A spec is JSON: registered <code>u2-*</code> tags, validated props, ' +
    '<code>bind</code> entries pointing into the render context, and <code>on</code> entries naming ' +
    'context commands — never code. Render it, edit the live UI, <code>dump()</code> it back to JSON ' +
    'with the current values folded in, and render that again: the round-trip is what makes an ' +
    'LLM-emitted layout persistable and the dashboard editor possible.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  const name = signal('Aspirin');
  const amount = signal(10);
  const saved = signal('nothing saved yet');
  const ctx = new SpecContext({
    data: {name, amount},
    commands: {save: () => saved.value = `name = ${JSON.stringify(name.peek())}, amount = ${amount.peek()}`},
  });

  main.append(el('h2', null, 'Spec'));
  const source = editor(DEMO, 20);
  main.append(source);

  const status = el('p', 'u2-gallery-status', 'Not rendered yet.');
  const host = el('div');
  host.style.border = 'var(--dg-border-width) solid var(--dg-border)';
  host.style.borderRadius = 'var(--dg-radius)';
  host.style.padding = 'var(--dg-space-m)';
  host.style.minHeight = '40px';
  const dump = editor('', 14);
  dump.readOnly = true;

  let instance = null;
  const show = (text) => {
    if (instance)
      instance.dispose();
    instance = null;
    host.textContent = '';
    status.style.color = '';
    try {
      instance = renderSpec(parseSpec(text), ctx);
      host.append(instance.root);
      status.textContent = 'Rendered. Nodes that could not be built show as red placeholders in the frame.';
    } catch (e) {
      status.textContent = `parseSpec: ${e.message}`;
      status.style.color = 'var(--dg-failure)';
    }
    refresh();
  };
  page.own(() => {
    if (instance)
      instance.dispose();
  });

  const buttons = divH([
    button('Render', () => show(source.value), {primary: true}),
    button('Dump', () => {
      if (!instance) {
        status.textContent = 'Nothing to dump — press Render first.';
        return;
      }
      const spec = instance.dump();
      dump.value = JSON.stringify(spec, null, 2);
      const folded = dumpedValue(spec.root, COMMENT);
      status.style.color = '';
      status.textContent = folded === undefined ?
        'Dumped.' :
        `Dumped — the unbound "${COMMENT}" input folded its live value in as ${JSON.stringify(folded)}, ` +
        'while the bound props stayed bindings.';
    }),
    button('Restore', () => show(dump.value)),
  ]);
  buttons.style.gap = 'var(--dg-space-m)';
  buttons.style.margin = 'var(--dg-space-m) 0';
  main.append(buttons, status, host);

  main.append(el('h2', null, 'Render context'));
  const context = Control.build(() => [
    readout('$.name', name),
    readout('$.amount', amount),
    readout('cmd:save', saved),
  ]);
  main.append(context.root);
  main.append(el('p', 'u2-gallery-status',
    'Both inputs took the context signals themselves, so typing moves these readouts and writing the ' +
    'signals from code moves the editors. The Save button carries no code: it names the command above.'));

  main.append(el('h2', null, 'Dump'));
  main.append(dump);
  main.append(el('p', 'u2-gallery-status',
    'Type into "Comment", press Dump, then press Restore: the restored tree comes up with the value ' +
    'you typed, and dumping it again yields the same JSON.'));

  main.append(el('h2', null, 'Manifest'));
  const manifest = el('div');
  main.append(el('p', 'u2-gallery-status',
    'registry.manifest() — checked in as manifest.json by "npm run manifest", and the only thing an ' +
    'LLM needs to read before emitting a spec.'), manifest);

  main.append(el('h2', null, 'LLM corpus'));
  const details = el('details');
  details.append(el('summary', null, 'Prompt → spec: the few-shot pair the manifest exists for'));
  details.append(el('p', 'u2-gallery-status',
    'Prompt: "Build a save-form with name and age, bound to the data, with a Save button."'));
  const corpus = el('pre', null, CORPUS);
  corpus.style.fontFamily = 'var(--dg-font-family-mono)';
  corpus.style.fontSize = 'var(--dg-font-size-small)';
  details.append(corpus);
  main.append(details);

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'A SpecInstance is one scope owning every component the spec built, including the bridges into the ' +
    'context signals. Disposing it — on re-render, on this button, or when the page goes away — takes ' +
    'all of them with it, and live scopes drop back to the page baseline.'));
  main.append(button('Dispose the rendered spec', () => {
    if (instance)
      instance.dispose();
    instance = null;
    host.textContent = '';
    status.textContent = 'Disposed.';
    refresh();
  }));

  show(source.value);

  try {
    const response = await fetch(new URL('../../manifest.json', import.meta.url));
    if (!response.ok)
      throw new Error(`${response.status} ${response.statusText}`);
    const json = await response.json();
    manifest.append(el('p', 'u2-gallery-status', `${json.schema} · ${json.components.length} components`));
    for (const meta of json.components) {
      const row = el('p', 'u2-gallery-status');
      row.append(el('code', null, meta.tag), ` — ${meta.description} · ${meta.props.length} props`);
      manifest.append(row);
    }
  } catch (e) {
    manifest.append(el('p', 'u2-gallery-status',
      `manifest.json could not be read (${e.message}) — run "npm run manifest".`));
  }
}
