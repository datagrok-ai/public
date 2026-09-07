import {signal, computed, Scope, bindText} from '../../src/index.js';
import {Section} from '../../src/components/containers/section.js';

function injectOnce(id, href) {
  if (document.getElementById(id)) return;
  const l = document.createElement('link');
  l.rel = 'stylesheet';
  l.href = new URL(href, import.meta.url).href;
  l.id = id;
  document.head.append(l);
}

injectOnce('u2-icons-css', '../../css/icons.css');
injectOnce('u2-section-css', '../../css/section.css');

function el(tag, cls, text) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (text !== undefined) e.textContent = text;
  return e;
}

function button(text, onClick) {
  const b = el('button', null, text);
  b.addEventListener('click', onClick);
  return b;
}

function paragraph(text) {
  const p = el('div');
  p.style.padding = '2px 0';
  p.textContent = text;
  return p;
}

export async function render(main) {
  main.append(el('h1', null, 'Section'));
  const intro = el('p');
  intro.innerHTML = 'One standalone header + collapsible body — the form-category look, unlike ' +
    '<code>Accordion</code>\'s managed pane set. Hover the header: the chevron appears on the ' +
    '<b>left</b> of the caption, riding the host\'s padding gutter so the title stays flush with ' +
    'the content. Click (or Enter/Space) toggles; collapsed content is hidden, never detached. ' +
    '<code>FuncCallForm</code> renders every parameter category as one of these.';
  main.append(intro);

  const scopeCount = el('span', null, String(Scope.liveCount));
  const countLine = el('p');
  countLine.append('Live scopes: ', scopeCount);
  const refresh = () => scopeCount.textContent = String(Scope.liveCount);
  main.append(countLine);

  // padding gives the hover chevron its gutter — a zero-padding overflow:hidden host clips it
  const host = el('div');
  host.style.padding = '0 var(--dg-space-xl)';
  host.style.maxWidth = '520px';
  main.append(host);

  const study = new Section({title: 'Study'});
  study.add(
    paragraph('Built once and kept: collapse hides the body, it never detaches.'),
    paragraph('Type below, collapse, expand — the value survives.'));
  const memo = el('input');
  memo.placeholder = 'Scratch…';
  study.add(memo);
  host.append(study.root);

  const collapsed = new Section({title: 'Starts collapsed', expanded: false});
  collapsed.add(paragraph('expanded: false at construction.'));
  host.append(collapsed.root);

  const plain = new Section({title: 'Not collapsible', collapsible: false});
  plain.add(paragraph('collapsible: false — a plain heading and body, no chevron, no click.'));
  host.append(plain.root);

  const state = el('code');
  bindText(study.scope, state, computed(() =>
    `Study: ${study.expanded.value}   Starts collapsed: ${collapsed.expanded.value}`));
  const readout = el('p', 'u2-gallery-status');
  readout.append(state);
  main.append(readout);

  main.append(el('h2', null, 'Programmatic expand / collapse'));
  main.append(el('p', 'u2-gallery-status',
    'The signal is two-way: writing it moves the chevron and the body, exactly as a click does.'));
  const row = el('div');
  row.append(
    button('Toggle "Study"', () => study.expanded.value = !study.expanded.value), ' ',
    button('Toggle "Starts collapsed"', () => collapsed.expanded.value = !collapsed.expanded.value));
  main.append(row);

  main.append(el('h2', null, 'Adopted signal'));
  const outside = signal(true);
  const adopted = new Section({title: 'Driven from outside', expanded: outside});
  adopted.add(paragraph('The Signal handed as expanded is adopted, not copied — the owner keeps driving it.'));
  const adoptedHost = el('div');
  adoptedHost.style.padding = '0 var(--dg-space-xl)';
  adoptedHost.style.maxWidth = '520px';
  main.append(adoptedHost);
  adoptedHost.append(adopted.root);
  main.append(button('outside.value = !outside.value', () => outside.value = !outside.value));

  main.append(el('h2', null, 'Disposal'));
  main.append(el('p', 'u2-gallery-status',
    'Dispose drops the header listeners and the display effect — live scopes drop back.'));
  main.append(button('Dispose sections', () => {
    for (const s of [study, collapsed, plain, adopted])
      s.dispose();
    refresh();
  }));
  refresh();
}
