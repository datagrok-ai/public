import {Scope} from '../../src/index.js';
import {VirtualTree} from '../../src/components/collections/tree.js';

const FOLDERS = 100;
const SUBFOLDERS = 10;
const LEAVES = 100;
const LOAD_DELAY = 150;

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

function ensureCss(name) {
  const id = `u2-${name}-css`;
  if (document.getElementById(id)) return;
  const link = el('link');
  link.id = id;
  link.rel = 'stylesheet';
  link.href = new URL(`../../css/${name}.css`, import.meta.url).href;
  document.head.append(link);
}

function later(value, ms) {
  return new Promise((resolve) => setTimeout(() => resolve(value), ms));
}

function makeRoots() {
  const roots = [];
  for (let f = 0; f < FOLDERS; f++) {
    const subfolders = [];
    for (let s = 0; s < SUBFOLDERS; s++) {
      const id = `f${f}/s${s}`;
      subfolders.push({
        id,
        label: `Subfolder ${s}`,
        children: () => later(Array.from({length: LEAVES},
          (_, i) => ({id: `${id}/l${i}`, label: `Item ${i}`, data: {folder: f, sub: s, leaf: i}})), LOAD_DELAY),
      });
    }
    if (f === 0) {
      subfolders.push({
        id: 'f0/broken',
        label: 'broken',
        children: () => later(null, LOAD_DELAY).then(() => {
          throw new Error('Access denied: /f0/broken');
        }),
      });
    }
    roots.push({id: `f${f}`, label: `Folder ${f}`, children: subfolders});
  }
  return roots;
}

export async function render(main) {
  ensureCss('list');
  ensureCss('tree');
  main.append(el('h1', null, 'Tree'));
  const intro = el('p');
  intro.innerHTML = 'A virtualized tree composed over <code>VirtualList</code>: expansion state and lazily ' +
    'loaded branches flatten into the row array the list renders, so 100k visible nodes cost the same as ' +
    '100k list rows. Click a twistie or use ← → to collapse/expand, ↑ ↓ to move, F2 or double-click to ' +
    `rename. Subfolder children resolve after ${LOAD_DELAY}ms — expand one to see the loading row; ` +
    'expand <code>Folder 0 › broken</code> to see a failed branch.';
  main.append(intro);

  const roots = makeRoots();
  const renameLog = el('p', 'u2-gallery-status', 'Last rename: —');
  const tree = new VirtualTree({
    onRename: (node, newLabel) => {
      renameLog.textContent = `Last rename: "${node.label}" → "${newLabel}"`;
      node.label = newLabel;
    },
  });
  tree.root.style.height = '360px';
  tree.setRoots(roots);
  main.append(tree.root);

  const readout = el('p', 'u2-gallery-status');
  tree.effect(() => {
    const node = tree.selectedNode.value;
    readout.textContent =
      `${(FOLDERS * SUBFOLDERS * LEAVES).toLocaleString()} nodes · selected = ${node ? node.label : 'none'}` +
      `${node && node.data ? ` (data: folder ${node.data.folder}, item ${node.data.leaf})` : ''}` +
      ` · expanded = ${tree.expanded.value.size} · live scopes = ${Scope.liveCount}`;
  });
  main.append(readout, renameLog);

  const disposed = el('p', 'u2-gallery-status');
  const row = el('div');
  row.append(
    button('Expand path to Folder 7 › Subfolder 3 › Item 42', () => tree.expandPath(['f7', 'f7/s3', 'f7/s3/l42'])),
    ' ',
    button('Collapse all', () => tree.expanded.value = new Set()),
    ' ',
    button('Dispose', () => {
      tree.dispose();
      disposed.textContent = `Disposed: live scopes = ${Scope.liveCount}. The inner VirtualList is a child ` +
        'component — the tree owns its dispose, so its effects, listeners and observers went with it.';
    }));
  main.append(row, disposed);
}
