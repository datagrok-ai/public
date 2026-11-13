import {test, testOwner} from './owned-test.js';
import assert from 'node:assert/strict';
import {register} from 'node:module';
import {fire, flush, resetDom} from './dom-shim.js';
import {layout} from './layout.mjs';
import {Scope} from '../src/core/scope.js';
import {VirtualTree} from '../src/components/collections/tree.js';
import {Project} from './platform-doubles.mjs';

register('./dg-stub.mjs', import.meta.url);
const grok = await import('datagrok-api/grok');
const DG = await import('datagrok-api/dg');
const {SpaceTree, recentSpaces, rememberSpace, entityName, spaceHasName, READ_ONLY} =
  await import('../src/dg/spaces/spaces.js');
const {SpaceInput} = await import('../src/dg/inputs/space-input.js');
const {saveToSpace, openFromSpace, SpaceBrowser} = await import('../src/dg/spaces/space-dialogs.js');
const {spaceSaveRow} = await import('../src/dg/spaces/space-save-row.js');

const space = (id, name) => new Project(name, {id, name, friendlyName: name, isSpace: true});
const item = (id, name) => new Project(name, {id, name, friendlyName: name, isSpace: false});

const CHEM = space('chem', 'Chemistry');
const KIN = space('kin', 'Kinase');
const Q4 = space('q4', 'Q4 screens');
const GPCR = space('gpcr', 'GPCR');
const ARCHIVE = space('arch', 'Archive');

/** A recording `grok.dapi.spaces` over `children` (subspaces by parent id, `null` for the roots) and `items`
 * (what each space holds); `taken` lists the names each space holds already. */
function fakeSpaces({children, items = {}, readOnly = [], taken = {}, linked = {}}) {
  const calls = {created: [], typesAsked: [], linkedAsked: 0};
  grok.dapi.spaces = {
    list: async () => children[null] ?? [],
    createRootSpace: async (name) => {
      calls.created.push([null, name]);
      return space(`new-${name}`, name);
    },
    id: (id) => ({
      children: {ofTypes: (types, includeLinked) => ({list: async () => {
        if (types === 'Project') {
          calls.linkedAsked += includeLinked ? 1 : 0;
          return [...children[id] ?? [], ...includeLinked ? linked[id] ?? [] : []];
        }
        calls.typesAsked.push(types);
        return items[id] ?? [];
      }})},
      subspaceExists: async (name) => (taken[id] ?? []).includes(name),
      addSubspace: async (name) => {
        calls.created.push([id, name]);
        return space(`new-${name}`, name);
      },
    }),
  };
  grok.dapi.permissions = {check: async (p) => !readOnly.includes(p.id)};
  return calls;
}

const TREE = {children: {null: [GPCR, CHEM, ARCHIVE], chem: [KIN], kin: [Q4]}, readOnly: ['arch']};

function ui(name, body) {
  test(name, async () => {
    const live = Scope.liveCount;
    try {
      await body();
    }
    finally {
      delete grok.dapi.spaces;
      delete grok.dapi.permissions;
      grok.userSettings.values = {};
      DG.ObjectHandler.defaultActions.length = 0;
      resetDom();
      await flush();
    }
    assert.equal(Scope.liveCount, live, 'live scopes back to baseline');
  });
}

async function until(cond, what, timeoutMs = 1500) {
  const start = Date.now();
  while (!cond()) {
    if (Date.now() - start > timeoutMs)
      assert.fail(`timed out waiting for ${what}`);
    await flush();
  }
}

const popup = () => document.querySelector('.u2-space-popup');
const treeRows = (root) => [...root.querySelectorAll('.u2-tree-row')];
const rowLabels = (root) => treeRows(root).map((r) => r.textContent.trim());
const rowOf = (root, label) => treeRows(root).find((r) => r.textContent.trim() === label);

/** Mounts a SpaceInput, opens it and waits for the tree rows. */
async function opened(options = {}) {
  const input = testOwner().build(() => new SpaceInput({label: 'Space', ...options}));
  document.body.append(input.root);
  fire(input.root.querySelector('.u2-space-box'), 'click');
  await until(() => popup() !== null && input.tree.peek() !== null, 'the spaces');
  layout(popup().querySelector('.u2-space-tree .u2-list'), {clientHeight: 200});
  fire(popup().querySelector('.u2-space-tree .u2-list'), 'scroll');
  await flush();
  return input;
}

ui('SpaceTree reads the roots, then each level of subspaces; paths, read-only reason, sorting and under()', async () => {
  const mine = space('mine', 'Tester');
  fakeSpaces({...TREE, children: {...TREE.children, null: [mine, ...TREE.children.null]}});
  grok.shell.dart.user = {...grok.shell.user, project: mine};
  const tree = await SpaceTree.load().finally(() => delete grok.shell.dart.user.project);
  assert.deepEqual(tree.roots.map((e) => e.project.friendlyName), ['Archive', 'Chemistry', 'GPCR'],
    'sorted by name, without the user\'s own namespace');
  assert.deepEqual(tree.path(Q4), ['Chemistry', 'Kinase', 'Q4 screens']);
  assert.equal(tree.get('arch').readOnly, READ_ONLY);
  assert.equal(tree.get('kin').readOnly, null);
  assert.deepEqual(tree.under('chem').map((e) => e.project.id), ['chem', 'kin', 'q4']);
  assert.deepEqual(tree.all().map((e) => e.project.id), ['arch', 'chem', 'kin', 'q4', 'gpcr']);
  assert.deepEqual(tree.path(space('elsewhere', 'Elsewhere')), ['Elsewhere'], 'a space outside the tree: its name');
});

test('entityName: the name the platform derives, which the clash check compares', () => {
  assert.equal(entityName('EGFR SAR review'), 'EGFRSARReview');
  assert.equal(entityName('  2024 kinase-panel_v2 '), 'KinasePanelV2', 'leading digits dropped, words capitalized');
  assert.equal(entityName('Q4 screens'), 'Q4Screens');
  assert.equal(entityName('—'), '');
});

ui('SpaceInput: empty shows the placeholder; a click on a space picks it and shows its path', async () => {
  fakeSpaces(TREE);
  const input = await opened({placeholder: 'Save to a space (optional)'});
  assert.equal(input.root.querySelector('.u2-space-text').textContent, 'Save to a space (optional)');
  assert.deepEqual(rowLabels(popup()), ['Archive', 'Chemistry', 'GPCR']);
  fire(rowOf(popup(), 'Chemistry').querySelector('.u2-tree-twistie'), 'click');
  await flush();
  assert.ok(popup() !== null, 'the twistie opens the branch, it does not pick');
  fire(rowOf(popup(), 'Kinase'), 'click');
  await flush();
  assert.equal(input.value, KIN);
  assert.equal(popup(), null, 'a pick closes the popup');
  assert.equal(input.root.querySelector('.u2-space-text').textContent, 'Chemistry / Kinase');
  input.dispose();
});

ui('SpaceInput: a read-only space is greyed with the reason and cannot be picked', async () => {
  fakeSpaces(TREE);
  const input = await opened();
  const archive = rowOf(popup(), 'Archive');
  assert.ok(archive.classList.contains('u2-tree-row-disabled'));
  assert.equal(archive.title, READ_ONLY);
  fire(archive, 'click');
  await flush();
  assert.equal(input.value, null);
  assert.ok(popup() !== null);
  input.dispose();
});

ui('SpaceInput: the search keeps the matching spaces with their ancestors, the recent list goes away', async () => {
  fakeSpaces(TREE);
  rememberSpace(GPCR);
  const input = await opened();
  assert.deepEqual([...popup().querySelectorAll('.u2-space-recent')].map((r) => r.firstChild.nextSibling.textContent),
    ['GPCR'], 'the space saved to last');
  const box = popup().querySelector('.u2-space-search input');
  box.value = 'q4';
  fire(box, 'input');
  await flush();
  fire(popup().querySelector('.u2-space-tree .u2-list'), 'scroll');
  await flush();
  assert.deepEqual(rowLabels(popup()), ['Chemistry', 'Kinase', 'Q4 screens']);
  assert.ok(popup().querySelector('.u2-space-recents').hidden);
  input.dispose();
});

ui('SpaceInput: a recent space picks in one click', async () => {
  fakeSpaces(TREE);
  rememberSpace(Q4);
  rememberSpace(KIN);
  assert.deepEqual(recentSpaces(), ['kin', 'q4'], 'the latest first');
  const input = await opened();
  fire(popup().querySelectorAll('.u2-space-recent')[1], 'click');
  await flush();
  assert.equal(input.value, Q4);
  input.dispose();
});

ui('SpaceInput: × clears the value without opening the tree; a root-limited input offers no × and no other roots',
  async () => {
    fakeSpaces(TREE);
    const input = testOwner().build(() => new SpaceInput({label: 'Space', value: KIN}));
    document.body.append(input.root);
    const clear = input.root.querySelector('.u2-space-clear');
    assert.equal(clear.hidden, false);
    fire(clear, 'click');
    await flush();
    assert.equal(input.value, null);
    assert.equal(popup(), null, 'the tree stays closed');
    assert.equal(clear.hidden, true, 'nothing to clear');
    input.dispose();
    resetDom();

    const limited = await opened({root: CHEM, value: KIN, nullable: false});
    assert.equal(limited.root.querySelector('.u2-space-clear').hidden, true);
    assert.deepEqual(rowLabels(popup()), ['Chemistry', 'Kinase', 'Q4 screens'], 'Chemistry alone, open to the value');
    limited.dispose();
  });

/** Clicks a "+" (the heading's, or a row's) and names the new space in the inline box. */
async function newSpace(plus, name) {
  fire(plus, 'click');
  await until(() => popup()?.querySelector('.u2-tree-rename') != null, 'the inline name box');
  const box = popup().querySelector('.u2-tree-rename');
  assert.equal(box.value, 'New space');
  box.value = name;
  fire(box, 'keydown', {key: 'Enter'});
}
const topPlus = () => popup().querySelector('.u2-space-head button');
const rowPlus = (label) => rowOf(popup(), label).querySelector('.u2-space-add');

ui('SpaceInput: the heading\'s "+" creates a top-level space at once; its row\'s "+" nests one inside it', async () => {
  const calls = fakeSpaces(TREE);
  const input = await opened();
  await newSpace(topPlus(), 'Screens');
  await until(() => input.value?.friendlyName === 'Screens', 'the new space picked');
  assert.deepEqual(calls.created, [[null, 'Screens']], 'created on naming, not on save');
  assert.ok(popup() !== null, 'the popup stays open');
  await newSpace(rowPlus('Screens'), 'Q1');
  await until(() => input.value?.friendlyName === 'Q1', 'the nested space picked');
  assert.deepEqual(calls.created, [[null, 'Screens'], ['new-Screens', 'Q1']]);
  assert.equal(input.root.querySelector('.u2-space-text').textContent, 'Screens / Q1');
  input.dispose();
});

ui('SpaceInput: a row\'s "+" creates inside that space without picking it first; a read-only space has none',
  async () => {
    const calls = fakeSpaces(TREE);
    const input = await opened();
    assert.equal(rowPlus('Archive'), null, 'read-only');
    await newSpace(rowPlus('GPCR'), 'Selectivity');
    await until(() => input.value?.friendlyName === 'Selectivity', 'the new space picked');
    assert.deepEqual(calls.created, [['gpcr', 'Selectivity']]);
    assert.equal(input.root.querySelector('.u2-space-text').textContent, 'GPCR / Selectivity');
    input.dispose();
  });

/** Clicks a "+" and ends the inline name box with `key`, the offered name left as it is. */
async function newSpaceKey(plus, key) {
  fire(plus, 'click');
  await until(() => popup()?.querySelector('.u2-tree-rename') != null, 'the inline name box');
  fire(popup().querySelector('.u2-tree-rename'), 'keydown', {key});
}

ui('SpaceInput: Enter on the offered name creates "New space"; a click away drops the row and a click picks again',
  async () => {
    const calls = fakeSpaces(TREE);
    const input = await opened();
    await newSpaceKey(topPlus(), 'Enter');
    await until(() => input.value?.friendlyName === 'New space', 'the offered name kept');
    assert.deepEqual(calls.created, [[null, 'New space']]);
    fire(rowPlus('GPCR'), 'click');
    await until(() => popup().querySelector('.u2-tree-rename') !== null, 'the inline name box');
    popup().querySelector('.u2-space-search input').focus();
    await until(() => popup().querySelector('.u2-tree-rename') === null &&
      treeRows(popup()).filter((r) => r.textContent.trim() === 'New space').length === 1, 'the pending row gone');
    assert.deepEqual(calls.created, [[null, 'New space']], 'nothing more created');
    fire(rowOf(popup(), 'GPCR'), 'click');
    await flush();
    assert.equal(input.value, GPCR, 'a click picks again');
    input.dispose();
  });

ui('SpaceInput: Escape while naming closes the popup with nothing created and the value kept', async () => {
  const calls = fakeSpaces(TREE);
  const input = await opened({value: KIN});
  await newSpaceKey(rowPlus('GPCR'), 'Escape');
  await until(() => popup() === null, 'the popup closed');
  assert.deepEqual(calls.created, []);
  assert.equal(input.value, KIN);
  fire(input.root.querySelector('.u2-space-box'), 'click');
  await until(() => popup() !== null && input.tree.peek() !== null, 'the popup again');
  layout(popup().querySelector('.u2-space-tree .u2-list'), {clientHeight: 200});
  fire(popup().querySelector('.u2-space-tree .u2-list'), 'scroll');
  await until(() => rowOf(popup(), 'GPCR') !== undefined, 'the spaces again');
  assert.equal(popup().querySelector('.u2-tree-rename'), null, 'no row left to name');
  input.dispose();
});

ui('SpaceInput: while a space is being created "+" waits, and a space picked meanwhile stays picked', async () => {
  const calls = fakeSpaces(TREE);
  let finish;
  grok.dapi.spaces.createRootSpace = (name) => new Promise((resolve) => finish = () => resolve(space('slow', name)));
  const input = await opened();
  await newSpace(topPlus(), 'Slow');
  await until(() => finish !== undefined, 'the create sent');
  fire(rowPlus('GPCR'), 'click');
  await flush();
  assert.equal(popup().querySelector('.u2-tree-rename'), null, 'no second name box while creating');
  fire(rowOf(popup(), 'GPCR'), 'click');
  await flush();
  assert.equal(input.value, GPCR);
  finish();
  await until(() => input.tree.peek().get('slow') !== null, 'the new space in the tree');
  assert.equal(input.value, GPCR, 'the pick kept');
  assert.deepEqual(calls.created, []);
  input.dispose();
});

ui('SpaceInput: a space the platform refuses to create says why, and nothing is picked', async () => {
  const calls = fakeSpaces(TREE);
  grok.dapi.spaces.createRootSpace = async () => {
    throw new Error('Not allowed');
  };
  const input = await opened();
  await newSpace(topPlus(), 'Top');
  await until(() => popup().querySelector('.u2-space-failed').textContent !== '', 'the reason');
  assert.equal(popup().querySelector('.u2-space-failed').textContent, 'Could not create "Top": Not allowed');
  assert.equal(input.value, null);
  assert.deepEqual(calls.created, []);
  input.dispose();
});

ui('VirtualTree.rename(id) opens the inline rename on that node; Enter reports the new label', async () => {
  const renamed = [];
  const tree = testOwner().build(() => new VirtualTree({onRename: (node, label) => renamed.push([node.id, label])}));
  document.body.append(tree.root);
  layout(tree.root.querySelector('.u2-list'), {clientHeight: 200});
  tree.setRoots([{id: 'a', label: 'Alpha', children: [{id: 'a1', label: 'Alpha one'}]}, {id: 'b', label: 'Beta'}]);
  assert.equal(await tree.rename('a1'), true);
  const box = tree.root.querySelector('.u2-tree-rename');
  assert.ok(box !== null, 'the rename box is open');
  assert.equal(box.value, 'Alpha one');
  assert.equal(tree.selectedNode.peek().id, 'a1');
  box.value = 'First';
  fire(box, 'keydown', {key: 'Enter'});
  assert.deepEqual(renamed, [['a1', 'First']]);
  assert.equal(await tree.rename('nope'), false);
  tree.dispose();
  const plain = testOwner().build(() => new VirtualTree());
  assert.equal(await plain.rename('a'), false, 'no onRename, no rename');
  plain.dispose();
});

const dialog = () => document.querySelector('.u2-dialog');
const dialogButton = (text) => [...(dialog()?.querySelectorAll('button') ?? [])].find((b) => b.textContent.trim() === text);
const off = (button) => button?.getAttribute('aria-disabled') === 'true';

async function type(input, text) {
  input.value = text;
  fire(input, 'input');
  await flush();
}

ui('saveToSpace(): SAVE is off while the space holds the name, back on for another name; the save gets the space', async () => {
  const calls = fakeSpaces({...TREE, taken: {kin: ['EGFRSAR']}});
  const saved = [];
  const result = saveToSpace({name: 'EGFR SAR', description: '', space: KIN,
    save: (into, values) => saved.push([into.id, values.name])});
  await until(() => off(dialogButton('SAVE')), 'the clash to switch SAVE off');
  assert.match(dialog().querySelector('.u2-space-clash').textContent, /Chemistry \/ Kinase already has "EGFR SAR"/);
  const name = [...dialog().querySelectorAll('input')].find((i) => i.value === 'EGFR SAR');
  await type(name, 'EGFR SAR v2');
  await until(() => dialogButton('SAVE') !== undefined && !off(dialogButton('SAVE')), 'SAVE back on');
  assert.equal(dialog().querySelector('.u2-space-clash').textContent, '');
  fire(dialogButton('SAVE'), 'click');
  assert.equal(await result, KIN);
  assert.deepEqual(saved, [['kin', 'EGFR SAR v2']]);
  assert.deepEqual(calls.created, []);
  assert.deepEqual(recentSpaces(), ['kin'], 'remembered for the Recent list');
});

ui('saveToSpace(): SAVE is off, saying why, until there is a space and a name', async () => {
  fakeSpaces(TREE);
  const saved = [];
  const result = saveToSpace({save: () => saved.push(1)});
  await until(() => off(dialogButton('SAVE')), 'SAVE off with no space');
  assert.equal(dialogButton('SAVE').getAttribute('aria-description'), 'Pick a space');
  fire(dialog().querySelector('.u2-space-box'), 'click');
  await until(() => popup() !== null && rowOf(popup(), 'GPCR') !== undefined, 'the spaces');
  fire(rowOf(popup(), 'GPCR'), 'click');
  await until(() => dialogButton('SAVE').getAttribute('aria-description') === 'Enter a name', 'SAVE off with no name');
  await type([...dialog().querySelectorAll('input')].find((i) => i.value === ''), 'Run 1');
  await until(() => !off(dialogButton('SAVE')), 'SAVE on');
  fire(dialogButton('CANCEL'), 'click');
  assert.equal(await result, null);
  assert.deepEqual(saved, []);
});

ui('saveToSpace(): CANCEL answers null and saves nothing', async () => {
  fakeSpaces(TREE);
  const saved = [];
  const result = saveToSpace({name: 'X', space: KIN, save: () => saved.push(1)});
  await until(() => dialog()?.textContent.includes('Chemistry / Kinase'), 'the dialog with the spaces read');
  fire(dialogButton('CANCEL'), 'click');
  assert.equal(await result, null);
  assert.deepEqual(saved, []);
});

const ITEMS = {
  chem: [item('p1', 'SAR review')],
  kin: [item('p2', 'Kinase panel'), item('r1', 'ADMET run'), space('kin-sub', 'Subspace')],
  q4: [item('p3', 'Q4 hits')],
  gpcr: [item('p4', 'GPCR assay')],
};

/** The names the item list shows, after laying it out. */
function itemNames(root) {
  const list = root.querySelector('.u2-list.u2-space-items');
  layout(list, {clientHeight: 300});
  fire(list, 'scroll');
  return [...list.querySelectorAll('.u2-space-item-name')].map((n) => n.textContent).sort();
}

async function shows(root, expected) {
  const want = JSON.stringify([...expected].sort());
  await until(() => JSON.stringify(itemNames(root)) === want, `items ${want} (have ${JSON.stringify(itemNames(root))})`);
}

async function browsed(options = {}) {
  const browser = testOwner().build(() => new SpaceBrowser(options));
  document.body.append(browser.root);
  await until(() => browser.spaces.peek() !== null, 'the spaces');
  layout(browser.root.querySelector('.u2-space-browser-tree .u2-list'), {clientHeight: 300});
  fire(browser.root.querySelector('.u2-space-browser-tree .u2-list'), 'scroll');
  return browser;
}

ui('SpaceBrowser: "All spaces" lists what every space holds, with where it lives; a space narrows to it and below', async () => {
  const calls = fakeSpaces({...TREE, items: ITEMS});
  const browser = await browsed();
  await shows(browser.root, ['SAR review', 'Kinase panel', 'ADMET run', 'Q4 hits', 'GPCR assay']);
  const row = [...browser.root.querySelectorAll('.u2-space-item')].find((r) => r.textContent.includes('Q4 hits'));
  assert.ok(row.textContent.includes('Chemistry / Kinase / Q4 screens'), 'the item shows its space path');
  const search = browser.root.querySelector('.u2-space-search input');
  assert.equal(search.placeholder, 'Search all spaces');
  assert.deepEqual([...new Set(calls.typesAsked)], [''], 'every type by default');

  fire(rowOf(browser.root, 'Chemistry'), 'click');
  await shows(browser.root, ['SAR review', 'Kinase panel', 'ADMET run', 'Q4 hits']);
  assert.equal(search.placeholder, 'Search Chemistry');

  await type(search, 'kinase');
  await shows(browser.root, ['Kinase panel']);
  await type(search, 'zzz');
  await until(() => browser.root.querySelector('.u2-space-empty-state').textContent === 'Nothing matches "zzz".',
    'the empty state');
  browser.dispose();
});

ui('SpaceBrowser: types go to the server, the filter keeps what it accepts', async () => {
  const calls = fakeSpaces({...TREE, items: ITEMS});
  const browser = await browsed({types: ['Project', 'FuncCall'], filter: (e) => e.friendlyName !== 'ADMET run'});
  await shows(browser.root, ['SAR review', 'Kinase panel', 'Q4 hits', 'GPCR assay']);
  assert.deepEqual([...new Set(calls.typesAsked)], ['Project,FuncCall']);
  browser.dispose();
});

ui('openFromSpace(): OPEN is off until an item is picked, then runs its default action and answers it', async () => {
  fakeSpaces({...TREE, items: ITEMS});
  const result = openFromSpace();
  await until(() => dialog() !== null && itemNames(dialog()).length === 5, 'the items');
  assert.ok(off(dialogButton('OPEN')));
  const row = [...dialog().querySelectorAll('.u2-space-items .u2-list-row')]
    .find((r) => r.textContent.includes('GPCR assay'));
  fire(row, 'click');
  await until(() => !off(dialogButton('OPEN')), 'OPEN on');
  fire(dialogButton('OPEN'), 'click');
  assert.deepEqual((await result).map((e) => e.id), ['p4']);
  assert.deepEqual(DG.ObjectHandler.defaultActions.map((e) => e.id), ['p4']);
});

/** A Space row mounted on the page, as the Save project dialog mounts it. */
function savingRow(home, name, projectId = null) {
  const row = spaceSaveRow(home, name, projectId);
  document.body.append(row.root);
  return row;
}
const rowText = (row) => row.root.querySelector('.u2-status-line-text').textContent;

ui('spaceSaveRow: a new project stays in the user\'s projects, or goes into the picked space with no Move or Link',
  async () => {
    fakeSpaces(TREE);
    const row = savingRow(null, 'Kinase panel');
    try {
      assert.deepEqual(await row.resolveTarget(), {spaceId: null, link: false, changed: false});
      row.target = KIN;
      await until(() => rowText(row) === 'Shared with whoever can see "Kinase".', 'the hint, the space quoted');
      assert.equal(row.root.querySelector('[data-u2="choice-chip"]').hidden, true, 'no Move or Link');
      assert.deepEqual(await row.resolveTarget(), {spaceId: 'kin', link: false, changed: true});
      assert.deepEqual(recentSpaces()[0], 'kin', 'remembered for the Recent list');
    }
    finally {
      row.detach();
    }
  });

ui('spaceSaveRow: a saved project stays in its space, moves to another, or is linked from it', async () => {
  fakeSpaces(TREE);
  const row = savingRow(KIN, 'Kinase panel');
  try {
    assert.deepEqual(await row.resolveTarget(), {spaceId: 'kin', link: false, changed: false});
    row.target = GPCR;
    await until(() => rowText(row) === 'Moves from "Kinase" to "GPCR".', 'the move hint');
    assert.deepEqual(await row.resolveTarget(), {spaceId: 'gpcr', link: false, changed: true});
    row.mode = 'link';
    await until(() => rowText(row) === 'Stays in "Kinase"; "GPCR" gets a link to it.', 'the link hint');
    assert.deepEqual(await row.resolveTarget(), {spaceId: 'gpcr', link: true, changed: true});
    row.target = null;
    await until(() => rowText(row) === 'Moves out of "Kinase" into your projects.', 'the hint for no space');
    assert.deepEqual(await row.resolveTarget(), {spaceId: null, link: false, changed: true});
  }
  finally {
    row.detach();
  }
});

ui('spaceSaveRow: a space holding a project of that name is a problem until the name changes; a link is not',
  async () => {
    fakeSpaces({...TREE, taken: {gpcr: ['KinasePanel']}});
    const row = savingRow(KIN, 'Kinase panel');
    let changes = 0;
    row.onChanged = () => changes++;
    try {
      row.target = GPCR;
      await until(() => row.problem !== null, 'the clash');
      assert.equal(row.problem, '"GPCR" already has a project named "Kinase panel". ' +
        'Rename the project or pick another space.');
      assert.ok(changes > 0, 'the dialog told');
      row.mode = 'link';
      await until(() => row.problem === null, 'a link checks no name');
      row.mode = 'move';
      await until(() => row.problem !== null, 'the clash again');
      row.setName('Kinase panel 2');
      await until(() => row.problem === null, 'another name accepted');
    }
    finally {
      row.detach();
    }
  });

ui('spaceHasName: the saved project, which the space links, does not count; another of its name does', async () => {
  const saved = new Project('Kinase panel', {id: 'p1', name: 'KinasePanel', friendlyName: 'Kinase panel'});
  const calls = fakeSpaces({...TREE, taken: {gpcr: ['KinasePanel']}, linked: {gpcr: [saved]}});
  assert.equal(await spaceHasName(GPCR, 'Kinase panel', 'p1'), false, 'its own link');
  assert.equal(await spaceHasName(GPCR, 'Kinase panel', 'p2'), true, 'another project');
  assert.equal(await spaceHasName(GPCR, 'Kinase panel'), true, 'a new project');
  const row = savingRow(KIN, 'Kinase panel', 'p2');
  try {
    row.target = GPCR;
    const asked = calls.linkedAsked;
    await until(() => row.problem !== null, 'the clash');
    assert.ok(calls.linkedAsked > asked, 'the row passed its project on');
  }
  finally {
    row.detach();
  }
});
