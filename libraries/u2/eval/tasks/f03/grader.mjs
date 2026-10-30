import {buttons, clickOk, controlOf, css, errors, expect, flush, liveScopes, mount, named, need, row, until, wait}
  from '../../lib/harness.mjs';

const ANN = {login: 'ann', friendlyName: 'Ann Lee'};
const BOB = {login: 'bob', friendlyName: 'Bob Stone'};
const CAMPAIGNS = [
  {code: 'KHL-12', name: 'PI3K quick wins', status: 'In Progress', kickoff: new Date(2026, 5, 29), owner: ANN,
    reference: {id: 'V012', smiles: 'CCO'}},
  {code: 'KHL-9', name: 'BRD4 follow-up', status: 'Submitted', kickoff: new Date(2026, 5, 29), owner: BOB,
    reference: {id: 'V009', smiles: 'CC(=O)Oc1ccccc1C(=O)O'}},
  {code: 'KHL-3', name: 'Kinase hits', status: 'In Progress', kickoff: new Date(2026, 6, 1), owner: ANN,
    reference: {id: 'V003', smiles: 'c1ccccc1'}},
  {code: 'KHL-4', name: 'GPCR panel', status: 'Submitted', kickoff: new Date(2026, 9, 2), owner: BOB,
    reference: {id: 'V004', smiles: 'CCN'}},
];
const DAYS = ['2026/6/29', '2026/6/29', '2026/7/1', '2026/10/2'];
const CHOSEN = CAMPAIGNS[1];

function deferred() {
  let resolve;
  const promise = new Promise((r) => resolve = r);
  return {promise, resolve};
}

async function build(mod) {
  const pending = deferred();
  const opened = [];
  const base = liveScopes();
  const home = need(mod, 'campaignsHome')(() => pending.promise, (c) => void opened.push(c));
  expect(typeof home?.dispose === 'function', 'campaignsHome() did not return a control');
  const root = mount(home);
  await flush();
  const resolve = () => pending.resolve(CAMPAIGNS.map((c) => ({...c})));
  return {home, root, base, opened, resolve};
}

const campaignRows = (root) => root.querySelectorAll('tr[data-u2-item]')
  .filter((tr) => !tr.dataset.u2Item.startsWith('section:'));
const saysLoading = (root) => /loading/i.test(root.textContent) || [root, ...root.querySelectorAll('*')]
  .some((el) => el.getAttribute('aria-busy') === 'true' || el.getAttribute('role') === 'status');

async function loaded(mod) {
  const b = await build(mod);
  b.resolve();
  expect(await until(() => campaignRows(b.root).length === CAMPAIGNS.length, 2000),
    `the campaigns never show after the load resolves (${campaignRows(b.root).length} rows)` +
    (errors().length ? `; logged: ${errors()[0]}` : ''));
  return b;
}

function catalog(root) {
  const el = root.querySelector('[data-u2="catalog-table"]');
  const table = el === null ? undefined : controlOf(el);
  expect(typeof table?.groupBy === 'function', 'the home shows no u2 CatalogTable');
  return {el, table};
}

const text = (el) => el?.textContent.replace(/\s+/g, ' ').trim() ?? '';
const firstCell = (tr) => tr.querySelector('td, th') ?? tr;
const rgb = (c) => `rgb(${c.join(', ')})`;
const differ = (a, b) => a.some((x, i) => Math.abs(x - b[i]) >= 2);
const px = (v) => /^([\d.]+)px$/.exec(v ?? '')?.[1];
const INTERNAL = /\.(u2|d4|ui|grok)-|\[data-u2/;

/** Each group header's label with the texts its rows show in `column`; ungrouped, one group of every row. */
async function groups(root, grouping, column) {
  const {el, table} = catalog(root);
  table.groupBy(grouping);
  await flush();
  const at = el.querySelectorAll('thead th').findIndex((th) => text(th).startsWith(column));
  expect(at >= 0, `no ${column} column`);
  const out = grouping === null ? [{label: null, cells: []}] : [];
  for (const item of table.getItems({all: true}).list) {
    if (item.key.startsWith('section:'))
      out.push({label: item.label, cells: []});
    else {
      const tr = el.querySelectorAll('tr').find((r) => r.dataset.u2Item === item.key);
      out.at(-1)?.cells.push(text(tr?.children[at]));
    }
  }
  expect(out.length > 0, `grouping by ${grouping} shows no group headers: every row fell in one group`);
  return out;
}

async function headersRead(mod, grouping, column, expected) {
  const {root} = await loaded(mod);
  const found = await groups(root, grouping, column);
  const labels = found.map((g) => g.label).join(', ');
  expect(found.length === expected, `grouping by ${grouping} makes ${found.length} groups (${labels}), not ${expected}`);
  for (const g of found) {
    for (const cell of g.cells)
      expect(cell === g.label, `the ${grouping} header reads "${g.label}" over rows whose ${column} reads "${cell}"`);
  }
}

/** Clicks the chosen campaign's name and answers its details: the drawing and what holds it, the list left out. */
async function showDetails(root) {
  const tr = campaignRows(root).find((r) => text(r).includes(CHOSEN.name));
  expect(tr !== undefined, `no row of ${CHOSEN.name}`);
  const name = named(tr, new RegExp(CHOSEN.name)) ?? tr.querySelectorAll('td').find((td) => text(td) === CHOSEN.name);
  await clickOk(name, `the name ${CHOSEN.name}`);
  await flush();
  const drawing = root.querySelectorAll('.eval-molecule').find((d) => d.dataset.smiles === CHOSEN.reference.smiles);
  expect(drawing !== undefined, `a click on ${CHOSEN.name} shows no drawing of its reference compound`);
  const list = catalog(root).el;
  let details = drawing;
  while (details.parentElement && details.parentElement !== root && !details.parentElement.contains(list))
    details = details.parentElement;
  return {drawing, details};
}

export const checks = {
  'loads': async (mod) => {
    const b = await build(mod);
    expect(await until(() => saysLoading(b.root), 500), 'nothing shows that it is loading');
    b.resolve();
    expect(await until(() => campaignRows(b.root).length === CAMPAIGNS.length, 2000),
      `the campaigns never show after the load resolves (${campaignRows(b.root).length} rows)`);
  },
  'groups': async (mod) => {
    const {root} = await loaded(mod);
    const sections = root.querySelectorAll('[data-u2-item]').filter((el) => el.dataset.u2Item.startsWith('section:'));
    const labels = sections.map((s) => text(s));
    expect(sections.length === 2, `${sections.length} group headers (${labels.join(' / ')}), not one per status`);
    expect(/In Progress.*2/.test(labels[0]) && /Submitted.*2/.test(labels[1]),
      `the group headers read ${labels.join(' / ')}, not each status with its count`);
  },
  'hover-distinct': async (mod) => {
    const {root} = await loaded(mod);
    const sections = root.querySelectorAll('[data-u2-item]').filter((el) => el.dataset.u2Item.startsWith('section:'));
    const rows = campaignRows(root);
    expect(sections.length > 0 && rows.length > 0, 'no group header and row to compare');
    const header = css.background(firstCell(sections[0]));
    const cells = rows[0].querySelectorAll('td, th');
    for (const [i, cell] of (cells.length ? cells : [rows[0]]).entries()) {
      const hovered = css.background(cell, {hover: cell});
      expect(differ(hovered, css.background(cell)),
        `cell ${i + 1} of the row under the pointer looks as at rest (${rgb(hovered)})`);
      expect(differ(hovered, header), `cell ${i + 1} of the hovered row has the group header's colour (${rgb(hovered)})`);
    }
  },
  'cells-format': async (mod) => {
    const {root} = await loaded(mod);
    const shown = (await groups(root, null, 'Kickoff'))[0].cells.sort();
    expect(JSON.stringify(shown) === JSON.stringify([...DAYS].sort()),
      `Kickoff cells read ${shown.join(', ')}, not the platform's year/month/day`);
    const owners = (await groups(root, null, 'Owner'))[0].cells;
    expect(owners.every((o) => o === ANN.friendlyName || o === BOB.friendlyName), `Owner cells read ${owners.join(', ')}`);
  },
  'kickoff-headers': (mod) => headersRead(mod, 'Kickoff', 'Kickoff', 3),
  'owner-headers': (mod) => headersRead(mod, 'Owner', 'Owner', 2),
  'details': async (mod) => {
    const {root, opened} = await loaded(mod);
    const {details} = await showDetails(root);
    for (const [caption, value] of [['Code', CHOSEN.code], ['Name', CHOSEN.name], ['Owner', BOB.friendlyName]]) {
      const r = row(caption, details);
      expect(r.text === value, `the details' ${caption} shows "${r.text}", not "${value}"`);
    }
    const open = buttons(details, /^open$/)[0] ?? buttons(root, /^open$/)[0];
    expect(open !== undefined, 'no OPEN button');
    await clickOk(open, 'OPEN');
    expect(opened.length === 1 && opened[0].code === CHOSEN.code,
      `open got ${JSON.stringify(opened.map((c) => c?.code))}, not ${CHOSEN.code}`);
  },
  'structure-in-value': async (mod) => {
    const {root} = await loaded(mod);
    const {drawing, details} = await showDetails(root);
    const structure = row('Structure', details);
    expect(structure.el.contains(drawing), 'the drawing is not in the Structure row');
    expect(drawing.closest('[data-u2-part="options"]') === null,
      'the drawing sits on the row\'s options rail (its side icons), not as its value');
  },
  'structure-has-room': async (mod) => {
    const {root} = await loaded(mod);
    const {drawing, details} = await showDetails(root);
    const structure = row('Structure', details);
    for (let el = drawing.parentElement; el && el !== structure.el.parentElement; el = el.parentElement) {
      for (const prop of ['height', 'max-height']) {
        const h = px(css.computed(el, prop));
        expect(h === undefined || +h >= 100, `the 100px drawing is held in a box ${h}px high (${el.className})`);
      }
    }
  },
  'dispose-while-loading': async (mod) => {
    const {home, root, base, resolve} = await build(mod);
    home.dispose();
    resolve();
    await wait(300);
    expect(liveScopes() === base, `${liveScopes() - base} scopes were built for a home closed while it loaded`);
    expect(campaignRows(root).length === 0, 'rows were built into the home after it was closed');
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
  'dispose-releases': async (mod) => {
    const {home, root, base} = await loaded(mod);
    await showDetails(root);
    home.dispose();
    await flush();
    expect(liveScopes() === base, `${liveScopes() - base} scopes outlive the closed home`);
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
  'package-css-clean': async (mod) => {
    await loaded(mod);
    for (const rule of css.parse(css.packageCss())) {
      for (const part of rule.parts)
        expect(!INTERNAL.test(part.selector), `package CSS styles u2's insides: ${part.selector}`);
      const forced = rule.decls.find((d) => d.important);
      expect(forced === undefined, `package CSS forces ${forced?.prop} with !important`);
    }
  },
};
