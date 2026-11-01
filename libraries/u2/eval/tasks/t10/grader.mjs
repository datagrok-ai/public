import {controlOf, expect, flush, mount, need} from '../../lib/harness.mjs';

const ANN = {login: 'ann', friendlyName: 'Ann Lee'};
const BOB = {login: 'bob', friendlyName: 'Bob Stone'};
const CAMPAIGNS = [
  {code: 'KHL-1', name: 'Kinase hits', kickoff: new Date(2026, 5, 29), owner: ANN},
  {code: 'KHL-2', name: 'PI3K quick wins', kickoff: new Date(2026, 5, 29), owner: BOB},
  {code: 'KHL-3', name: 'BRD4 follow-up', kickoff: new Date(2026, 6, 1), owner: ANN},
  {code: 'KHL-4', name: 'GPCR panel', kickoff: new Date(2026, 9, 2), owner: BOB},
];
const DAYS = {'KHL-1': '2026/6/29', 'KHL-2': '2026/6/29', 'KHL-3': '2026/7/1', 'KHL-4': '2026/10/2'};

async function build(mod) {
  const table = need(mod, 'campaignCatalog')(CAMPAIGNS.map((c) => ({...c})));
  const root = mount(table);
  await flush();
  expect(typeof table.groupBy === 'function', 'campaignCatalog() did not return a CatalogTable');
  return {table, root};
}

const text = (el) => el?.textContent.replace(/\s+/g, ' ').trim() ?? '';

/** Each group header's label with the texts its rows show in `column`; ungrouped, one group of every row. */
async function groups(table, root, grouping, column) {
  table.groupBy(grouping);
  await flush();
  const at = root.querySelectorAll('thead th').findIndex((th) => text(th).startsWith(column));
  expect(at >= 0, `no ${column} column`);
  const out = grouping === null ? [{label: null, cells: []}] : [];
  for (const item of controlOf(root).getItems({all: true}).list) {
    if (item.key.startsWith('section:'))
      out.push({label: item.label, cells: []});
    else {
      const tr = root.querySelectorAll('tr').find((r) => r.dataset.u2Item === item.key);
      out.at(-1)?.cells.push(text(tr?.children[at]));
    }
  }
  expect(out.length > 0, `grouping by ${grouping} shows no group headers: every row fell in one group`);
  return out;
}

async function headersRead(mod, grouping, column, expected) {
  const {table, root} = await build(mod);
  const found = await groups(table, root, grouping, column);
  const labels = found.map((g) => g.label).join(', ');
  expect(found.length === expected,
    `grouping by ${grouping} makes ${found.length} groups (${labels}), not ${expected}`);
  for (const g of found) {
    for (const cell of g.cells)
      expect(cell === g.label, `the ${grouping} header reads "${g.label}" over rows whose ${column} reads "${cell}"`);
  }
}

export const checks = {
  'cells-format': async (mod) => {
    const {table, root} = await build(mod);
    const shown = (await groups(table, root, null, 'Kickoff'))[0].cells.sort();
    expect(JSON.stringify(shown) === JSON.stringify(Object.values(DAYS).sort()),
      `Kickoff cells read ${shown.join(', ')}, not the platform's year/month/day`);
    const owners = (await groups(table, root, null, 'Owner'))[0].cells;
    expect(owners.every((o) => o === ANN.friendlyName || o === BOB.friendlyName),
      `Owner cells read ${owners.join(', ')}`);
  },
  'kickoff-headers': (mod) => headersRead(mod, 'Kickoff', 'Kickoff', 3),
  'owner-headers': (mod) => headersRead(mod, 'Owner', 'Owner', 2),
};
