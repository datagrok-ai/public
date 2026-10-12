/* eslint-disable max-len */
/* The substructure filter's query round trip (crux-sketch spike query-roundtrip; the product owner's report of
 * 2026-10-07: "every complex query feature in sketcher is correctly saved in datagrok filter, and when sketcher is
 * re-opened query should be same, and also, the things that are being filtered should be correct ... should be checked
 * on 3 datasets: mol1K, SPGI, smiles.csv"). In a table view's filter panel, as the user works: the card's sketcher opened
 * in its dialog, a molecule drawn (Filter as you draw on, as an account has it by default), then a query feature added by
 * the sketcher's own controls (Crux: its toolbar, its canvas menu, its label editor and its dialogs, by the events a user's
 * pointer and keys give; Ketcher: its structure as its tools leave it, loaded as KET, Ketcher's own format), OK. Then:
 * the rows the filter passes are exactly those RDKit matches with the feature's intended SMARTS; the sketcher reopened
 * from the card holds the same query (the molblock the host reads, compared by what the search makes of it, and the
 * sketcher's own SMARTS); and switched to the other sketcher and back, it holds it still. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';

// RDKit's own generics (QueryOps), the meanings Crux's generics have
const RD_M = '[!#0&!#2&!#5&!#6&!#7&!#8&!#9&!#10&!#14&!#15&!#16&!#17&!#18&!#33&!#34&!#35&!#36&!#52&!#53&!#54&!#85&!#86&!#1]';

type CruxStep = ['tool', string, 'atom' | 'bond', number] | ['menu', 'atom' | 'bond', number, ...string[]] |
  ['label', number, string] | ['atomCustom', number, string] | ['bondCustom', number, string];
/** [the feature, the molecule drawn first, the steps that add the feature, its intended SMARTS] */
type CruxFeature = [string, string, CruxStep[], string];

const CRUX: CruxFeature[] = [
  ['toluene, its methyl bond marked aromatic (the report)', 'Cc1ccccc1', [['tool', 'bond.aromatic', 'bond', 0]],
    '[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'],
  ['aromatic bond', 'CN', [['tool', 'bond.aromatic', 'bond', 0]], '[#6]:[#7]'],
  ['any bond', 'CO', [['menu', 'bond', 0, 'bond-query', 'bond-query.any']], '[#6]~[#8]'],
  ['single or double bond', 'CO', [['menu', 'bond', 0, 'bond-query', 'bond-query.single-or-double']], '[#6]-,=[#8]'],
  ['single or aromatic bond', 'CN', [['menu', 'bond', 0, 'bond-query', 'bond-query.single-or-aromatic']], '[#6]-,:[#7]'],
  ['double or aromatic bond', 'CN', [['menu', 'bond', 0, 'bond-query', 'bond-query.double-or-aromatic']], '[#6]=,:[#7]'],
  ['ring bond', 'CC', [['menu', 'bond', 0, 'topology', 'topology.ring']], '[#6]-&@[#6]'],
  ['chain bond', 'CO', [['menu', 'bond', 0, 'topology', 'topology.chain']], '[#6]-&!@[#8]'],
  ['custom bond query', 'CC', [['bondCustom', 0, '-,=;@']], '[#6]-,=;@[#6]'],
  ['atom list', 'CN', [['label', 1, '[N,O]']], '[#6]-[#7,#8]'],
  ['NOT list', 'CN', [['label', 1, '[!N,O]']], '[#6]-[!#7&!#8]'],
  ['generic A', 'CC', [['label', 1, 'A']], '[#6]-[!#1]'],
  ['generic AH', 'CC', [['label', 1, 'AH']], '[#6]-*'],
  ['generic Q', 'CC', [['label', 1, 'Q']], '[#6]-[!#6&!#1]'],
  ['generic QH', 'CC', [['label', 1, 'QH']], '[#6]-[!#6]'],
  ['generic X', 'CC', [['label', 1, 'X']], '[#6]-[F,Cl,Br,I,At]'],
  ['generic XH', 'CC', [['label', 1, 'XH']], '[#6]-[F,Cl,Br,I,At,#1]'],
  ['generic M', 'CC', [['label', 1, 'M']], `[#6]-${RD_M}`],
  ['generic MH', 'CC', [['label', 1, 'MH']], `[#6]-${RD_M.replace('&!#1]', ']')}`],
  ['generic *', 'CC', [['label', 1, '*']], '[#6]-*'],
  ['charge', 'CN', [['menu', 'atom', 1, 'charge', 'charge.plus']], '[#6]-[#7+]'],
  ['ring bond count 2', 'CO', [['menu', 'atom', 0, 'query', 'query.ring-bond-count', 'query.ring-bond-count.2']], '[#6;x2]-[#8]'],
  ['ring bond count as drawn', 'CO', [['menu', 'atom', 0, 'query', 'query.ring-bond-count', 'query.ring-bond-count.as-drawn']], '[#6;x0]-[#8]'],
  ['substitution count 2', 'CO', [['menu', 'atom', 0, 'query', 'query.substitution', 'query.substitution.2']], '[#6;D2]-[#8]'],
  ['substitution count as drawn', 'CO', [['menu', 'atom', 0, 'query', 'query.substitution', 'query.substitution.as-drawn']], '[#6;D1]-[#8]'],
  ['unsaturated', 'CO', [['menu', 'atom', 0, 'query', 'query.unsaturated']], '[#6;$(*=,:,#*)]-[#8]'],
  ['H count 1 or more', 'CO', [['menu', 'atom', 0, 'query', 'query.h-count', 'query.h-count.1']], '[#6;h{1-}]-[#8]'],
  ['aromatic atom', 'CO', [['menu', 'atom', 0, 'query', 'query.aromaticity', 'query.aromaticity.aromatic']], '[#6;a]-[#8]'],
  ['aliphatic atom', 'CO', [['menu', 'atom', 0, 'query', 'query.aromaticity', 'query.aromaticity.aliphatic']], '[#6;A]-[#8]'],
  ['implicit H count 3', 'CO', [['menu', 'atom', 0, 'query', 'query.implicit-h', 'query.implicit-h.3']], '[#6;h3]-[#8]'],
  ['ring membership 1', 'CO', [['menu', 'atom', 0, 'query', 'query.ring-membership', 'query.ring-membership.1']], '[#6;R1]-[#8]'],
  ['ring size 6', 'CO', [['menu', 'atom', 0, 'query', 'query.ring-size', 'query.ring-size.6']], '[#6;r6]-[#8]'],
  ['connectivity 4', 'CC', [['menu', 'atom', 0, 'query', 'query.connectivity', 'query.connectivity.4']], '[#6;X4]-[#6]'],
  ['custom atom query', 'CO', [['atomCustom', 0, '[#6;$([#6]=[#8])]']], '[#6;$([#6]=[#8])]-[#8]'],
];

/** Ketcher's structure as its tools leave it, as KET (its own format, which it reads without change): two atoms 1.5
 * apart and their bond, unless the toluene. [the feature, atoms, bonds, intended SMARTS] */
type KetFeature = [string, any[], any[], string];
const at2 = (a: any, b: any, bond: any = {}): [any[], any[]] =>
  [[{location: [0, 0, 0], ...a}, {location: [1.5, 0, 0], ...b}], [{type: 1, atoms: [0, 1], ...bond}]];
const TOLUENE: [any[], any[]] = [
  [{label: 'C', location: [3.897, 1.5, 0]}, ...[[0, 1.5], [-1.299, 0.75], [-1.299, -0.75], [0, -1.5], [1.299, -0.75],
    [1.299, 0.75]].map(([x, y]) => ({label: 'C', location: [x + 1.299, y, 0]}))],
  [{type: 4, atoms: [0, 6]}, {type: 2, atoms: [1, 2]}, {type: 1, atoms: [2, 3]}, {type: 2, atoms: [3, 4]},
    {type: 1, atoms: [4, 5]}, {type: 2, atoms: [5, 6]}, {type: 1, atoms: [6, 1]}]];
const KETCHER: KetFeature[] = [
  ['toluene, its methyl bond marked aromatic (the report)', ...TOLUENE, '[#6]:[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'],
  ['aromatic bond', ...at2({label: 'C'}, {label: 'N'}, {type: 4}), '[#6]:[#7]'],
  ['any bond', ...at2({label: 'C'}, {label: 'O'}, {type: 8}), '[#6]~[#8]'],
  ['single or double bond', ...at2({label: 'C'}, {label: 'O'}, {type: 5}), '[#6]-,=[#8]'],
  ['single or aromatic bond', ...at2({label: 'C'}, {label: 'N'}, {type: 6}), '[#6]-,:[#7]'],
  ['double or aromatic bond', ...at2({label: 'C'}, {label: 'N'}, {type: 7}), '[#6]=,:[#7]'],
  ['ring bond', ...at2({label: 'C'}, {label: 'C'}, {topology: 1}), '[#6]-&@[#6]'],
  ['chain bond', ...at2({label: 'C'}, {label: 'O'}, {topology: 2}), '[#6]-&!@[#8]'],
  ['custom bond query', ...at2({label: 'C'}, {label: 'C'}, {type: undefined, customQuery: '-,=;@'}), '[#6]-,=;@[#6]'],
  ['atom list', ...at2({label: 'C'}, {type: 'atom-list', elements: ['N', 'O'], notList: false}), '[#6]-[#7,#8]'],
  ['NOT list', ...at2({label: 'C'}, {type: 'atom-list', elements: ['N', 'O'], notList: true}), '[#6]-[!#7&!#8]'],
  ['generic A', ...at2({label: 'C'}, {label: 'A'}), '[#6]-[!#1]'],
  ['generic AH', ...at2({label: 'C'}, {label: 'AH'}), '[#6]-*'],
  ['generic Q', ...at2({label: 'C'}, {label: 'Q'}), '[#6]-[!#6&!#1]'],
  ['generic QH', ...at2({label: 'C'}, {label: 'QH'}), '[#6]-[!#6]'],
  ['generic X', ...at2({label: 'C'}, {label: 'X'}), '[#6]-[F,Cl,Br,I,At]'],
  ['generic XH', ...at2({label: 'C'}, {label: 'XH'}), '[#6]-[F,Cl,Br,I,At,#1]'],
  ['generic M', ...at2({label: 'C'}, {label: 'M'}), `[#6]-${RD_M}`],
  ['generic MH', ...at2({label: 'C'}, {label: 'MH'}), `[#6]-${RD_M.replace('&!#1]', ']')}`],
  ['generic *', ...at2({label: 'C'}, {label: '*'}), '[#6]-*'],
  ['charge', ...at2({label: 'C'}, {label: 'N', charge: 1}), '[#6]-[#7+]'],
  ['ring bond count 2', ...at2({label: 'C', ringBondCount: 2}, {label: 'O'}), '[#6;x2]-[#8]'],
  ['ring bond count as drawn', ...at2({label: 'C', ringBondCount: -2}, {label: 'O'}), '[#6;x0]-[#8]'],
  ['substitution count 2', ...at2({label: 'C', substitutionCount: 2}, {label: 'O'}), '[#6;D2]-[#8]'],
  ['substitution count as drawn', ...at2({label: 'C', substitutionCount: -2}, {label: 'O'}), '[#6;D1]-[#8]'],
  ['unsaturated', ...at2({label: 'C', unsaturatedAtom: true}, {label: 'O'}), '[#6;$(*=,:,#*)]-[#8]'],
  ['H count 1 or more', ...at2({label: 'C', hCount: 2}, {label: 'O'}), '[#6;h{1-}]-[#8]'],
  ['aromatic atom', ...at2({label: 'C', queryProperties: {aromaticity: 'aromatic'}}, {label: 'O'}), '[#6;a]-[#8]'],
  ['aliphatic atom', ...at2({label: 'C', queryProperties: {aromaticity: 'aliphatic'}}, {label: 'O'}), '[#6;A]-[#8]'],
  ['implicit H count 3', ...at2({label: 'C', implicitHCount: 3}, {label: 'O'}), '[#6;h3]-[#8]'],
  ['ring membership 1', ...at2({label: 'C', queryProperties: {ringMembership: 1}}, {label: 'O'}), '[#6;R1]-[#8]'],
  ['ring size 6', ...at2({label: 'C', queryProperties: {ringSize: 6}}, {label: 'O'}), '[#6;r6]-[#8]'],
  ['connectivity 4', ...at2({label: 'C', queryProperties: {connectivity: 4}}, {label: 'C'}), '[#6;X4]-[#6]'],
  ['custom atom query', ...at2({label: 'C', queryProperties: {customQuery: '#6;$([#6]=[#8])'}}, {label: 'O'}),
    '[#6;$([#6]=[#8])]-[#8]'],
];

const ketOf = (atoms: any[], bonds: any[]) =>
  JSON.stringify({root: {nodes: [{$ref: 'mol0'}]}, mol0: {type: 'molecule', atoms, bonds}});

const DATASETS: [string, string, string][] = [
  ['mol1K', 'mol1K.csv', 'molecule'],
  ['SPGI', 'tests/spgi-100.csv', 'Structure'],
  ['smiles.csv', 'smiles.csv', 'canonical_smiles'],
];

interface CruxElement extends HTMLElement {
  readonly molfile: string;
  readonly smarts: string;
  readonly isPending: boolean;
  readonly positions: {atoms: ({x: number, y: number} | null)[], bonds: ({x: number, y: number} | null)[]};
  setValue(text: string, options?: {format?: string}): void;
}

let rdkit: any;

/** The row indices of `col` RDKit matches with `smarts`. */
function matches(col: DG.Column, smarts: string): number[] {
  const q = rdkit.get_qmol(smarts);
  const out: number[] = [];
  try {
    for (let i = 0; i < col.length; i++) {
      let m = null;
      try {
        m = rdkit.get_mol(col.get(i));
      } catch {
        m = null;
      }
      if (m === null)
        continue;
      try {
        if (m.get_substruct_match(q) !== '{}')
          out.push(i);
      } finally {
        m.delete();
      }
    }
  } finally {
    q.delete();
  }
  return out;
}

function rows(bits: DG.BitSet): number[] {
  const out: number[] = [];
  for (let i = bits.findNext(-1, true); i !== -1; i = bits.findNext(i, true))
    out.push(i);
  return out;
}

/** What the search makes of a molecule string: the SMARTS of its query (getQueryMolSafe). */
function searchReads(molecule: string): string {
  const q = getQueryMolSafe(molecule, '', rdkit);
  try {
    return q?.get_smarts() ?? '';
  } finally {
    q?.delete();
  }
}

/** `promise`, or an error saying `what` after 30 s: a step that never ends fails its feature, not the whole test. */
function within<T>(promise: Promise<T>, what: string | (() => string), ms = 30000): Promise<T> {
  return Promise.race([promise, delay(ms).then((): T => {
    throw new Error(typeof what === 'string' ? what : what());
  })]);
}

const tables = new Map<string, DG.DataFrame>();

/** A table view of `file` (read once, a clone each time), its filter panel holding the substructure card of `column` alone. */
async function openTable(file: string, column: string): Promise<{df: DG.DataFrame, filter: SubstructureFilter}> {
  grok.shell.closeAll();
  if (!tables.has(file)) {
    const read = await grok.dapi.files.readCsv(`System:AppData/Chem/${file}`);
    await grok.data.detectSemanticTypes(read);
    tables.set(file, read);
  }
  const df = tables.get(file)!.clone();
  df.name = file.split('/').pop()!.replace('.csv', '');
  const tv = grok.shell.addTableView(df);
  const fg = tv.getFiltersGroup({createDefaultFilters: false});
  fg.updateOrAdd({type: DG.FILTER_TYPE.SUBSTRUCTURE, column: column, columnName: column}, false);
  let filter: SubstructureFilter | undefined;
  await awaitCheck(() => {
    // the panel's filter is the package's own class, not this test bundle's: found by what it has
    filter = fg.filters.find((f: any) => f?.sketcher instanceof DG.chem.Sketcher && f.columnName === column) as unknown as SubstructureFilter;
    return filter !== undefined && filter.laidOut;
  }, `no substructure filter card for ${column}`, 30000);
  await awaitCheck(() => !filter!.calculating && filter!.currentSearches.size === 0, 'the pre-search did not end', 60000);
  return {df, filter: filter!};
}

/** Opens the card's sketcher dialog with `type` the session's sketcher, as a click on its thumbnail or Sketch link does. */
async function openDialog(filter: SubstructureFilter, type: string): Promise<DG.chem.SketcherBase> {
  grok.chem.currentSketcherType = type;
  const s = filter.sketcher;
  (s.extSketcherDiv.querySelector('canvas, .sketch-link') as HTMLElement ?? s.extSketcherDiv).click();
  await awaitCheck(() => s.sketcherDialogOpened, 'the sketcher dialog did not open', 5000);
  const impl = await within(s.sketcherReady(), `${type} did not get ready in the dialog`);
  return impl;
}

function dialogOf(filter: SubstructureFilter): HTMLElement {
  return Array.from(document.querySelectorAll<HTMLElement>('.d4-dialog')).find((d) => d.contains(filter.sketcher.host))!;
}

function pressDialogButton(filter: SubstructureFilter, name: 'OK' | 'CANCEL'): void {
  const button = Array.from(dialogOf(filter).querySelectorAll<HTMLElement>('button, .ui-btn'))
    .find((b) => b.textContent!.trim().toUpperCase() === name);
  button!.click();
}

async function settle(filter: SubstructureFilter): Promise<void> {
  await delay(50);
  await awaitCheck(() => !filter.calculating && filter.currentSearches.size === 0,
    `the search did not end (searches ${filter.currentSearches.size}, calculating ${filter.calculating})`, 20000);
}

// ---------------------------------------------------------------- Crux, by its own controls

function cruxOf(filter: SubstructureFilter): CruxElement {
  return filter.sketcher.sketcher!.root.querySelector<CruxElement>('crux-sketch')!;
}

function shadow(el: CruxElement, testid: string): HTMLElement {
  const found = el.shadowRoot!.querySelector<HTMLElement>(`[data-testid="${testid}"]`);
  if (!found)
    throw new Error(`Crux shows no ${testid}`);
  return found;
}

function pointAt(el: CruxElement, kind: 'atom' | 'bond', i: number): {clientX: number, clientY: number} {
  const p = (kind === 'atom' ? el.positions.atoms : el.positions.bonds)[i];
  if (!p)
    throw new Error(`Crux draws no ${kind} ${i}`);
  const r = el.getBoundingClientRect();
  return {clientX: r.left + p.x, clientY: r.top + p.y};
}

/** A user's click (left) on the canvas at a point, as the pointer's events reach it. */
function click(el: CruxElement, at: {clientX: number, clientY: number}): void {
  const canvas = shadow(el, 'canvas');
  const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true, ...at,
    pointerId: 1, pointerType: 'mouse', isPrimary: true, button: 0, buttons});
  canvas.dispatchEvent(new PointerEvent('pointerdown', init(1)));
  canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
}

/** A user's right-click on the canvas at a point: the canvas's own menu opens there. */
function rightClick(el: CruxElement, at: {clientX: number, clientY: number}): void {
  const canvas = shadow(el, 'canvas');
  const init = (buttons: number): PointerEventInit => ({bubbles: true, cancelable: true, composed: true, ...at,
    pointerId: 1, pointerType: 'mouse', isPrimary: true, button: 2, buttons});
  canvas.dispatchEvent(new PointerEvent('pointerdown', init(2)));
  canvas.dispatchEvent(new PointerEvent('pointerup', init(0)));
  canvas.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, composed: true, ...at, button: 2}));
}

async function cruxSteps(el: CruxElement, steps: CruxStep[]): Promise<void> {
  for (const step of steps) {
    if (step[0] === 'tool') {
      shadow(el, `toolbar.${step[1]}`).click();
      click(el, pointAt(el, step[2], step[3]));
    } else if (step[0] === 'menu') {
      const [, kind, i, ...items] = step;
      rightClick(el, pointAt(el, kind, i));
      for (const item of items) {
        await awaitCheck(() => el.shadowRoot!.querySelector(`[data-testid="menu.context.${item}"]`) !== null,
          `Crux's canvas menu shows no ${item}`, 3000);
        shadow(el, `menu.context.${item}`).click();
      }
    } else if (step[0] === 'label') {
      rightClick(el, pointAt(el, 'atom', step[1]));
      shadow(el, 'menu.context.edit-label').click();
      await awaitCheck(() => el.shadowRoot!.querySelector('[data-testid="canvas.label-editor"]') !== null,
        'Crux opened no label editor', 3000);
      const input = shadow(el, 'canvas.label-editor') as HTMLInputElement;
      input.value = step[2];
      input.dispatchEvent(new InputEvent('input', {bubbles: true, composed: true}));
      input.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true, cancelable: true, composed: true}));
    } else {
      const atom = step[0] === 'atomCustom';
      rightClick(el, pointAt(el, atom ? 'atom' : 'bond', step[1]));
      shadow(el, `menu.context.${atom ? 'atom-properties' : 'bond-properties'}`).click();
      const d = atom ? 'dialog.atom-properties' : 'dialog.bond-properties';
      await awaitCheck(() => el.shadowRoot!.querySelector(`[data-testid="${d}"]`) !== null, `Crux opened no ${d}`, 3000);
      const custom = shadow(el, atom ? `${d}.query.custom` : `${d}.custom`) as HTMLInputElement;
      custom.click();
      const smarts = shadow(el, atom ? `${d}.query.smarts` : `${d}.smarts`) as HTMLInputElement;
      smarts.value = step[2];
      smarts.dispatchEvent(new InputEvent('input', {bubbles: true, composed: true}));
      shadow(el, `${d}.apply`).click();
    }
    await awaitCheck(() => !el.isPending, 'Crux did not settle', 3000);
  }
}

/** The query the dialog's sketcher itself holds, not the host's getters (which give back the string the host set until the
 * user's next edit): its own molblock, as the search reads it, and its own SMARTS. */
async function heldQuery(filter: SubstructureFilter): Promise<{search: string, smarts: string}> {
  const impl: any = filter.sketcher.sketcher;
  if (impl.root.querySelector('crux-sketch')) {
    const el = cruxOf(filter);
    return {search: searchReads(el.molfile), smarts: el.smarts};
  }
  if (impl._sketcher?.getMolfile) {
    // Ketcher: the adapter's own export of what Ketcher shows, from its last change (the host's getters give the string
    // the host set), and Ketcher's SMARTS, which the adapter asks for in turn with its other conversions (Ketcher's
    // conversion service drops a request made beside another one unanswered)
    // what the host set loads asynchronously, and its change exports what Ketcher shows (before it, the adapter holds the
    // host's string)
    await within(impl._loading ?? Promise.resolve(), 'Ketcher did not load what the host set');
    await delay(300);
    await awaitCheck(() => impl._molV2000 !== null, 'Ketcher exported nothing', 10000);
    return {search: searchReads(impl._molV2000), smarts: await impl.getSmarts()};
  }
  return {search: searchReads(impl._sketcher.getMolFile()), smarts: ''};
}

/** The user's flow on one feature: the card's dialog opened with `sketcher`, the molecule drawn first, then the feature,
 * OK; the rows, the sketcher reopened, a switch to `other` and back. The problems found, one line each. */
async function roundTrip(file: string, column: string, sketcher: string, other: string, name: string,
  draw: (filter: SubstructureFilter) => Promise<void>, intended: string, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const want = matches(df.col(column)!, intended);
  try {
    at.stage = 'opening the dialog';
    await openDialog(filter, sketcher);
    at.stage = 'drawing';
    await draw(filter);
    at.stage = 'the search of the drawing';
    await settle(filter);
    at.stage = 'reading the drawing';
    const drawn = await heldQuery(filter);
    pressDialogButton(filter, 'OK');
    at.stage = 'the search after OK';
    await settle(filter);
    const passing = rows(df.filter);
    if (JSON.stringify(passing) !== JSON.stringify(want))
      problems.push(`${name}: ${passing.length} rows pass, ${intended} matches ${want.length}`);
    if (searchReads(filter.currentMolecule) !== drawn.search)
      problems.push(`${name}: the filter keeps ${searchReads(filter.currentMolecule)}, drawn ${drawn.search}`);
    at.stage = 'reopening';
    await openDialog(filter, sketcher);
    at.stage = 'reading the reopened sketcher';
    const reopened = await heldQuery(filter);
    if (reopened.search !== drawn.search || reopened.smarts !== drawn.smarts)
      problems.push(`${name}: reopened, ${JSON.stringify(reopened)}; drawn, ${JSON.stringify(drawn)}`);
    for (const type of [other, sketcher]) {
      at.stage = `switching to ${type}`;
      grok.chem.currentSketcherType = type;
      filter.sketcher.sketcherType = type;
      await within(filter.sketcher.sketcherReady(), `${type} did not get ready after the switch`);
      // Ketcher's drawing loads asynchronously after it is ready
      await delay(type === 'Ketcher' ? 1000 : 0);
    }
    at.stage = 'reading the sketcher switched back';
    const back = await heldQuery(filter);
    if (back.search !== drawn.search || back.smarts !== drawn.smarts)
      problems.push(`${name}: through ${other} and back, ${JSON.stringify(back)}; drawn, ${JSON.stringify(drawn)}`);
    pressDialogButton(filter, 'CANCEL');
  } catch (e: any) {
    problems.push(`${name}: ${e?.message ?? e}`);
  }
  return problems;
}

/** Draws on Crux in the filter's dialog: the molecule, then the feature by Crux's own controls. */
const cruxDraws = (base: string, steps: CruxStep[]) => async (filter: SubstructureFilter): Promise<void> => {
  const el = cruxOf(filter);
  el.setValue(base);
  await settle(filter);
  await cruxSteps(el, steps);
};

/** Draws on Ketcher in the filter's dialog: the molecule (each query atom an element, each query bond single), then the
 * structure with the feature, as its tools leave it. */
const ketcherDraws = (atoms: any[], bonds: any[]) => async (filter: SubstructureFilter): Promise<void> => {
  const ketcher = (filter.sketcher.sketcher as any)._sketcher;
  const plainAtoms = atoms.map((a) => ({label: a.type === 'atom-list' || !/^[A-Z][a-z]?$/.test(a.label) ||
    ['A', 'Q', 'X', 'M'].includes(a.label) ? 'C' : a.label, location: a.location}));
  const plainBonds = bonds.map((b) => ({type: (b.type ?? 1) > 3 ? 1 : (b.type ?? 1), atoms: b.atoms}));
  const load = (ket: string) => Promise.race([ketcher.setMolecule(ket),
    delay(15000).then(() => {throw new Error('Ketcher did not load the structure');})]);
  await load(ketOf(plainAtoms, plainBonds));
  await settle(filter);
  await load(ketOf(atoms, bonds));
  // the change's exports run in turn after it (V2000 at once, Indigo's then)
  await delay(1000);
};

/** How long a test of the category may take before it starts no other feature: under its own timeout, so that it ends
 * by itself, saying what it did not reach, and leaves no feature running into the tests after it. */
const TEST_BUDGET_MS = 15 * 60000;
const FEATURE_MS = 90000;

/** Each feature's round trip in turn, each given FEATURE_MS; the problems found. Two features over their time end the
 * test (a sketcher that stopped answering), as does the test's budget; either way the problems name the features not
 * reached, and a failing test names each feature over 30 s and the page's JS heap, so that a slow run says where its
 * time went. */
async function eachFeature(features: [string, (at: {stage: string}) => Promise<string[]>][]): Promise<string[]> {
  const problems: string[] = [];
  const slow: string[] = [];
  const heap = (): number => Math.round(((performance as any).memory?.usedJSHeapSize ?? 0) / 1e6);
  const heapBefore = heap();
  const start = Date.now();
  let i = 0;
  try {
    for (; i < features.length; i++) {
      if (Date.now() - start > TEST_BUDGET_MS - FEATURE_MS ||
        problems.filter((p) => p.includes(`not done in ${FEATURE_MS / 1000} s`)).length >= 2)
        break;
      const [name, run] = features[i];
      const at = {stage: ''};
      const t0 = Date.now();
      problems.push(...await within(run(at), () => `${name}: not done in ${FEATURE_MS / 1000} s, at ${at.stage}`,
        FEATURE_MS).catch((e) => [String(e?.message ?? e)]));
      const ms = Date.now() - t0;
      if (ms > 30000)
        slow.push(`${name} ${Math.round(ms / 1000)} s`);
    }
  } finally {
    for (const d of (DG.Dialog as any).getOpenDialogs?.() ?? [])
      d.close();
    grok.shell.closeAll();
  }
  if (i < features.length) {
    problems.push(`not reached in ${Math.round((Date.now() - start) / 1000)} s: ` +
      features.slice(i).map(([n]) => n).join(', '));
  }
  if (problems.length) {
    problems.push(`features over 30 s: ${slow.join('; ') || 'none'}; ` +
      `JS heap ${heapBefore} MB before, ${heap()} MB after`);
  }
  return problems;
}

category('query round trip', () => {
  let sketcherType: string;

  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
    rdkit = chemCommonRdKit.getRdKitModule();
    sketcherType = grok.chem.currentSketcherType;
  });

  after(async () => {
    grok.chem.currentSketcherType = sketcherType;
    grok.shell.closeAll();
  });

  for (const [ds, file, column] of DATASETS) {
    const cruxTest = `${ds}: every query feature drawn in Crux filters its rows, ` +
      'and comes back reopened and through Ketcher';
    test(cruxTest, async () => {
      const problems = await eachFeature(CRUX.map(([name, base, steps, intended]) => [name, (at: {stage: string}) =>
        roundTrip(file, column, 'Crux', 'Ketcher', name, cruxDraws(base, steps), intended, at)]));
      expect(problems.length, 0, problems.join('\n'));
    }, {timeout: 1200000});

    const ketcherTest = `${ds}: every query feature drawn in Ketcher filters its rows, ` +
      'and comes back reopened and through Crux';
    test(ketcherTest, async () => {
      const problems = await eachFeature(KETCHER.map(([name, atoms, bonds, intended]) => [name, (at: {stage: string}) =>
        roundTrip(file, column, 'Ketcher', 'Crux', name, ketcherDraws(atoms, bonds), intended, at)]));
      expect(problems.length, 0, problems.join('\n'));
    }, {timeout: 1200000});
  }
});
