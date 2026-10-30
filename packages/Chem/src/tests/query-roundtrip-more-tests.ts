/* eslint-disable max-len */
/* More of the substructure filter's query round trips (the product owner, 2026-10-13 evening: "go through some datagrok
 * scenarios. check more back and forth tests for complex queries, query structures, or even non realistic drawings and
 * back and forth in datagrok queries"). The categories below are "query back and forth: ..." (one filter, a prefix, runs
 * them all; not under "query round trip", whose page they would grow past what this machine holds), each case a test of
 * its own:
 * - complex queries: several query features in one drawing, by each sketcher's own means (Crux's controls, Ketcher's
 *   structure as its tools leave it, OpenChemLib's editor given a query molfile);
 * - entered as text: SMARTS, query molfiles and CXSMILES pasted (or typed and entered) in the sketcher's field;
 * - non-realistic drawings: what no chemist means, drawn anyway (pentavalent carbon, a lone R label, a ring of 20...);
 * - routes: reopened and OK unchanged, switched through every sketcher, edited after a reopening, kept through a layout,
 *   cleared; and the molecule cell editor's round trip.
 * Every case is judged by RDKit, never by eye: the rows the filter passes are RDKit's matches of the query intended; the
 * sketcher's own drawing, as the search reads it, matches the same rows (it shows the query the filter searches: not
 * empty, not changed); reopened from the card, switched to another sketcher and back, and reopened after a Cancel, the
 * sketcher holds the query it was given. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {after, waitUntil, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {PackageFunctions} from '../package';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {SubstructureFilter} from '../widgets/chem-substructure-filter';
import {CruxElement, CruxStep, DATASETS, cruxOf, cruxSteps, heldQuery, ketOf, ketcherExported, notationBalloons,
  openDialog, openTable, pasteInField, pressDialogButton, rows, searchReads, settle, shadow, tables, useRdKit,
  within} from './query-roundtrip-tests';

type SketcherType = 'Crux' | 'Ketcher' | 'OpenChemLib';
type Ds = 'mol1K' | 'SPGI' | 'smiles.csv';

/** How a case draws its query in the dialog's sketcher. */
type Draw =
  /** Crux: `crux` read as Crux's own value (as a paste into its canvas; in `format` when given), then `steps` by its
   * controls */
  {crux: string, format?: string, steps?: CruxStep[]} |
  /** Crux: `cruxSteps` by its controls, on what it shows */
  {cruxSteps: CruxStep[]} |
  /** Ketcher: its structure as its tools leave it (KET) */
  {ket: [any[], any[]]} |
  /** Ketcher: a molfile Ketcher loads as its own */
  {ketcher: string} |
  /** OpenChemLib: a molfile its editor reads as its own */
  {ocl: string} |
  /** the host's field: `field` pasted (the default), or typed and entered */
  {field: string, by?: 'paste' | 'enter'};

/** The query intended: a SMARTS (RDKit's `get_qmol`), or a molecule string RDKit reads as a molecule (`get_mol`), its
 * matches those of the molecule as a query (a SMILES entered: the search reads a plain SMILES as a molecule). */
type Intended = string | {mol: string};

interface Case {
  name: string;
  ds: Ds;
  sketcher: SketcherType;
  draw: Draw;
  intended: Intended;
  /** the sketchers switched to in turn once the dialog is reopened, and back (by default Ketcher from Crux, else Crux) */
  through?: SketcherType[];
  /** the sketchers passed through that hold the query too (each sketcher holds it once switched back to) */
  checkThrough?: SketcherType[];
  skip?: string;
}

let rdkit: any;
const FILE: {[ds in Ds]: [string, string]} =
  Object.fromEntries(DATASETS.map(([ds, file, column]) => [ds, [file, column]])) as any;

// ---------------------------------------------------------------- RDKit's judgement

/** Each dataset's molecules read by RDKit once (a test's table is a clone of the one read, its rows in that order). */
const datasetMols = new Map<string, any[]>();

function molsOf(file: string, col: DG.Column): any[] {
  if (!datasetMols.has(file)) {
    const mols: any[] = [];
    for (let i = 0; i < col.length; i++) {
      let m = null;
      try {
        m = rdkit.get_mol(col.get(i));
      } catch {
        m = null;
      }
      mols.push(m);
    }
    datasetMols.set(file, mols);
  }
  return datasetMols.get(file)!;
}

function rowsMatching(mols: any[], q: any): number[] {
  const out: number[] = [];
  for (let i = 0; i < mols.length; i++) {
    if (mols[i] && mols[i].get_substruct_match(q) !== '{}')
      out.push(i);
  }
  return out;
}

function intendedRows(mols: any[], intended: Intended): number[] {
  const q = typeof intended === 'string' ? rdkit.get_qmol(intended) : rdkit.get_mol(intended.mol);
  if (!q)
    throw new Error(`RDKit cannot read the query intended, ${JSON.stringify(intended)}`);
  try {
    return rowsMatching(mols, q);
  } finally {
    q.delete();
  }
}

/** The rows a molecule string matches as the search reads it (getQueryMolSafe); every row for an empty one (the filter
 * then filters nothing); null when the search cannot read it. */
function drawingRows(mols: any[], molString: string): number[] | null {
  if (!molString || DG.chem.Sketcher.isEmptyMolfile(molString))
    return mols.map((_, i) => i);
  const q = getQueryMolSafe(molString, '', rdkit);
  if (!q)
    return null;
  try {
    return rowsMatching(mols, q);
  } finally {
    q.delete();
  }
}

const describe = (intended: Intended): string => typeof intended === 'string' ? intended : intended.mol;
const sameRows = (a: number[] | null, b: number[]): boolean => a !== null && a.length === b.length && a.every((x, i) => x === b[i]);

// ---------------------------------------------------------------- what a sketcher holds

/** The drawing the dialog's sketcher itself shows (not the host's string): its molfile, as the search reads it, and the
 * sketcher's own SMARTS (Crux's, Ketcher's; none for OpenChemLib). */
interface Held {molfile: string, search: string, smarts: string}

function isCrux(impl: any): boolean {
  return !!impl?.root?.querySelector('crux-sketch');
}

function isKetcher(impl: any): boolean {
  return !!impl?._sketcher?.getMolfile;
}

async function held(filter: SubstructureFilter): Promise<Held> {
  const impl: any = filter.sketcher.sketcher;
  if (isCrux(impl)) {
    const el = cruxOf(filter);
    await waitUntil(() => !el.isPending, 'Crux did not settle', 5000);
    const molfile = (el as any).isEmpty ? '' : el.molfile;
    return {molfile, search: searchReads(molfile), smarts: el.smarts};
  }
  if (isKetcher(impl)) {
    // an empty Ketcher given nothing to load exports nothing
    if (impl._molV2000 === null && impl._sketcher.editor.struct().atoms.size === 0 && impl.importedMoleculesCounter === 0)
      return {molfile: '', search: '', smarts: ''};
    const h = await heldQuery(filter);
    return {molfile: impl._molV2000 ?? '', search: h.search, smarts: h.smarts};
  }
  const molfile: string = impl._sketcher.getMolFile();
  return {molfile, search: searchReads(molfile), smarts: ''};
}

const heldText = (h: Held): string => JSON.stringify({search: h.search, smarts: h.smarts});

/** Ketcher's own atoms, their labels and query properties (a failure's diagnosis); '' for another sketcher. */
function ketcherAtoms(filter: SubstructureFilter): string {
  const impl: any = filter.sketcher.sketcher;
  if (!isKetcher(impl))
    return '';
  try {
    return `; Ketcher's atoms: ${JSON.stringify([...impl._sketcher.editor.struct().atoms.values()].map((a: any) =>
      [a.label, a.atomList?.ids ?? null, a.queryProperties ?? null]))}`;
  } catch (e: any) {
    return `; Ketcher's atoms unread: ${e?.message ?? e}`;
  }
}

// ---------------------------------------------------------------- drawing

/** Ketcher's structure as its tools leave it (KET) for `smiles` laid out by RDKit (Kekulé, as Ketcher's templates draw
 * rings), each atom and bond given the properties listed by its index (an atom list: `{type: 'atom-list', elements}`). */
function ketFrom(smiles: string, atomProps: {[i: number]: any} = {}, bondProps: {[i: number]: any} = {}): [any[], any[]] {
  const mol = rdkit.get_mol(smiles);
  let lines: string[];
  try {
    mol.set_new_coords();
    lines = mol.get_molblock().split('\n');
  } finally {
    mol.delete();
  }
  const n = parseInt(lines[3].substring(0, 3));
  const nb = parseInt(lines[3].substring(3, 6));
  const atoms: any[] = [];
  for (let i = 0; i < n; i++) {
    const l = lines[4 + i];
    const location = [parseFloat(l.substring(0, 10)), parseFloat(l.substring(10, 20)), 0];
    const p = atomProps[i] ?? {};
    atoms.push(p.type === 'atom-list' ? {location, notList: false, ...p} :
      {label: l.substring(31, 34).trim(), location, ...p});
  }
  const bonds: any[] = [];
  for (let j = 0; j < nb; j++) {
    const l = lines[4 + n + j];
    bonds.push({type: parseInt(l.substring(6, 9)), atoms: [parseInt(l.substring(0, 3)) - 1, parseInt(l.substring(3, 6)) - 1],
      ...(bondProps[j] ?? {})});
  }
  return [atoms, bonds];
}

/** Ketcher has loaded what it was last given and exported it. */
async function ketcherSettled(impl: any): Promise<void> {
  await within(impl._loading ?? Promise.resolve(), 'Ketcher did not load what it was given');
  await ketcherExported(impl, 3000);
}

async function draw(filter: SubstructureFilter, d: Draw): Promise<void> {
  const impl: any = filter.sketcher.sketcher;
  if ('crux' in d) {
    const el = cruxOf(filter);
    el.setValue(d.crux, d.format ? {format: d.format} : undefined);
    await settle(filter);
    const error = (el as any).lastError;
    if (error)
      throw new Error(`Crux refused ${d.crux}: ${error.message}`);
    await cruxSteps(el, d.steps ?? []);
  } else if ('cruxSteps' in d)
    await cruxSteps(cruxOf(filter), d.cruxSteps);
  else if ('ket' in d) {
    const [atoms, bonds] = d.ket;
    const plainAtoms = atoms.map((a) => ({label: a.type === 'atom-list' || !/^[A-Z][a-z]?$/.test(a.label) ||
      ['A', 'Q', 'X', 'M'].includes(a.label) ? 'C' : a.label, location: a.location}));
    const plainBonds = bonds.map((b) => ({type: (b.type ?? 1) > 3 ? 1 : (b.type ?? 1), atoms: b.atoms}));
    const load = (ket: string) => within(impl._sketcher.setMolecule(ket), 'Ketcher did not load the structure', 15000);
    await load(ketOf(plainAtoms, plainBonds));
    await settle(filter);
    const exported = impl._exportId;
    await load(ketOf(atoms, bonds));
    await ketcherExported(impl, 3000, exported);
  } else if ('ketcher' in d) {
    await ketcherSettled(impl);
    const exported = impl._exportId;
    await within(impl._sketcher.setMolecule(d.ketcher), 'Ketcher did not load the structure', 15000);
    await ketcherExported(impl, 3000, exported);
  } else if ('ocl' in d)
    impl._sketcher.setMolFile(d.ocl);
  else {
    if (d.by === 'enter') {
      filter.sketcher.molInput.value = d.field;
      filter.sketcher.molInput.dispatchEvent(new KeyboardEvent('keydown', {key: 'Enter', bubbles: true, cancelable: true}));
    } else
      pasteInField(filter, d.field);
    if (isKetcher(filter.sketcher.sketcher))
      await ketcherSettled(filter.sketcher.sketcher);
  }
  await settle(filter);
}

/** Switches the open dialog's sketcher to `type`, as its ≡ menu does, and waits until it shows the host's molecule. */
async function switchTo(filter: SubstructureFilter, type: SketcherType): Promise<void> {
  grok.chem.currentSketcherType = type;
  filter.sketcher.sketcherType = type;
  await within(filter.sketcher.sketcherReady(), `${type} did not get ready after the switch`);
  if (type === 'Ketcher')
    await ketcherSettled(filter.sketcher.sketcher);
}

function closeEverything(): void {
  for (const d of (DG.Dialog as any).getOpenDialogs?.() ?? [])
    d.close();
  DG.Balloon.closeAll();
  grok.shell.closeAll();
}

// ---------------------------------------------------------------- the round trip

/** The user's flow on one case: the card's dialog opened with the case's sketcher, the query drawn, OK; the rows; the
 * sketcher reopened; switched through the others and back; Cancel (the rows stay); reopened once more. The problems
 * found, one line each. */
async function queryCase(c: Case, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE[c.ds];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const mols = molsOf(file, df.col(column)!);
  const want = intendedRows(mols, c.intended);
  const what = `${describe(c.intended)} matches ${want.length}`;
  const checkDrawing = (stage: string, h: Held): void => {
    const r = drawingRows(mols, h.molfile);
    if (!sameRows(r, want)) {
      problems.push(`${stage}: the sketcher's drawing ${!h.molfile ? 'is empty' : r === null ?
        `is not read by the search (${h.molfile.length} chars)` : `matches ${r.length} rows (${h.search})`}; ${what}` +
        // its lines joined by ' / ': the runner shows a message from its last line on
        (stage === 'drawn' ? `; its molfile: ${h.molfile.split(/\r?\n/).join(' / ')}${ketcherAtoms(filter)}` : ''));
    }
  };
  at.stage = 'opening the dialog';
  await openDialog(filter, c.sketcher);
  at.stage = 'drawing';
  await draw(filter, c.draw);
  at.stage = 'reading the drawing';
  const drawn = await held(filter);
  checkDrawing('drawn', drawn);
  pressDialogButton(filter, 'OK');
  at.stage = 'the search after OK';
  await settle(filter);
  const passing = rows(df.filter);
  if (!sameRows(passing, want))
    problems.push(`after OK: ${passing.length} rows pass (the filter's ${JSON.stringify(searchReads(filter.currentMolecule))}); ${what}`);
  DG.Balloon.closeAll();
  at.stage = 'reopening';
  await openDialog(filter, c.sketcher);
  if (isKetcher(filter.sketcher.sketcher))
    await ketcherSettled(filter.sketcher.sketcher);
  at.stage = 'reading the reopened sketcher';
  const reopened = await held(filter);
  if (heldText(reopened) !== heldText(drawn))
    problems.push(`reopened: ${heldText(reopened)}; drawn: ${heldText(drawn)}`);
  const warned = notationBalloons();
  if (warned)
    problems.push(`reopened with a warning: ${warned}`);
  const through = c.through ?? (c.sketcher === 'Crux' ? ['Ketcher' as SketcherType] : ['Crux' as SketcherType]);
  for (const type of through) {
    at.stage = `switching to ${type}`;
    await switchTo(filter, type);
    if (c.checkThrough?.includes(type))
      checkDrawing(`in ${type}`, await held(filter));
  }
  at.stage = `switching back to ${c.sketcher}`;
  await switchTo(filter, c.sketcher);
  const back = await held(filter);
  if (heldText(back) !== heldText(drawn))
    problems.push(`through ${through.join(', ')} and back: ${heldText(back)}; drawn: ${heldText(drawn)}`);
  pressDialogButton(filter, 'CANCEL');
  at.stage = 'the search after Cancel';
  await settle(filter);
  const afterCancel = rows(df.filter);
  if (!sameRows(afterCancel, want))
    problems.push(`after Cancel: ${afterCancel.length} rows pass; ${what}`);
  at.stage = 'reopening after Cancel';
  await openDialog(filter, c.sketcher);
  if (isKetcher(filter.sketcher.sketcher))
    await ketcherSettled(filter.sketcher.sketcher);
  const again = await held(filter);
  if (heldText(again) !== heldText(drawn))
    problems.push(`reopened after Cancel: ${heldText(again)}; drawn: ${heldText(drawn)}`);
  pressDialogButton(filter, 'CANCEL');
  await settle(filter);
  return problems;
}

/** Runs `flow` in a minute and a half at most, and fails the test with the problems it found. */
async function expectNoProblems(flow: (at: {stage: string}) => Promise<string[]>): Promise<void> {
  const at = {stage: ''};
  let problems: string[];
  try {
    problems = await within(flow(at), () => `not done in 90 s, at ${at.stage}`, 90000);
  } catch (e: any) {
    problems = [`${at.stage}: ${e?.message ?? e}`];
  } finally {
    closeEverything();
  }
  // on one line: the runner shows a failure's message by its last line
  expect(problems.length, 0, problems.join(' || '));
}

function caseTests(cases: Case[]): void {
  for (const c of cases) {
    test(`${c.ds}: ${c.sketcher}: ${c.name}`, async () => {
      await expectNoProblems((at) => queryCase(c, at));
    }, {timeout: 120000, skipReason: c.skip});
  }
}

const roundTripBefore = async (): Promise<string> => {
  rdkit = await useRdKit();
  return grok.chem.currentSketcherType;
};

function roundTripAfter(sketcherType: string): void {
  grok.chem.currentSketcherType = sketcherType;
  closeEverything();
}

/** Frees the molecules read for the datasets (the last category's `after`). */
function freeDatasetMols(): void {
  for (const mols of datasetMols.values()) {
    for (const m of mols)
      m?.delete();
  }
  datasetMols.clear();
}

// ---------------------------------------------------------------- the molfiles the cases enter (RDKit's layouts)

/** Chlorobenzene, its chlorine the list [Cl,Br,I] (`M  ALS`), its ipso carbon two ring bonds (`M  RBC`): V2000. */
const V2000_ALS_RBC = [
  '', '  list-rbc', '',
  '  7  7  0  0  0  0  0  0  0  0999 V2000',
  '    2.5714    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.0714    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.3214   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.1786   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.9286    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.1786    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.3214    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0', '  2  3  2  0', '  3  4  1  0', '  4  5  2  0', '  5  6  1  0', '  6  7  2  0', '  7  2  1  0',
  'M  ALS   1  3 F Cl  Br  I   ', 'M  RBC  1   2   2', 'M  END', '',
].join('\n');
const V2000_ALS_RBC_SMARTS = '[#17,#35,#53]-[#6;x2]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';

/** Quinoline, its nitrogen the NOT list [!C,S] (`M  ALS ... T`): V2000, Kekulé. */
const V2000_NOT_QUINOLINE = [
  '', '  not-list', '',
  ' 10 11  0  0  0  0  0  0  0  0999 V2000',
  '   -2.5981   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -2.5981    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.0000    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990    1.5000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0',
  '    2.5981    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    2.5981   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.0000   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  2  3  1  0', '  3  4  2  0', '  4  5  1  0', '  5  6  2  0', '  6  7  1  0', '  7  8  2  0',
  '  8  9  1  0', '  9 10  2  0', ' 10  1  1  0', '  9  4  1  0',
  'M  ALS   5  2 T C   S   ', 'M  END', '',
].join('\n');
const NOT_QUINOLINE_SMARTS = '[#6]1:[#6]:[#6]:[#6]2:[!#6&!#16]:[#6]:[#6]:[#6]:[#6]:2:[#6]:1';

/** Cyclohexane, its six bonds "any" (type 8): V2000. */
const V2000_ANY_RING = [
  '', '  any-ring', '',
  '  6  6  0  0  0  0  0  0  0  0999 V2000',
  '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  8  0', '  2  3  8  0', '  3  4  8  0', '  4  5  8  0', '  5  6  8  0', '  6  1  8  0',
  'M  END', '',
].join('\n');
const ANY_RING_SMARTS = '[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]~1';

/** SPGI's triazole, Kekulé, with four V3000 query features: an atom list ([C,N]), a substitution count (SUBST=3), a ring
 * bond count (RBCNT=-1: none) and an unsaturated atom (UNSAT=1). */
const V3000_TRIAZOLE = [
  '', '  v3000-queries', '',
  '  0  0  0  0  0  0  0  0  0  0999 V3000',
  'M  V30 BEGIN CTAB', 'M  V30 COUNTS 8 8 0 0 0', 'M  V30 BEGIN ATOM',
  'M  V30 1 N -0.750000 -1.593739 0.000000 0', 'M  V30 2 N 0.750000 -1.593739 0.000000 0',
  'M  V30 3 C 1.213525 -0.167155 0.000000 0 SUBST=3', 'M  V30 4 [C,N] 2.640110 0.296371 0.000000 0',
  'M  V30 5 N 0.000000 0.714523 0.000000 0', 'M  V30 6 C 0.000000 2.214523 0.000000 0 RBCNT=-1',
  'M  V30 7 C -1.213525 -0.167155 0.000000 0', 'M  V30 8 C -2.640110 0.296371 0.000000 0 UNSAT=1',
  'M  V30 END ATOM', 'M  V30 BEGIN BOND',
  'M  V30 1 1 1 2', 'M  V30 2 2 2 3', 'M  V30 3 1 3 4', 'M  V30 4 1 3 5', 'M  V30 5 1 5 6', 'M  V30 6 1 5 7',
  'M  V30 7 1 7 8', 'M  V30 8 2 7 1',
  'M  V30 END BOND', 'M  V30 END CTAB', 'M  END', '',
].join('\n');
const V3000_TRIAZOLE_SMARTS = '[#7]1:[#7]:[#6;D3](-[#6,#7]):[#7](-[#6;x0]):[#6]:1-[#6;$(*=,:,#*)]';

/** Quinoline, its nitrogen the V3000 NOT list "NOT [C,S]". */
const V3000_NOT_QUINOLINE = [
  '', '  v3000-not-list', '',
  '  0  0  0  0  0  0  0  0  0  0999 V3000',
  'M  V30 BEGIN CTAB', 'M  V30 COUNTS 10 11 0 0 0', 'M  V30 BEGIN ATOM',
  'M  V30 1 C -2.598076 -0.750000 0.000000 0', 'M  V30 2 C -2.598076 0.750000 0.000000 0',
  'M  V30 3 C -1.299038 1.500000 0.000000 0', 'M  V30 4 C -0.000000 0.750000 0.000000 0',
  'M  V30 5 "NOT [C,S]" 1.299038 1.500000 0.000000 0', 'M  V30 6 C 2.598076 0.750000 0.000000 0',
  'M  V30 7 C 2.598076 -0.750000 0.000000 0', 'M  V30 8 C 1.299038 -1.500000 0.000000 0',
  'M  V30 9 C -0.000000 -0.750000 0.000000 0', 'M  V30 10 C -1.299038 -1.500000 0.000000 0',
  'M  V30 END ATOM', 'M  V30 BEGIN BOND',
  'M  V30 1 2 1 2', 'M  V30 2 1 2 3', 'M  V30 3 2 3 4', 'M  V30 4 1 4 5', 'M  V30 5 2 5 6', 'M  V30 6 1 6 7',
  'M  V30 7 2 7 8', 'M  V30 8 1 8 9', 'M  V30 9 2 9 10', 'M  V30 10 1 10 1', 'M  V30 11 1 9 4',
  'M  V30 END BOND', 'M  V30 END CTAB', 'M  END', '',
].join('\n');

/** Benzene drawn Kekulé, three of its bonds then marked aromatic. */
const PARTLY_AROMATIC = [
  '', '  partly-aromatic', '',
  '  6  6  0  0  0  0  0  0  0  0999 V2000',
  '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  4  0', '  2  3  4  0', '  3  4  4  0', '  4  5  1  0', '  5  6  2  0', '  6  1  1  0',
  'M  END', '',
].join('\n');

/** Three carbons in a chain joined by aromatic bonds (no ring). */
const CHAIN_AROMATIC = [
  '', '  chain-aromatic', '',
  '  3  2  0  0  0  0  0  0  0  0999 V2000',
  '   -1.2990   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990   -0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  4  0', '  2  3  4  0',
  'M  END', '',
].join('\n');

/** Ethane, and an oxygen drawn on its second carbon (the same coordinates), bonded to nothing. */
const OVERLAPPING = [
  '', '  overlapping', '',
  '  3  1  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0',
  'M  END', '',
].join('\n');

/** Nitrobenzene drawn with a pentavalent N, N(=O)=O, its ring Kekulé (RDKit's layout). */
const NITRO_PENTAVALENT = [
  '', '  nitro-pentavalent', '',
  '  9  9  0  0  0  0  0  0  0  0999 V2000',
  '    3.7500   -1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '    3.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
  '    3.7500    1.2990    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  2  3  2  0', '  2  4  1  0', '  4  5  2  0', '  5  6  1  0', '  6  7  2  0', '  7  8  1  0',
  '  8  9  2  0', '  9  4  1  0',
  'M  END', '',
].join('\n');

/** Methyl azide drawn with a pentavalent N, C-N=N#N. */
const AZIDE_PENTAVALENT = [
  '', '  azide-pentavalent', '',
  '  4  3  0  0  0  0  0  0  0  0999 V2000',
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
  '    2.5981    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
  '    3.8971    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  1  0', '  2  3  2  0', '  3  4  3  0',
  'M  END', '',
].join('\n');

/** Pyridine N-oxide drawn with a pentavalent N, N=O, its ring Kekulé. */
const N_OXIDE_PENTAVALENT = [
  '', '  n-oxide-pentavalent', '',
  '  7  7  0  0  0  0  0  0  0  0999 V2000',
  '    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  2  3  2  0', '  3  4  1  0', '  4  5  2  0', '  5  6  1  0', '  6  7  2  0', '  7  2  1  0',
  'M  END', '',
].join('\n');

/** Indole as OpenChemLib writes [nH] (its V2000, Kekulé): an H count query (hhh 2, at least one H) on the ring NH. */
const OCL_INDOLE_HCOUNT = [
  '', 'Actelion Java MolfileCreator 1.0', '',
  '  9 10  0  0  0  0  0  0  0  0999 V2000',
  '    0.4158    1.9563   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.5352    1.6473   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -0.7431    0.6691   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000   -0.0000   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000   -1.0000   -0.0000 N   0  0  0  2  0  0  0  0  0  0  0  0',
  '    0.9511   -1.3090   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.5388   -0.5000   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.9511    0.3090   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.1590    1.2872   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0  0  0  0', '  2  3  1  0  0  0  0', '  3  4  2  0  0  0  0', '  4  5  1  0  0  0  0',
  '  5  6  1  0  0  0  0', '  6  7  2  0  0  0  0', '  7  8  1  0  0  0  0', '  4  8  1  0  0  0  0',
  '  8  9  2  0  0  0  0', '  1  9  1  0  0  0  0',
  'M  END', '',
].join('\n');

// ---------------------------------------------------------------- the cases

const BENZENE = '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1';
const NITRO_COMPLEX = '[#8-]-[#7+](=[#8])-[#6]1:[#6;x2]:[#6;h{1-}]:[#6,#7]:[#6]:[#6]:1';
const TRIAZOLE_COMPLEX = '[#7]1:[#7]:[#6](-[#6,#7]):[#7](-[#6;x0]):[#6]:1-[#6;a]';
const PYRROLIDINE_COMPLEX = '[#6]-[#6]1(-[#7,#8;D2])-[#6]-[#6]-[#7;x2]-[#6;R1]-1';
const ALANINE = '[#6]-[#6](-[#7])-[#6](=[#8])-[#8]';
const query = (q: string, ...items: string[]): string[] => ['query', `query.${q}`, ...items.map((i) => `query.${q}.${i}`)];
/** Indole, its ring N with at least one H (OpenChemLib's [nH]). */
const INDOLE_NH = '[#6]1:[#6]:[#6]:[#6]2:[#7&h{1-}]:[#6]:[#6]:[#6]:2:[#6]:1';
/** [the molecule, its dataset, Crux's drawing (the last bond to the N added by Crux's bond tool), the molfile Ketcher and
 * OpenChemLib are given, the query intended (charge-separated)]. No pyridine N-oxide: none of the datasets holds one
 * (`a pentavalent N drawn matches the molecule holding it charge-separated` has it). */
const PENTAVALENT_N: [string, Ds, Draw, string, string][] = [
  ['a nitrobenzene', 'smiles.csv', {crux: 'O=Nc1ccccc1', steps: [['tool', 'bond.double', 'atom', 1], ['label', 8, 'O']]},
    NITRO_PENTAVALENT, `[#8-]-[#7+](=[#8])-${BENZENE}`],
  ['methyl azide', 'mol1K', {crux: 'CN=N', steps: [['tool', 'bond.triple', 'atom', 2], ['label', 3, 'N']]},
    AZIDE_PENTAVALENT, '[#6]-[#7]=[#7+]=[#7-]'],
];
/** [a molblock drawn with a pentavalent N, a molecule holding the group charge-separated] */
const PENTAVALENT_N_MOLECULES: [string, string][] = [[NITRO_PENTAVALENT, '[O-][N+](=O)c1ccc(C)cc1'],
  [AZIDE_PENTAVALENT, 'CCN=[N+]=[N-]'], [N_OXIDE_PENTAVALENT, '[O-][n+]1ccc(C)cc1']];

const COMPLEX: Case[] = [
  // Crux, by its own controls
  {name: 'nitrobenzene: a list in the ring, a ring bond count, an H count, charges', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {crux: '[O-][N+](=O)c1ccccc1', steps: [['label', 6, '[C,N]'], ['menu', 'atom', 4, ...query('ring-bond-count', '2')],
      ['menu', 'atom', 5, ...query('h-count', '1')]]}, intended: NITRO_COMPLEX},
  {name: 'a charged atom list', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {crux: 'CN', steps: [['label', 1, '[N,O]'], ['menu', 'atom', 1, 'charge', 'charge.plus']]}, intended: '[#6]-[#7,#8;+]'},
  {name: 'a NOT list on a fused ring', ds: 'mol1K', sketcher: 'Crux',
    draw: {crux: 'c1ccc2ncccc2c1', steps: [['label', 4, '[!C,S]']]}, intended: NOT_QUINOLINE_SMARTS},
  {name: 'any bonds around a ring', ds: 'SPGI', sketcher: 'Crux',
    draw: {crux: 'C1CCCCC1', steps: [0, 1, 2, 3, 4, 5].map((b): CruxStep => ['menu', 'bond', b, 'bond-query', 'bond-query.any'])},
    intended: ANY_RING_SMARTS},
  {name: 'aromatic, ring size and ring bond count queries on a fused system', ds: 'mol1K', sketcher: 'Crux',
    draw: {crux: 'c1ccc2ccccc2c1', steps: [['menu', 'atom', 0, ...query('aromaticity', 'aromatic')],
      ['menu', 'atom', 3, ...query('ring-bond-count', '3')], ['menu', 'atom', 5, ...query('ring-size', '6')]]},
    intended: '[#6;a]1:[#6]:[#6]:[#6;x3]2:[#6]:[#6;r6]:[#6]:[#6]:[#6]:2:[#6]:1'},
  {name: 'SPGI\'s pyrrolidine with four features (a list with a substitution count, a ring bond count, a ring membership)',
    ds: 'SPGI', sketcher: 'Crux', draw: {crux: 'CC1(N)CCNC1', steps: [['label', 2, '[N,O]'],
      ['menu', 'atom', 2, ...query('substitution', '2')], ['menu', 'atom', 5, ...query('ring-bond-count', '2')],
      ['menu', 'atom', 6, ...query('ring-membership', '1')]]}, intended: PYRROLIDINE_COMPLEX},
  {name: 'SPGI\'s triazole: a list, no ring bonds, an aromatic atom', ds: 'SPGI', sketcher: 'Crux',
    draw: {crux: 'n1nc(C)n(C)c1C', steps: [['label', 3, '[C,N]'], ['menu', 'atom', 5, ...query('ring-bond-count', 'as-drawn')],
      ['menu', 'atom', 7, ...query('aromaticity', 'aromatic')]]}, intended: TRIAZOLE_COMPLEX},
  {name: 'two fragments, one a list', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {crux: 'c1ccccc1.Cl', steps: [['label', 6, '[Cl,Br]']]}, intended: `${BENZENE}.[#17,#35]`},
  {name: 'an R label (RDKit: any atom)', ds: 'mol1K', sketcher: 'Crux',
    draw: {crux: '[*:1]c1ccccc1O'}, intended: '*-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]'},
  {name: 'an attachment point (RDKit: none)', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {crux: 'Cc1ccccc1', steps: [['tool', 'rgroup.attachment', 'atom', 0]]}, intended: `[#6]-${BENZENE}`},
  {name: 'an explicit H (RDKit: at least one H)', ds: 'mol1K', sketcher: 'Crux',
    draw: {crux: '[H]OC(=O)c1ccccc1'}, intended: '[#8;!H0]-[#6](=[#8])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'},
  {name: 'an isotope', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: '[13CH3]c1ccccc1'}, intended: {mol: '[13CH3]c1ccccc1'}},
  {name: 'a radical', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: '[CH2]c1ccccc1'}, intended: {mol: '[CH2]c1ccccc1'}},
  {name: 'a wedge (the search ignores chirality, as RDKit does by default)', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {crux: 'C[C@H](N)C(=O)O'}, intended: ALANINE},
  {name: 'an AND group (the search ignores it, as RDKit does by default)', ds: 'mol1K', sketcher: 'Crux',
    draw: {crux: 'C[C@H](N)C(=O)O |&1:1|', format: 'cxsmiles'}, intended: ALANINE},
  // Ketcher, its structure as its tools leave it
  {name: 'nitrobenzene: a list in a Kekulé ring, a ring bond count, an H count, charges', ds: 'smiles.csv', sketcher: 'Ketcher',
    draw: {ket: [] as any}, intended: NITRO_COMPLEX},
  {name: 'a NOT list on a fused Kekulé ring', ds: 'mol1K', sketcher: 'Ketcher', draw: {ket: [] as any}, intended: NOT_QUINOLINE_SMARTS},
  {name: 'any bonds around a ring', ds: 'SPGI', sketcher: 'Ketcher', draw: {ket: [] as any}, intended: ANY_RING_SMARTS},
  {name: 'two fragments, one a list', ds: 'smiles.csv', sketcher: 'Ketcher', draw: {ket: [] as any}, intended: `${BENZENE}.[#17,#35]`},
  {name: 'SPGI\'s triazole: a list, no ring bonds, an aromatic atom', ds: 'SPGI', sketcher: 'Ketcher', draw: {ket: [] as any},
    intended: TRIAZOLE_COMPLEX},
  // OpenChemLib, its editor given a query molfile
  {name: 'any bonds around a ring', ds: 'SPGI', sketcher: 'OpenChemLib', draw: {ocl: V2000_ANY_RING}, intended: ANY_RING_SMARTS},
  {name: 'a NOT list on a fused Kekulé ring', ds: 'mol1K', sketcher: 'OpenChemLib', draw: {ocl: V2000_NOT_QUINOLINE},
    intended: NOT_QUINOLINE_SMARTS},
  // the product owner (2026-10-14): an H count on a ring NH in a Kekulé V2000 found no rows (RDKit's V2000 reader gave the
  // atom of an hhh no implicit H: read sanitized, the N was a radical and its ring not aromatic)
  {name: 'an H count on the ring NH of a Kekulé indole (its V2000)', ds: 'mol1K', sketcher: 'OpenChemLib',
    draw: {ocl: OCL_INDOLE_HCOUNT}, intended: INDOLE_NH, through: []},
  // the product owner (2026-10-14): a pentavalent N drawn (a nitro N(=O)=O, an azide N=N#N, an N-oxide N=O) found no rows,
  // the molecules holding it charge-separated; the query intended is the charge-separated one. Not switched through the
  // other sketchers: each writes the N in its own form, which the search reads alike
  ...PENTAVALENT_N.flatMap(([what, ds, crux, molfile, intended]): Case[] => [
    {name: `${what} drawn with a pentavalent N`, ds, sketcher: 'Crux', draw: crux, intended, through: []},
    {name: `${what} drawn with a pentavalent N`, ds, sketcher: 'Ketcher', draw: {ketcher: molfile}, intended, through: []},
    {name: `${what} drawn with a pentavalent N`, ds, sketcher: 'OpenChemLib', draw: {ocl: molfile}, intended, through: []},
  ]),
];

/** Ketcher's structures, laid out by RDKit once its module is in (`before`). */
function ketcherStructures(): void {
  const [nitro, quinoline, ring, fragments, triazole] = COMPLEX.filter((c) => 'ket' in c.draw);
  nitro.draw = {ket: ketFrom('[O-][N+](=O)c1ccccc1', {0: {charge: -1}, 1: {charge: 1},
    6: {type: 'atom-list', elements: ['C', 'N']}, 4: {ringBondCount: 2}, 5: {hCount: 2}})};
  quinoline.draw = {ket: ketFrom('c1ccc2ncccc2c1', {4: {type: 'atom-list', elements: ['C', 'S'], notList: true}})};
  ring.draw = {ket: ketFrom('C1CCCCC1', {}, Object.fromEntries([0, 1, 2, 3, 4, 5].map((b) => [b, {type: 8}])))};
  fragments.draw = {ket: ketFrom('c1ccccc1.Cl', {6: {type: 'atom-list', elements: ['Cl', 'Br']}})};
  triazole.draw = {ket: ketFrom('n1nc(C)n(C)c1C', {3: {type: 'atom-list', elements: ['C', 'N']}, 5: {ringBondCount: -2},
    7: {queryProperties: {aromaticity: 'aromatic'}}})};
}

const AMIDE = '[NX3][CX3]=[OX1]';
const TEXT: Case[] = [
  {name: 'a plain SMARTS', ds: 'mol1K', sketcher: 'Crux', draw: {field: '[#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'},
    intended: '[#7]1:[#6]:[#6]:[#6]:[#6]:[#6]:1'},
  {name: 'a SMARTS with a list', ds: 'smiles.csv', sketcher: 'Crux', draw: {field: '[Cl,Br,I]-c1ccccc1'}, intended: '[Cl,Br,I]-c1ccccc1'},
  {name: 'a recursive SMARTS', ds: 'SPGI', sketcher: 'Crux', draw: {field: '[#6;$([#6]=[#8])]-[#7]'}, intended: '[#6;$([#6]=[#8])]-[#7]'},
  {name: 'a SMARTS whose ; binds last ([#6,#7;R])', ds: 'mol1K', sketcher: 'Crux', draw: {field: '[#6,#7;R]-[#8]'}, intended: '[#6,#7;R]-[#8]'},
  {name: 'a SMARTS whose & binds first ([#6,#7&R])', ds: 'mol1K', sketcher: 'Crux', draw: {field: '[#6,#7&R]-[#8]'}, intended: '[#6,#7&R]-[#8]'},
  {name: 'a SMARTS with ring closures and a ring bond (@)', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {field: '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]@[#7]'}, intended: '[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]@[#7]'},
  {name: 'a SMARTS of two recursive alternatives ($(),$())', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {field: '[$([#7;H2]-c),$([#8;H1]-c)]'}, intended: '[$([#7;H2]-c),$([#8;H1]-c)]'},
  {name: 'a SMARTS of primitives alone (Daylight\'s amide, [NX3][CX3]=[OX1])', ds: 'SPGI', sketcher: 'Crux',
    draw: {field: AMIDE}, intended: AMIDE},
  {name: 'a SMARTS whose only query syntax is an any bond (c1ccccc1~O)', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {field: 'c1ccccc1~O'}, intended: 'c1ccccc1~O'},
  {name: 'a V3000 molblock with a list, SUBST, RBCNT and UNSAT', ds: 'SPGI', sketcher: 'Crux', draw: {field: V3000_TRIAZOLE},
    intended: V3000_TRIAZOLE_SMARTS},
  {name: 'a V3000 molblock with a NOT list', ds: 'mol1K', sketcher: 'Crux', draw: {field: V3000_NOT_QUINOLINE},
    intended: NOT_QUINOLINE_SMARTS},
  {name: 'a V2000 molblock with an M  ALS list and an M  RBC count', ds: 'smiles.csv', sketcher: 'Crux',
    draw: {field: V2000_ALS_RBC}, intended: V2000_ALS_RBC_SMARTS},
  {name: 'a CXSMILES with a query label (Q_e)', ds: 'smiles.csv', sketcher: 'Crux', draw: {field: 'Cc1ccccc1 |$;;;;Q_e;;$|'},
    intended: {mol: 'Cc1ccccc1 |$;;;;Q_e;;$|'}},
  {name: 'a SMARTS typed and entered', ds: 'smiles.csv', sketcher: 'Crux', draw: {field: '[#6]-[#7,#8]', by: 'enter'},
    intended: '[#6]-[#7,#8]'},
  {name: 'a SMARTS with a list', ds: 'smiles.csv', sketcher: 'Ketcher', draw: {field: '[Cl,Br,I]-c1ccccc1'}, intended: '[Cl,Br,I]-c1ccccc1'},
  {name: 'a recursive SMARTS', ds: 'SPGI', sketcher: 'Ketcher', draw: {field: '[#6;$([#6]=[#8])]-[#7]'}, intended: '[#6;$([#6]=[#8])]-[#7]',
    skip: 'Ketcher (Indigo\'s SMARTS reading): the recursive $() is dropped, its atom a plain C with no custom query, so the drawing is [#6]-[#7] (100 of SPGI\'s rows, 88 wanted)'},
  {name: 'a SMARTS of primitives alone (Daylight\'s amide)', ds: 'SPGI', sketcher: 'Ketcher', draw: {field: AMIDE}, intended: AMIDE},
  {name: 'a V3000 molblock with a list, SUBST, RBCNT and UNSAT', ds: 'SPGI', sketcher: 'Ketcher', draw: {field: V3000_TRIAZOLE},
    intended: V3000_TRIAZOLE_SMARTS},
  {name: 'a V2000 molblock with an M  ALS list and an M  RBC count', ds: 'smiles.csv', sketcher: 'OpenChemLib',
    draw: {field: V2000_ALS_RBC}, intended: V2000_ALS_RBC_SMARTS},
  {name: 'a SMARTS of primitives alone (Daylight\'s amide)', ds: 'SPGI', sketcher: 'OpenChemLib', draw: {field: AMIDE}, intended: AMIDE,
    skip: 'OpenChemLib: it holds no connectivity (X) query, and the query molblock RDKit writes for one (its SMARTS in data S-groups) it reads as halogen lists and isotopes (search [#7]-[#9,#17,#35,#53,#85]-[3#6]-...)'},
];

const label = (atom: number, text: string): CruxStep => ['label', atom, text];
const UNREALISTIC: Case[] = [
  {name: 'a pentavalent carbon', ds: 'mol1K', sketcher: 'Crux', draw: {crux: 'CC(C)(C)(C)C'},
    intended: '[#6]-[#6](-[#6])(-[#6])(-[#6])-[#6]'},
  {name: 'a hypervalent nitrogen', ds: 'mol1K', sketcher: 'Crux', draw: {crux: 'CN(C)(C)(C)C'},
    intended: '[#6]-[#7](-[#6])(-[#6])(-[#6])-[#6]'},
  {name: 'a carbon with charge +3', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: 'C[C+3]C'}, intended: '[#6]-[#6+3]-[#6]'},
  {name: 'a lone R label', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: '[*:1]'}, intended: '*'},
  {name: 'a one-atom query', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: 'C', steps: [label(0, 'N')]}, intended: '[#7]'},
  {name: 'generics alone (A-A)', ds: 'mol1K', sketcher: 'Crux', draw: {crux: 'CC', steps: [label(0, 'A'), label(1, 'A')]},
    intended: '[!#1]-[!#1]'},
  {name: 'generics alone (Q~*)', ds: 'SPGI', sketcher: 'Crux', draw: {crux: 'CC', steps: [label(0, 'Q'), label(1, '*'),
    ['menu', 'bond', 0, 'bond-query', 'bond-query.any']]}, intended: '[!#6&!#1]~*'},
  {name: 'disconnected atoms (Cl.Cl)', ds: 'mol1K', sketcher: 'Crux', draw: {crux: 'Cl.Cl'}, intended: '[#17].[#17]'},
  {name: 'a 20-membered ring', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: 'C1CCCCCCCCCCCCCCCCCCC1'},
    intended: '[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-1'},
  {name: 'aromatic bonds outside a ring', ds: 'mol1K', sketcher: 'Crux', draw: {crux: CHAIN_AROMATIC}, intended: '[#6]:[#6]:[#6]'},
  {name: 'a Kekulé ring partly marked aromatic', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: PARTLY_AROMATIC}, intended: BENZENE},
  {name: 'a dummy atom (*)', ds: 'smiles.csv', sketcher: 'Crux', draw: {crux: '*c1ccccc1'}, intended: `*-${BENZENE}`,
    skip: 'RDKit: a SMILES * is a dummy atom ([#0] in its SMARTS), a molfile * (RDKit writes it R) any atom; Crux reads both as RDKit does, so its SMARTS of the drawing turns from [#0] to * on the reopening (the search, by the molfile, reads * both times)'},
  {name: 'atoms drawn on each other', ds: 'SPGI', sketcher: 'Crux', draw: {crux: OVERLAPPING}, intended: '[#6]-[#6].[#8]'},
  {name: 'an impossible H count (H5)', ds: 'mol1K', sketcher: 'Crux', draw: {crux: 'CC', steps: [['atomCustom', 0, '[#6;H5]']]},
    intended: '[#6;H5]-[#6]'},
  {name: 'a pentavalent carbon', ds: 'mol1K', sketcher: 'Ketcher', draw: {ket: [
    [{label: 'C', location: [0, 0, 0]}, ...[0, 72, 144, 216, 288].map((deg) =>
      ({label: 'C', location: [1.5 * Math.cos(deg * Math.PI / 180), 1.5 * Math.sin(deg * Math.PI / 180), 0]}))],
    [1, 2, 3, 4, 5].map((i) => ({type: 1, atoms: [0, i]}))]}, intended: '[#6]-[#6](-[#6])(-[#6])(-[#6])-[#6]'},
];

/** A query drawn; then a SMARTS typed in the field and entered over it: no warning (the host compared the text with the
 * drawing's SMILES, which a query has none of), the drawing and the rows those of the SMARTS. */
async function enteredOverQuery(type: SketcherType, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE['smiles.csv'];
  const smarts = '[#6]-[#7,#8;!R]';
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const mols = molsOf(file, df.col(column)!);
  const want = intendedRows(mols, smarts);
  at.stage = 'drawing a query';
  await openDialog(filter, type);
  await draw(filter, type === 'Crux' ? {crux: 'CN', steps: [label(1, '[N,O]')]} :
    {ket: ketFrom('CN', {1: {type: 'atom-list', elements: ['N', 'O']}})});
  DG.Balloon.closeAll();
  at.stage = 'entering a SMARTS over it';
  await draw(filter, {field: smarts, by: 'enter'});
  // a warning's balloon lands after the call that shows it
  await delay(500);
  const warned = notationBalloons();
  if (warned)
    problems.push(`entered over a query: a warning, ${warned}`);
  const h = await held(filter);
  if (!sameRows(drawingRows(mols, h.molfile), want))
    problems.push(`entered over a query: the drawing (${h.search}) is not ${smarts}`);
  pressDialogButton(filter, 'OK');
  await settle(filter);
  if (!sameRows(rows(df.filter), want))
    problems.push(`after OK: ${rows(df.filter).length} rows pass; ${smarts} matches ${want.length}`);
  return problems;
}

/** A CXSMILES pasted in the field of a molecule sketcher (a dialog's, as a script's or a form's): the sketcher shows that
 * molecule, its enhanced stereo kept, as a molecule (not a query), with no warning. */
async function pastedInSketcher(type: SketcherType, cxsmiles: string, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  grok.chem.currentSketcherType = type;
  const host = new DG.chem.Sketcher();
  const dialog = ui.dialog('a molecule').add(host).show();
  try {
    at.stage = `${type} getting ready`;
    await within(host.sketcherReady(), `${type} did not get ready`);
    DG.Balloon.closeAll();
    at.stage = 'pasting';
    const data = new DataTransfer();
    data.setData('text/plain', cxsmiles);
    host.molInput.dispatchEvent(new ClipboardEvent('paste', {clipboardData: data, bubbles: true, cancelable: true}));
    if (isKetcher(host.sketcher))
      await ketcherSettled(host.sketcher);
    if (host._smarts !== null)
      problems.push(`the host took the CXSMILES for SMARTS`);
    const impl: any = host.sketcher;
    if (isCrux(impl)) {
      const el = impl.root.querySelector('crux-sketch');
      await waitUntil(() => !el.isPending, 'Crux did not settle', 5000);
      if (el.hasQuery)
        problems.push(`Crux shows a query: ${el.smarts}`);
    }
    at.stage = 'reading';
    const smiles = host.getSmiles();
    if (canonical(smiles) !== canonical(cxsmiles))
      problems.push(`the sketcher holds ${canonical(smiles)}; pasted ${canonical(cxsmiles)}`);
    const warned = notationBalloons();
    if (warned)
      problems.push(`a warning: ${warned}`);
  } finally {
    dialog.close();
  }
  return problems;
}

// ---------------------------------------------------------------- routes

/** The searches the filter starts and the change events its sketcher host fires, from now on. */
function countActivity(filter: SubstructureFilter): {searches: () => number, changes: () => number, stop: () => void} {
  let searches = 0;
  let changes = 0;
  const original = filter.getFilterBitset;
  filter.getFilterBitset = function(this: SubstructureFilter) {
    searches++;
    return original.call(this);
  };
  const sub = filter.sketcher.onChanged.subscribe(() => changes++);
  return {searches: () => searches, changes: () => changes, stop: () => {
    filter.getFilterBitset = original;
    sub.unsubscribe();
  }};
}

/** Drawn, OK; reopened and OK unchanged: the rows the same, no search started, the host's change events few and ended. */
async function reopenUnchanged(type: SketcherType, d: Draw, intended: string, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE['smiles.csv'];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const want = intendedRows(molsOf(file, df.col(column)!), intended);
  at.stage = 'drawing';
  await openDialog(filter, type);
  await draw(filter, d);
  pressDialogButton(filter, 'OK');
  await settle(filter);
  if (!sameRows(rows(df.filter), want))
    problems.push(`after OK: ${rows(df.filter).length} rows pass; ${intended} matches ${want.length}`);
  const before = filter.currentMolecule;
  const activity = countActivity(filter);
  try {
    at.stage = 'reopening';
    await openDialog(filter, type);
    if (isKetcher(filter.sketcher.sketcher))
      await ketcherSettled(filter.sketcher.sketcher);
    await waitUntil(() => filter.sketcher.sketcher?.isRenderPending !== true,
      'the reopened sketcher did not settle', 5000);
    pressDialogButton(filter, 'OK');
    at.stage = 'the search after the second OK';
    await settle(filter);
    await delay(1000);
    const changes = activity.changes();
    await delay(1000);
    if (activity.changes() !== changes)
      problems.push(`the host's change events go on after the dialog closed: ${changes}, then ${activity.changes()}`);
    if (changes > 3)
      problems.push(`reopened and OK unchanged: ${changes} change events`);
    if (activity.searches() !== 0)
      problems.push(`reopened and OK unchanged: ${activity.searches()} searches started`);
    if (searchReads(filter.currentMolecule) !== searchReads(before))
      problems.push(`the filter's query changed: ${searchReads(before)} to ${searchReads(filter.currentMolecule)}`);
    if (!sameRows(rows(df.filter), want))
      problems.push(`after the second OK: ${rows(df.filter).length} rows pass; ${intended} matches ${want.length}`);
  } finally {
    activity.stop();
  }
  return problems;
}

/** Drawn, OK; reopened, a second feature added, OK: the rows those of both; reopened, both there. */
async function editAfterReopen(type: SketcherType, first: Draw, second: Draw, intended: string,
  at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE['smiles.csv'];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const mols = molsOf(file, df.col(column)!);
  const want = intendedRows(mols, intended);
  at.stage = 'drawing the first feature';
  await openDialog(filter, type);
  await draw(filter, first);
  pressDialogButton(filter, 'OK');
  await settle(filter);
  at.stage = 'reopening';
  await openDialog(filter, type);
  if (isKetcher(filter.sketcher.sketcher))
    await ketcherSettled(filter.sketcher.sketcher);
  at.stage = 'adding the second feature';
  await draw(filter, second);
  const edited = await held(filter);
  if (!sameRows(drawingRows(mols, edited.molfile), want))
    problems.push(`edited: the drawing (${edited.search}) is not ${intended}`);
  pressDialogButton(filter, 'OK');
  at.stage = 'the search after the edit';
  await settle(filter);
  if (!sameRows(rows(df.filter), want))
    problems.push(`after the edit's OK: ${rows(df.filter).length} rows pass; ${intended} matches ${want.length}`);
  at.stage = 'reopening after the edit';
  await openDialog(filter, type);
  if (isKetcher(filter.sketcher.sketcher))
    await ketcherSettled(filter.sketcher.sketcher);
  const reopened = await held(filter);
  if (heldText(reopened) !== heldText(edited))
    problems.push(`reopened after the edit: ${heldText(reopened)}; edited: ${heldText(edited)}`);
  pressDialogButton(filter, 'CANCEL');
  await settle(filter);
  return problems;
}

/** Drawn, OK; the view's layout saved, the table opened afresh and the layout applied: the restored card filters the same
 * rows, and its sketcher opens with the query drawn. */
async function throughLayout(c: Case, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE[c.ds];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  const want = intendedRows(molsOf(file, df.col(column)!), c.intended);
  at.stage = 'drawing';
  await openDialog(filter, c.sketcher);
  await draw(filter, c.draw);
  const drawn = await held(filter);
  pressDialogButton(filter, 'OK');
  await settle(filter);
  if (!sameRows(rows(df.filter), want))
    problems.push(`after OK: ${rows(df.filter).length} rows pass; ${describe(c.intended)} matches ${want.length}`);
  at.stage = 'saving the layout';
  const layout = grok.shell.tv.saveLayout();
  grok.shell.closeAll();
  const df2 = tables.get(file)!.clone();
  df2.name = df.name;
  const tv2 = grok.shell.addTableView(df2);
  at.stage = 'applying the layout';
  tv2.loadLayout(layout);
  let restored: SubstructureFilter | undefined;
  await waitUntil(() => {
    restored = tv2.getFiltersGroup({createDefaultFilters: false}).filters
      .find((f: any) => f?.sketcher instanceof DG.chem.Sketcher && f.columnName === column) as unknown as SubstructureFilter;
    return restored !== undefined && restored.laidOut;
  }, 'the layout restored no substructure filter card', 20000);
  at.stage = 'the restored card\'s search';
  await waitUntil(() => sameRows(rows(df2.filter), want) && !restored!.calculating && restored!.currentSearches.size === 0,
    'the restored card\'s search did not end', 20000).catch(() => problems.push(
    `the restored card passes ${rows(df2.filter).length} rows; ${describe(c.intended)} matches ${want.length}`));
  at.stage = 'opening the restored card\'s sketcher';
  await openDialog(restored!, c.sketcher);
  if (isKetcher(restored!.sketcher.sketcher))
    await ketcherSettled(restored!.sketcher.sketcher);
  const reopened = await held(restored!);
  if (heldText(reopened) !== heldText(drawn))
    problems.push(`opened from the restored card: ${heldText(reopened)}; drawn: ${heldText(drawn)}`);
  pressDialogButton(restored!, 'CANCEL');
  await settle(restored!);
  return problems;
}

/** A query drawn, OK; then cleared (in the sketcher, or by the card's Clear button), OK: every row passes, the card shows
 * the Sketch link, and the sketcher reopens empty, with no warning. */
async function cleared(how: 'sketcher' | 'card', type: SketcherType, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  const [file, column] = FILE['smiles.csv'];
  at.stage = 'opening the table';
  const {df, filter} = await openTable(file, column);
  at.stage = 'drawing';
  await openDialog(filter, type);
  await draw(filter, type === 'Crux' ? {crux: 'CN', steps: [label(1, '[N,O]')]} : {field: '[#6]-[#7,#8]'});
  pressDialogButton(filter, 'OK');
  await settle(filter);
  if (df.filter.trueCount === df.rowCount)
    problems.push('the query filtered nothing');
  at.stage = 'clearing';
  if (how === 'card') {
    // the card draws the query asynchronously and shows its Clear button with it
    await waitUntil(() => filter.sketcher.extSketcherDiv.contains(filter.sketcher.clearSketcherButton),
      'the card shows no Clear button for the query', 10000);
    filter.sketcher.clearSketcherButton.click();
  }
  else {
    await openDialog(filter, type);
    if (type === 'Crux')
      shadow(cruxOf(filter), 'toolbar.clear').click();
    else
      throw new Error(`no clearing in ${type}`);
    await settle(filter);
    pressDialogButton(filter, 'OK');
  }
  await settle(filter);
  await waitUntil(() => df.filter.trueCount === df.rowCount,
    `cleared: ${df.filter.trueCount} of ${df.rowCount} rows pass`, 10000)
    .catch((e) => problems.push(e?.message ?? String(e)));
  await waitUntil(() => filter.sketcher.extSketcherDiv.querySelector('.sketch-link') !== null, 'no Sketch link', 5000)
    .catch(() => problems.push('the card shows no Sketch link after the clearing'));
  DG.Balloon.closeAll();
  at.stage = 'reopening';
  await openDialog(filter, type);
  const reopened = await held(filter);
  if (!DG.chem.Sketcher.isEmptyMolfile(reopened.molfile))
    problems.push(`reopened after the clearing: ${heldText(reopened)}`);
  const warned = notationBalloons();
  if (warned)
    problems.push(`reopened with a warning: ${warned}`);
  pressDialogButton(filter, 'CANCEL');
  await settle(filter);
  return problems;
}

// ---------------------------------------------------------------- the molecule cell editor

/** RDKit's canonical CXSMILES of a molecule string, without coordinates: what two strings of one molecule share. */
function canonical(molecule: string): string {
  const m = rdkit.get_mol(molecule);
  if (!m)
    return `unreadable: ${molecule}`;
  try {
    return m.get_cxsmiles(JSON.stringify({CX_ALL_BUT_COORDS: true}));
  } finally {
    m.delete();
  }
}

/** `molecule` with its atom `atom` (0-based) made `element`, as RDKit's V3000 (which keeps enhanced stereo). */
function withElement(molecule: string, atom: number, element: string): string {
  const m = rdkit.get_mol(molecule);
  try {
    return molfileWithElement(m.get_v3Kmolblock(), atom, element);
  } finally {
    m.delete();
  }
}

/** A molfile, V2000 or V3000, with its atom `atom` (0-based) made `element`. */
function molfileWithElement(molfile: string, atom: number, element: string): string {
  const lines = molfile.split('\n');
  if (molfile.includes('V3000')) {
    const start = lines.findIndex((l) => l.trim() === 'M  V30 BEGIN ATOM') + 1;
    const parts = lines[start + atom].split(' ');
    parts[4] = element;
    lines[start + atom] = parts.join(' ');
  } else {
    const l = lines[4 + atom];
    lines[4 + atom] = l.substring(0, 31) + element.padEnd(3, ' ') + l.substring(34);
  }
  return lines.join('\n');
}

/** The cell editor's host once its sketcher is ready. */
async function cellEditorHost(): Promise<DG.chem.Sketcher> {
  let host: DG.chem.Sketcher | null = null;
  await waitUntil(() => {
    const dialogs = DG.Dialog.getOpenDialogs();
    const box = dialogs[dialogs.length - 1]?.root.querySelector('.chem-sketcher-host');
    const child = box?.firstElementChild;
    const impl = child ? DG.Widget.find(child) as DG.chem.SketcherBase | null : null;
    host = impl?.host ?? null;
    return host !== null;
  }, 'the cell editor shows no sketcher', 15000);
  await within(host!.sketcherReady(), 'the cell editor\'s sketcher did not get ready');
  if (isKetcher(host!.sketcher))
    await ketcherSettled(host!.sketcher);
  return host!;
}

/** The molecule the cell editor's sketcher itself shows, as it writes it for the cell: its SMILES (Crux's own, Ketcher's
 * Indigo export, OpenChemLib's), or its molfile for a molblock column. */
async function shownMolecule(host: DG.chem.Sketcher, units: string): Promise<string> {
  const impl: any = host.sketcher;
  if (isCrux(impl)) {
    const el = impl.root.querySelector('crux-sketch') as CruxElement;
    await waitUntil(() => !el.isPending, 'Crux did not settle', 5000);
    return units === DG.chem.Notation.MolBlock ? (el as any).molV3000 : (el as any).smiles;
  }
  if (isKetcher(impl)) {
    await waitUntil(() => impl._smiles !== null || impl._molV3000 !== null,
      'Ketcher exported nothing', 10000).catch(() => {});
    return units === DG.chem.Notation.MolBlock ? impl._molV3000 ?? impl._molV2000 : impl._smiles ?? impl._molV2000;
  }
  return units === DG.chem.Notation.MolBlock ? impl._sketcher.getMolFileV3() : impl._sketcher.getSmiles();
}

/** Edits the cell editor's drawing by the sketcher's own means: atom `atom` made `element` (Crux: its element tool;
 * Ketcher and OpenChemLib: their own molfile of what they show, so edited, loaded as their own). */
async function editAtom(host: DG.chem.Sketcher, atom: number, element: string): Promise<void> {
  const impl: any = host.sketcher;
  if (isCrux(impl)) {
    const el = impl.root.querySelector('crux-sketch') as CruxElement;
    await cruxSteps(el, [['tool', `element.${element.toLowerCase()}`, 'atom', atom]]);
  } else if (isKetcher(impl)) {
    await waitUntil(() => impl._molV3000 !== null, 'Ketcher exported no V3000', 10000);
    const exported = impl._exportId;
    await within(impl._sketcher.setMolecule(molfileWithElement(impl._molV3000, atom, element)),
      'Ketcher did not load the edit', 15000);
    await ketcherExported(impl, 5000, exported);
  } else
    impl._sketcher.setMolFile(molfileWithElement(impl._sketcher.getMolFile(), atom, element));
}

async function pressOK(host: DG.chem.Sketcher): Promise<void> {
  const dialog = DG.Dialog.getOpenDialogs().find((d) => d.root.contains(host.root))!;
  const done = new Promise<void>((resolve) => {
    const sub = dialog.onAfterOK.subscribe(() => {
      sub.unsubscribe();
      resolve();
    });
  });
  dialog.root.querySelector<HTMLElement>('[name="button-OK"]')!.click();
  await within(done, 'the cell editor\'s OK did not end');
}

/** Double-click on a molecule cell (the cell editor), the molecule shown as the cell holds it; an atom edited, OK: the cell
 * holds the edited molecule; opened again, the sketcher shows the cell's molecule. Compared as RDKit's canonical
 * CXSMILES (stereo, enhanced stereo, isotopes, charges, R labels). */
async function cellEditor(type: SketcherType, value: string, units: string, atom: number, element: string,
  edited: string, at: {stage: string}): Promise<string[]> {
  const problems: string[] = [];
  grok.chem.currentSketcherType = type;
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('molecule', [value])]);
  const col = df.col('molecule')!;
  col.semType = DG.SEMTYPE.MOLECULE;
  col.meta.units = units;
  const recent = localStorage.getItem(DG.chem.Sketcher.RECENT_KEY);
  localStorage.setItem(DG.chem.Sketcher.RECENT_KEY, '[]');
  const tv = grok.shell.addTableView(df);
  try {
    at.stage = 'opening the cell editor';
    await PackageFunctions.editMoleculeCell(tv.grid.cell('molecule', 0));
    let host = await cellEditorHost();
    const shown = await shownMolecule(host, units);
    if (canonical(shown) !== canonical(value))
      problems.push(`the cell editor shows ${canonical(shown)}; the cell holds ${canonical(value)}`);
    at.stage = 'editing';
    await editAtom(host, atom, element);
    at.stage = 'OK';
    await pressOK(host);
    const want = canonical(edited);
    const written = df.get('molecule', 0);
    if (canonical(written) !== want)
      problems.push(`after the edit the cell holds ${canonical(written)}; the edit makes ${want}`);
    at.stage = 'opening the cell editor again';
    await PackageFunctions.editMoleculeCell(tv.grid.cell('molecule', 0));
    host = await cellEditorHost();
    const again = await shownMolecule(host, units);
    if (canonical(again) !== canonical(written))
      problems.push(`opened again, the cell editor shows ${canonical(again)}; the cell holds ${canonical(written)}`);
  } finally {
    for (const d of DG.Dialog.getOpenDialogs())
      d.close();
    tv.close();
    if (recent === null)
      localStorage.removeItem(DG.chem.Sketcher.RECENT_KEY);
    else
      localStorage.setItem(DG.chem.Sketcher.RECENT_KEY, recent);
  }
  return problems;
}

/** [the molecule's kind, the cell's SMILES, the atom edited, the element it is made, the SMILES edited] */
const CELL_MOLECULES: [string, string, number, string, string][] = [
  ['stereo', 'N[C@@H](C)C(=O)O', 5, 'N', 'N[C@@H](C)C(=O)N'],
  ['enhanced stereo (an AND group)', 'C[C@H](O)[C@@H](C)CC |&1:1,3|', 6, 'O', 'C[C@H](O)[C@@H](C)CO |&1:1,3|'],
  ['an isotope', '[13CH3]c1ccc(O)cc1', 5, 'N', '[13CH3]c1ccc(N)cc1'],
  ['charges', 'C[N+](C)(C)CCC(=O)[O-]', 4, 'O', 'C[N+](C)(C)OCC(=O)[O-]'],
  ['R labels', '[*:1]c1ccc([*:2])cc1C', 8, 'N', '[*:1]c1ccc([*:2])cc1N'],
];

// ---------------------------------------------------------------- the categories

category('query back and forth: complex queries', () => {
  let sketcherType: string;
  before(async () => {
    sketcherType = await roundTripBefore();
    ketcherStructures();
  });
  after(async () => roundTripAfter(sketcherType));
  caseTests(COMPLEX);
  // the product owner (2026-10-14): a pentavalent N drawn found no rows; read as the search reads a query molblock, each
  // matches a molecule that holds the group charge-separated (as RDKit and the datasets hold it)
  test('a pentavalent N drawn (nitro, azide, N-oxide) matches the molecule holding it charge-separated', async () => {
    const missed: string[] = [];
    for (const [molBlock, smiles] of PENTAVALENT_N_MOLECULES) {
      const q = getQueryMolSafe(molBlock, '', rdkit);
      const mol = rdkit.get_mol(smiles);
      try {
        if (!q || mol.get_substruct_match(q) === '{}')
          missed.push(`${q?.get_smarts() ?? 'no query'} does not match ${smiles}`);
      } finally {
        q?.delete();
        mol.delete();
      }
    }
    expect(missed.length, 0, missed.join('; '));
  });
});

category('query back and forth: entered as text', () => {
  let sketcherType: string;
  before(async () => {
    sketcherType = await roundTripBefore();
  });
  after(async () => roundTripAfter(sketcherType));
  caseTests(TEXT);
  for (const type of ['Crux', 'Ketcher'] as SketcherType[]) {
    test(`smiles.csv: ${type}: a SMARTS typed and entered over a query drawn, with no warning`,
      async () => await expectNoProblems((at) => enteredOverQuery(type, at)), {timeout: 120000});
  }
  for (const type of ['Crux', 'Ketcher'] as SketcherType[]) {
    test(`a molecule sketcher: ${type}: a CXSMILES with an AND group pasted stays that molecule, not a query`,
      async () => await expectNoProblems((at) => pastedInSketcher(type, 'C[C@H](N)C(=O)O |&1:1|', at)), {timeout: 120000});
  }
});

category('query back and forth: non-realistic drawings', () => {
  let sketcherType: string;
  before(async () => {
    sketcherType = await roundTripBefore();
  });
  after(async () => roundTripAfter(sketcherType));
  caseTests(UNREALISTIC);
  for (const how of ['sketcher', 'card'] as const) {
    test(`smiles.csv: Crux: a query cleared ${how === 'sketcher' ? 'in the sketcher' : 'by the card\'s Clear'} filters nothing and reopens empty`,
      async () => await expectNoProblems((at) => cleared(how, 'Crux', at)), {timeout: 120000});
  }
  test('smiles.csv: Ketcher: a query cleared by the card\'s Clear filters nothing and reopens empty',
    async () => await expectNoProblems((at) => cleared('card', 'Ketcher', at)), {timeout: 120000});
});

category('query back and forth: routes', () => {
  let sketcherType: string;
  before(async () => {
    sketcherType = await roundTripBefore();
  });
  after(async () => {
    roundTripAfter(sketcherType);
    freeDatasetMols();
  });

  const LIST = '[#6]-[#7,#8]';
  const reopenDraws: [SketcherType, () => Draw][] = [
    ['Crux', () => ({crux: 'CN', steps: [label(1, '[N,O]')]})],
    ['Ketcher', () => ({ket: ketFrom('CN', {1: {type: 'atom-list', elements: ['N', 'O']}})})],
    ['OpenChemLib', () => ({ocl: [...V2000_ALS_RBC.split('\n').slice(0, 3), '  2  1  0  0  0  0  0  0  0  0999 V2000',
      '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
      '    1.5000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0',
      '  1  2  1  0', 'M  ALS   2  2 F N   O   ', 'M  END', ''].join('\n')})],
  ];
  for (const [type, d] of reopenDraws) {
    test(`smiles.csv: ${type}: reopened and OK unchanged, the filter neither searches again nor loops`,
      async () => await expectNoProblems((at) => reopenUnchanged(type, d(), LIST, at)), {timeout: 120000});
  }

  // switched through every sketcher and back, each holding the query
  caseTests([
    {name: 'a list in a ring, through OpenChemLib and Ketcher and back, each holding it', ds: 'smiles.csv', sketcher: 'Crux',
      draw: {crux: 'Clc1ccccc1', steps: [label(0, '[Cl,Br,I]')]}, intended: `[#17,#35,#53]-${BENZENE}`,
      through: ['OpenChemLib', 'Ketcher'], checkThrough: ['OpenChemLib', 'Ketcher']},
    {name: 'a SMARTS entered, through OpenChemLib and Ketcher and back, Ketcher holding it', ds: 'smiles.csv', sketcher: 'Crux',
      draw: {field: '[Cl,Br,I]-c1ccccc1'}, intended: '[Cl,Br,I]-c1ccccc1', through: ['OpenChemLib', 'Ketcher'], checkThrough: ['Ketcher']},
    {name: 'a SMARTS entered, through OpenChemLib, which holds it', ds: 'smiles.csv', sketcher: 'Crux',
      draw: {field: '[Cl,Br,I]-c1ccccc1'}, intended: '[Cl,Br,I]-c1ccccc1', through: ['OpenChemLib'], checkThrough: ['OpenChemLib'],
      skip: 'OpenChemLib: a SMARTS reaches it as the V2000 query molblock RDKit writes, which has no aromatic-atom query (c becomes #6, its implicit bonds single or aromatic); no molfile OpenChemLib reads holds one'},
    {name: 'a list in a ring, through Crux and OpenChemLib and back, each holding it', ds: 'smiles.csv', sketcher: 'Ketcher',
      draw: {ketcher: V2000_ALS_RBC}, intended: V2000_ALS_RBC_SMARTS, through: ['Crux', 'OpenChemLib'], checkThrough: ['Crux', 'OpenChemLib']},
  ]);

  // edited after a reopening
  const RBC_LIST = '[#6;x2]-[#7,#8]';
  test('smiles.csv: Crux: a feature added after a reopening, both come back', async () => await expectNoProblems((at) =>
    editAfterReopen('Crux', {crux: 'CN', steps: [label(1, '[N,O]')]},
      {cruxSteps: [['menu', 'atom', 0, ...query('ring-bond-count', '2')]]}, RBC_LIST, at)), {timeout: 120000});
  test('smiles.csv: Crux: a feature added to a SMARTS entered, after a reopening, both come back', async () =>
    await expectNoProblems((at) => editAfterReopen('Crux', {field: LIST},
      {cruxSteps: [['menu', 'atom', 0, ...query('ring-bond-count', '2')]]}, RBC_LIST, at)), {timeout: 120000});
  test('smiles.csv: Ketcher: a feature added after a reopening, both come back', async () => await expectNoProblems((at) =>
    editAfterReopen('Ketcher', {ket: ketFrom('CN', {1: {type: 'atom-list', elements: ['N', 'O']}})},
      {ket: ketFrom('CN', {0: {ringBondCount: 2}, 1: {type: 'atom-list', elements: ['N', 'O']}})}, RBC_LIST, at)),
  {timeout: 120000});

  // kept through a layout saved and applied
  const layoutCases: Case[] = [
    {name: 'a query drawn', ds: 'SPGI', sketcher: 'Crux', draw: {crux: 'CC1(N)CCNC1', steps: [label(2, '[N,O]'),
      ['menu', 'atom', 5, ...query('ring-bond-count', '2')]]}, intended: '[#6]-[#6]1(-[#7,#8])-[#6]-[#6]-[#7;x2]-[#6]-1'},
    {name: 'a SMARTS entered', ds: 'SPGI', sketcher: 'Crux', draw: {field: AMIDE}, intended: AMIDE},
    // a SMARTS whose syntax said SMARTS before isSmarts' fix as well: the layout's restore alone (set as a molblock, the
    // card's sketcher opened empty)
    {name: 'a SMARTS with a list entered', ds: 'smiles.csv', sketcher: 'Crux', draw: {field: '[Cl,Br,I]-c1ccccc1'},
      intended: '[Cl,Br,I]-c1ccccc1'},
    {name: 'a query drawn', ds: 'smiles.csv', sketcher: 'Ketcher', draw: {ketcher: V2000_ALS_RBC}, intended: V2000_ALS_RBC_SMARTS},
  ];
  for (const c of layoutCases) {
    test(`${c.ds}: ${c.sketcher}: ${c.name}, kept through the view's layout saved and applied`,
      async () => await expectNoProblems((at) => throughLayout(c, at)), {timeout: 120000});
  }

  // the molecule cell editor
  const cellCases: [SketcherType, number[]][] = [['Crux', [0, 1, 2, 3, 4]], ['Ketcher', [0, 1]], ['OpenChemLib', [2, 3]]];
  for (const [type, which] of cellCases) {
    for (const i of which) {
      const [kind, smiles, atom, element, edited] = CELL_MOLECULES[i];
      test(`the cell editor: ${type}: a molecule with ${kind}, edited, read back the same`, async () =>
        await expectNoProblems((at) => cellEditor(type, smiles, DG.chem.Notation.Smiles, atom, element, edited, at)),
      {timeout: 120000});
    }
  }
  test('the cell editor: Crux: a V3000 molblock with enhanced stereo, edited, read back the same', async () => {
    await useRdKit();
    const m = rdkit.get_mol(CELL_MOLECULES[1][1]);
    const v3000 = m.get_v3Kmolblock();
    m.delete();
    await expectNoProblems((at) => cellEditor('Crux', v3000, DG.chem.Notation.MolBlock, 6, 'O',
      withElement(v3000, 6, 'O'), at));
  }, {timeout: 120000});
});
