/**
 * Category `Foot-guns: Crux`: the checker module's rules (utils/footguns/rules.ts) on the foot-gun catalogue's values;
 * the dialog footer's middle slot (OK and Cancel never move, a block disables OK, a long text is cut and kept whole in
 * the tooltip); the context panel's Checks section; the fixes: MOL-07 (an edit of only what canonical SMILES drop is
 * written), Chem's urea formation template taking two reactants (dr014), TPL-08 (a variable used twice), RXN-13 (the
 * Reaction panel's steps as the value writes them), MAP-02 and MAP-03 (the automatic mapping refusing a template,
 * naming an unreadable component), and what the enumerator returns for TPL-06, TPL-11 and TRF-07 (checked here rather
 * than through the app); crux-sketch's refusal of a ChemDraw reaction in the molecule editor (MOL-17); the reaction
 * cell editor's slot on a reaction Crux Sketch does not open (RXN-06). RDKit mode flags nothing and changes nothing,
 * but for the reaction editor's checks, which run in both modes.
 */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {Subject} from 'rxjs';
import {before, category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {SketcherHandle} from '../../vendor/crux-sketch/crux-sketch';
import {Flag, checkReaction, checkReactionEdit, checkReadLosses, checkSketcherWrite,
  checkTemplateMapping, checkTransform, checkValue, checksModule} from '../utils/footguns/checks';
import {enumeratorFlags} from '../utils/footguns/enumerator';
import {mappingFlags, moleculeFlags, reactantCountFlags, reactionFlags, sketcherWriteFlags}
  from '../utils/footguns/rules';
import {FOOTER_SLOT_CLASS, footerSlot} from '../utils/footguns/footer-slot';
import {unchangedByEdit, watchMoleculeEditor, watchReactionEditor} from '../utils/footguns/editors';
import type {ReactionCellEditor} from '../utils/reactions/reaction-cell-editor';
import {drawingLosses} from '../utils/reactions/reaction-losses';
import {CHECKS_LIST_CLASS, checksPanel, checksShown} from '../widgets/checks-panel';
import {MAPPING_CLASS, reactionPanel, setReactionMapperForTests} from '../widgets/reaction-panel';
import {resolveReactionVariables, runTwoComponentReaction} from '../utils/reactions/reactions';
import {NamedReaction} from '../utils/reactions/types';
import {cloneConfig, DEFAULT_CONFIG} from '../utils/reaction-enumeration/config';
import {enumerate} from '../utils/reaction-enumeration/enumerate';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {_package} from '../package-test';
import {inEngine} from './engine-utils';
import {loadFileAsBytes} from './utils';
import {_importSdf} from '../open-chem/sdf-importer';
import {CruxSketcher, setCruxSketchLoaderForTests} from '../crux/crux-sketcher';

/** `run(module)` in Crux mode, with the checks' module (Crux's). */
const inCrux = <T>(run: (module: RDModule) => Promise<T> | T): Promise<T> =>
  inEngine(ChemicalEngine.Crux, async () => {
    const module = checksModule();
    expect(module !== null, true, 'the checks run on Crux\'s module in Crux mode');
    return run(module!);
  });

/** `run()` in RDKit mode. */
const inRdkit = <T>(run: () => Promise<T> | T): Promise<T> => inEngine(ChemicalEngine.RDKit, async () => run());

/** The flags as one line each, for the messages. */
const said = (flags: Flag[]): string => flags.map((f) => `${f.level}: ${f.text}`).join(' | ') || '(none)';

/** Whether a flag of `level` holds every word (any case). */
function has(flags: Flag[], level: string, words: string[]): boolean {
  return flags.some((f) => f.level === level && words.every((w) => f.text.toLowerCase().includes(w.toLowerCase())));
}

/** Expects a flag of `level` holding every word. */
function expectFlag(name: string, flags: Flag[], level: string, words: string[]): void {
  expect(has(flags, level, words), true, `${name}: a ${level} with ${JSON.stringify(words)}; said: ${said(flags)}`);
}

/** Expects no flag of these levels. */
function expectNone(name: string, flags: Flag[], levels: string[]): void {
  expect(flags.some((f) => levels.includes(f.level)), false, `${name}: no ${levels.join('/')}; said: ${said(flags)}`);
}

/** A molblock of `s` as a query (Crux's), or of a molecule read unsanitized. */
function molblockOf(module: RDModule, s: string, query = false): string {
  const m = query ? module.get_qmol(s) : module.get_mol(s, JSON.stringify({sanitize: false}));
  try {
    return m.get_molblock(JSON.stringify({kekulize: false}));
  } finally {
    m.delete();
  }
}

/** A dialog with OK and Cancel, shown; `body` its content. */
function okCancelDialog(title: string, body: HTMLElement = ui.div('A body')): DG.Dialog {
  return ui.dialog(title).add(body).onOK(() => {}).show();
}

/** OK's and Cancel's places on the page (CSS px). */
function buttonRects(dialog: DG.Dialog): {ok: DOMRect, cancel: DOMRect} {
  const buttons = [...dialog.root.querySelectorAll('.d4-dialog-footer button, .d4-command-bar button')] as
    HTMLButtonElement[];
  const ok = dialog.getButton('OK');
  const cancel = buttons.find((b) => b !== ok && /cancel/i.test(b.textContent ?? ''))!;
  return {ok: ok.getBoundingClientRect(), cancel: cancel.getBoundingClientRect()};
}

/** Whether two rects are the same place, to half a pixel. */
const samePlace = (a: DOMRect, b: DOMRect): boolean =>
  ['left', 'top', 'width', 'height'].every((k) => Math.abs((a as any)[k] - (b as any)[k]) <= 0.5);

/** The next frame (layout done). */
const frame = () => new Promise<void>((resolve) => requestAnimationFrame(() => setTimeout(resolve, 0)));

/** Makes `handle`'s four getters (smiles, molfile, molV3000, smarts) throw `error()`, as a stopped engine's or a
 * drawing the engine cannot write do: own properties over its class's getters, so nothing else of it changes (no
 * Proxy: its methods, destroy among them, keep their own `this`); the function returned puts them back. */
function throwingGetters(handle: object, error: () => unknown): () => void {
  const names = ['smiles', 'molfile', 'molV3000', 'smarts'];
  const saved = names.map((n) => [n, Object.getOwnPropertyDescriptor(handle, n)] as const);
  for (const n of names) {
    Object.defineProperty(handle, n, {configurable: true, get: () => {
      throw error();
    }});
  }
  return () => {
    for (const [n, d] of saved) {
      if (d)
        Object.defineProperty(handle, n, d);
      else
        delete (handle as {[k: string]: unknown})[n];
    }
  };
}

/** What the MOL-17 test uses of Crux's element: crux-sketch's handle. */
type CruxAdd = HTMLElement & Pick<SketcherHandle, 'addMolecule' | 'lastError'>;

/** Resolves once a sketcher has nothing on its way (`isRenderPending`) on two checks in a row, as the "Crux sketcher"
 * category settles. */
async function settled(impl: DG.chem.SketcherBase): Promise<void> {
  let held = 0;
  await waitUntil(() => {
    held = impl.isRenderPending ? 0 : held + 1;
    return held >= 2;
  }, `${impl.type} is still pending`);
}

const T = {template: true};
const R = {template: false};

category('Foot-guns: Crux', () => {
  // run alone, the category meets an RDKit module that no earlier category has started
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('rules: a molecule Chem can\'t read back is a block (MOL-01, MOL-02, UI-02)', async () => inCrux((m) => {
    for (const units of ['smiles', 'molblock']) {
      for (const d of ['C(C)(C)(C)(C)C', 'CN(C)(C)C']) {
        const written = units === 'smiles' ? d : molblockOf(m, d);
        expectFlag(`MOL-01 ${units} ${d}`, sketcherWriteFlags(m, written, {units}), 'block', ['valence']);
      }
      const ring = units === 'smiles' ? 'c1cccc1' : molblockOf(m, 'c1cccc1');
      expectFlag(`MOL-02 ${units}`, sketcherWriteFlags(m, ring, {units}), 'block', ['aromatic', 'Kekulé']);
    }
    // pyrrole written without its [nH] (an example of the BDD outline, checked here instead)
    expectFlag('MOL-02 c1ccnc1', sketcherWriteFlags(m, 'c1ccnc1', {units: 'smiles'}), 'block', ['aromatic', 'Kekulé']);
    expectNone('pyrrole', sketcherWriteFlags(m, 'c1cc[nH]c1', {units: 'smiles'}), ['block', 'warn', 'note']);
    expectNone('neopentane', sketcherWriteFlags(m, 'CC(C)(C)C', {units: 'smiles'}), ['block', 'warn', 'note']);
    expectFlag('UI-02 every fragment named', sketcherWriteFlags(m,
      'C(C)(C)(C)(C)C.N(C)(C)(C)(C)C.c1cccc1.[CH3:1]C(C)(C)(C)(C)C', {units: 'smiles'}), 'block',
    ['valence on C (atom 1)', 'N (atom 7)', 'Kekulé']);
    expectFlag('MOL-01 cell', moleculeFlags(m, 'C(C)(C)(C)(C)C', {units: 'smiles'}), 'block', ['valence']);
  }));

  test('rules: a query in a SMILES column is a block, in a molblock column a warn (MOL-03, MOL-04)', async () =>
    inCrux((m) => {
      for (const q of ['[C,N]1=CC=CC=C1', 'C~C', '[#6;R]1CCCCC1']) {
        const drawing = molblockOf(m, q, true);
        expectFlag(`MOL-03 ${q} (Crux's sketcher says query)`, sketcherWriteFlags(m, '', {units: 'smiles'}, '',
          drawing),
        'block', ['query', 'SMILES']);
        expectFlag(`MOL-03 ${q} (another sketcher's empty SMILES)`,
          sketcherWriteFlags(m, '', {units: 'smiles'}, null, drawing), 'block', ['query', 'SMILES']);
      }
      const q = molblockOf(m, '[C,N]1=CC=CC=C1', true);
      expectNone('MOL-03 molblock column', sketcherWriteFlags(m, q, {units: 'molblock'}, null, q), ['block']);
      expectFlag('MOL-04 footer', sketcherWriteFlags(m, q, {units: 'molblock'}, null, q, molblockOf(m, 'c1ccccc1')),
        'warn', ['query']);
      expectFlag('MOL-04 cell', moleculeFlags(m, q, {units: 'molblock'}), 'warn', ['query']);
      expectNone('an empty drawing', sketcherWriteFlags(m, '', {units: 'smiles'}, null, ''), ['block', 'warn']);
    }));

  // The product owner (2026-11-19): opening the molecule editor on row 2 or 3 of mol1K.sdf (a molblock column) warned
  // "A query (an atom query) in a molecule column". The check reads a molblock as SMARTS, and SMARTS writes a
  // molecule's own stereo centres, charges, H counts and isotopes as `[#6@H]`, `[#7+]`, `[13#6]`: molecule features.
  test('rules: a molecule\'s stereo, charges, H counts and isotopes in a molblock column are no query (MOL-04)',
    async () => inCrux(async (m) => {
      const M = {units: 'molblock'};
      const noQuery = (name: string, flags: Flag[]) => expect(flags.some((f) => f.rule === 'queryInMoleculeColumn'),
        false, `${name}: no query warn; said: ${said(flags)}`);
      const molecules = [
        '[H][C@]1(O[C@H]2C[C@H](N)[C@H](O)[C@H](C)O2)C[C@](O)(C(=O)CO)Cc2c(O)c3c(c(O)c21)C(=O)c1c(OC)cccc1C3=O',
        'Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](COP(=O)(O)O)O2)c(=O)[nH]c1=O',
        'C[13CH2]O', '[2H]OC', 'C[N+](C)(C)[O-]', '[NH4+].[Cl-]', 'C[S@@](=O)C', '[Fe+2].[Cl-].[Cl-]',
        'C[C@](F)(Cl)Br |&1:1|', 'c1cc[nH]c1',
      ];
      for (const s of molecules) {
        const mol = m.get_mol(s);
        const blocks = [mol.get_molblock(), mol.get_v3Kmolblock()];
        mol.delete();
        for (const [i, block] of blocks.entries()) {
          const name = `${s} (${i ? 'V3000' : 'V2000'})`;
          noQuery(`${name} cell`, moleculeFlags(m, block, M));
          noQuery(`${name} editor opened`, sketcherWriteFlags(m, block, M, null, block, block));
        }
      }
      // the file itself, as the editor and the Checks section see each of its cells
      const [df] = _importSdf(await loadFileAsBytes('mol1K.sdf'));
      const col = df.columns.byName('molecule');
      expect(col.length, 1000, 'mol1K.sdf read as 1,000 molblocks');
      const flagged = col.toList().map((v, i) => [i + 1, v] as [number, string])
        .filter(([, v]) => moleculeFlags(m, v, M).some((f) => f.rule === 'queryInMoleculeColumn')).map(([i]) => i);
      expect(flagged.length, 0, `mol1K.sdf rows flagged as queries: ${flagged.slice(0, 20).join(', ')}`);
      // a query drawn into the column is still one
      expectFlag('MOL-04 a charged query atom', moleculeFlags(m, molblockOf(m, '[N+;H0](C)(C)(C)C', true), M),
        'warn', ['query']);
    }));

  test('rules: what an edit adds that a canonical SMILES misreads (MOL-05, 06, 08, 10..13)', async () => inCrux((m) => {
    const S = {units: 'smiles'};
    const M = {units: 'molblock'};
    expectFlag('MOL-05', sketcherWriteFlags(m, '[CH3:1][OH:2]', S, null, '', 'CO'), 'warn', ['map', 'duplicate']);
    expectFlag('MOL-06 R1', sketcherWriteFlags(m, '[*:1]c1ccccc1', S, null, '', 'c1ccccc1'), 'warn', ['R1']);
    expectFlag('MOL-06 AP', sketcherWriteFlags(m, '*c1ccccc1 |$_AP1;;;;;;$|', S, null, '', 'c1ccccc1'), 'warn',
      ['attachment']);
    expectNone('MOL-06 a core cell', moleculeFlags(m, '[*:1]c1ccc([*:2])cc1', S), ['block', 'warn']);
    expectFlag('MOL-08', moleculeFlags(m, 'C[C@H](O)[C@@H](C)F |o1:1,3|', S), 'note', ['enhanced stereo', 'canonical']);
    expectNone('MOL-08 no warn', moleculeFlags(m, 'C[C@H](O)[C@@H](C)F |o1:1,3|', S), ['block', 'warn']);
    expectFlag('MOL-09', moleculeFlags(m, 'CCO |Sg:n:1,2::ht|', S), 'note', ['S-groups', 'canonical']);
    expectNone('MOL-09 no warn', moleculeFlags(m, 'CCO |Sg:n:1,2::ht|', S), ['block', 'warn']);
    const lines = molblockOf(m, 'CCO').split('\n');
    for (let i = 4; i < 7; i++)
      lines[i] = lines[i].slice(0, 20) + (0.5 * (i - 3)).toFixed(4).padStart(10) + lines[i].slice(30);
    const threeD = lines.join('\n');
    expectNone('MOL-10 unedited', sketcherWriteFlags(m, threeD, M, null, threeD, threeD), ['block', 'warn', 'note']);
    expectFlag('MOL-10 edited flat', sketcherWriteFlags(m, molblockOf(m, 'CCN'), M, null, '', threeD), 'warn', ['3D']);
    expectFlag('MOL-11', sketcherWriteFlags(m, molblockOf(m, 'C'.repeat(1200)), M, null, '', molblockOf(m, 'CCO')),
      'warn', ['999', 'V3000']);
    expectNone('MOL-12', sketcherWriteFlags(m, '[H]OC([H])([H])[H]', S, null, '', 'CO'), ['block', 'warn']);
    expectFlag('MOL-13', sketcherWriteFlags(m, 'CCO.c1ccccc1', S, null, '', 'CCO'), 'note', ['2', 'fragments']);
  }));

  test('rules: templates (TPL-01..17, TRF-01)', async () => inCrux((m) => {
    const tpl = (t: string) => reactionFlags(m, t, T);
    expectFlag('TPL-01', tpl('[C:1](=[O:2])[OH].[NH2:4]>>[C:1](=[O:2])[NH:5]'), 'warn', ['map number 5',
      'product side']);
    expectFlag('TPL-02', tpl('[c:1]1[c:2][c:3][c:4][c:5][c:6]1[CH2][O]C(=O)[NH:7]>>[NH2:7]'), 'note', ['removed']);
    expectFlag('TPL-03', tpl('[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]'), 'warn', ['map number 1', 'twice']);
    expectFlag('TPL-04', tpl('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3][C:1]'), 'warn', ['map number 1', 'twice']);
    expectFlag('TPL-05', tpl('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[N:3]C'), 'note', ['no map number', 'new']);
    expectFlag('TPL-06 chain', tpl('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C'), 'warn',
      ['route', '2']);
    expectFlag('TPL-07', tpl('[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]'), 'block', ['placeholder', '${n}']);
    expectFlag('TPL-09', tpl('[c:1][Cl,Br:2]>>[c:1][F,I:2]'), 'warn', ['product', 'ignored']);
    expectFlag('TPL-10', tpl('[NH2:1].[O]=[C]=[O].[NH2:2]>>[NH:1]C(=O)[NH:2]'), 'warn', ['3 reactants',
      'no mapped atom']);
    expectFlag('TPL-11', tpl('[C:1](=[O:2])[OH].[NH2:3]>CN(C)C=O>[C:1](=[O:2])[NH:3]'), 'warn', ['agents',
      'enumeration']);
    expectFlag('TPL-14', tpl('[C:1](=[O:2])[OH]>>'), 'warn', ['no', 'product']);
    // the inputs of templates.feature's outline that the lines above lack (the BDD keeps one example a level; these
    // checks live here)
    expectFlag('TPL-06 a --**-- route', tpl('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]--**--' +
      '[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C'), 'warn', ['route', '2']);
    expectFlag('TPL-07 ${rGroupNumber}', tpl('[#8]=[#6](-[#7:1])-[#8]-[#6]-[#6]1-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-' +
      '[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2-1>>[*:1][*${rGroupNumber}]'), 'block', ['placeholder', '${rGroupNumber}']);
    expectFlag('TPL-09 a recursive SMARTS', tpl('[C:1][OH:2]>>[C:1][O;$(OC):2]C'), 'warn', ['product', 'ignored']);
    expectFlag('TPL-10 an agent', tpl('[C:1](=[O:2])[OH].[NH2:3].CN(C)C=O>>[C:1](=[O:2])[NH:3]'), 'warn',
      ['3 reactants', 'no mapped atom']);
    expectFlag('TPL-14 no reactant', tpl('>>[C:1](=[O:2])O'), 'warn', ['no', 'reactant']);
    expectFlag('TPL-15', tpl('[C;$(C[$(C[$(C[$(C[$(C[$(C[$(C[$(C[$(C)])])])])])])])]):1]>>[C:1]O'), 'block',
      ['recursive', '8']);
    expectFlag('TPL-17', tpl('CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|'), 'block', ['f:']);
    expectNone('TPL-17 fixed', tpl('CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:2.3|'), ['block']);
    for (const t of ['[c:1][Br].[#5]([OH])([OH])[c:2]>>[c:1][c:2]',
      '[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]',
      '[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]'])
      expectNone(`clean ${t}`, tpl(t), ['block', 'warn']);
    expectFlag('TRF-01', reactantCountFlags('[NH2:1].[O]=[C]=[O].[NH2:2]>>[NH:1]C(=O)[NH:2]', 2), 'warn',
      ['3 reactants', '2']);
    expectNone('TRF-01 two', reactantCountFlags('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]', 2), ['warn']);
    expectNone('a Transform runs agents', checkTransform('[C:1](=[O:2])[OH].[NH2:3]>CN(C)C=O>[C:1](=[O:2])[NH:3]', 2),
      ['warn']);
  }));

  test('rules: reactions (RXN-01..04) and mappings (MAP-01, MAP-02)', async () => inCrux((m) => {
    const rx = (t: string) => reactionFlags(m, t, R);
    expectFlag('RXN-01', rx('CN(C)C.CI>>C[N+](C)(C)C.[I-] |f:7.9|'), 'block', ['f:', 'Chem reads no reaction']);
    expectFlag('RXN-02', rx('CC>O>CO>N'), 'block', ['no reaction']);
    expectFlag('RXN-03 valence', rx('C(C)(C)(C)(C)C>>CO'), 'warn', ['reactant 1', 'valence']);
    expectFlag('RXN-03 Kekulé', rx('c1cccc1>>C1CCCC1'), 'warn', ['reactant 1', 'Kekulé']);
    expectFlag('RXN-04', rx('>CN(C)C=O>'), 'warn', ['no reactant']);
    expectNone('a clean reaction', rx('CCO.CC(=O)O>>CCOC(C)=O.O'), ['block', 'warn', 'note']);
    expectNone('a route in a reaction column', rx('CCO>>CC=O>>CC(=O)O'), ['block', 'warn', 'note']);
    const chain = 'CC(=O)O.NC>CN(C)C=O>CC(=O)NC>>CC(=O)N(C)C';
    expectFlag('RXN-05', rx(chain), 'note', ['agents', 'chain']);
    expectNone('RXN-05 no warn', rx(chain), ['block', 'warn']);
    expectFlag('MAP-01 under 0.5', mappingFlags(0.3), 'warn', ['check this mapping']);
    expectNone('MAP-01 over 0.5', mappingFlags(0.94), ['warn']);
    for (const t of ['[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]',
      '[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]'])
      expectFlag(`MAP-02 ${t}`, checkTemplateMapping(t), 'warn', ['template', 'map numbers']);
    expectNone('MAP-02 a reaction', checkTemplateMapping('CCO.CC(=O)O>>CCOC(C)=O.O'), ['warn']);
  }));

  test('rules: a template in a column taken for something else is checked as one (TPL-13, RXN-12); MOL-17',
    async () => inCrux((m) => {
      const route = {template: false, units: 'route'};
      const smiles = {template: false, units: 'smiles'};
      expectFlag('RXN-12 the route', reactionFlags(m,
        '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C', route), 'warn', ['route']);
      expectFlag('RXN-12 a template', reactionFlags(m, '[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]', route), 'warn',
        ['map number 1', 'twice']);
      expectFlag('TPL-13 a template', reactionFlags(m, '[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]', smiles),
        'warn',
        ['map number 1', 'twice']);
      const vinyl = reactionFlags(m, '[C:1]=[C:2].[c:3][Br]>>[C:1]=[C:2][c:3]', smiles);
      expectFlag('TPL-13 the note', vinyl, 'note', ['SMILES']);
      expectNone('TPL-13 a clean template', vinyl, ['block', 'warn']);
      const mapped = '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[NH:5][CH3:6]';
      expectNone('a mapper\'s mapped reaction is no template', reactionFlags(m, mapped, smiles), ['block', 'warn',
        'note']);
      expectFlag('MOL-17', checkReadLosses(['text', 'arrow']), 'warn', ['arrow']);
      expectNone('MOL-17 a paste that kept everything', checkReadLosses([]), ['warn']);
    }));

  test('rules: the reaction cell editor\'s OK (RXN-15, RXN-16, RXN-04, TPL-12, TPL-03)', async () => inCrux(() => {
    const S = {template: false, units: 'smiles'};
    const T = {template: true, units: 'smarts'};
    expectFlag('RXN-15 space', checkReactionEdit('CC >> CO', 'CC >> CO', S, false), 'block', ['reaction', 'space']);
    const dative = '[N:1]->[B:2]>>[N:1].[B:2]';
    expectFlag('RXN-15 dative', checkReactionEdit(dative, dative, T, false), 'block', ['reaction', '->']);
    const amide = 'CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1';
    expectNone('opened, unedited', checkReactionEdit(amide, amide, S, false), ['block', 'warn', 'note']);
    expectFlag('RXN-16', checkReactionEdit(amide, 'CC(=O)O.NCc1ccccc1.CC(=O)NCc1ccccc1', S, true), 'block', ['arrow']);
    expectFlag('RXN-04', checkReactionEdit('CC(=O)O.NC>>CC(=O)NC', 'CC(=O)O.NC>>', S, true), 'warn', ['no product']);
    const tpl = '[NX3;H2,H1;!$(N=*);!$(N[#6]=[O,S,N]):1].[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]';
    expectFlag('TPL-12', checkReactionEdit(tpl, '[N:1](C).[#6:2][C:3](=[O:4])[OH,O-]>>[N:1][C:3](=[O:4])[#6:2]', T,
      true),
    'warn', ['H2,H1', '!$(N=*)', '!$(N[#6]=[O,S,N])']);
    expectNone('TPL-12 a component left as written', checkReactionEdit(tpl, `${tpl}C`, T, true)
      .filter((f) => f.rule === 'queryLost'), ['warn']);
    expectFlag('TPL-03', checkReactionEdit('[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]',
      '[C:1](=[O:2])[OH].[N;H2:1]>>[C:1](=[O:2])[N:3]', T, true), 'warn', ['map number 1', 'twice']);
    // an RXN block (the rule reads its header: the editor's original was one): methane to methanol
    const mol = (symbol: string) => ['$MOL', '', '     RDKit          2D', '',
      '  1  0  0  0  0  0  0  0  0  0999 V2000',
      `    0.0000    0.0000    0.0000 ${symbol.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`, 'M  END'].join('\n');
    const block = ['$RXN', '', '      RDKit', '', '  1  1', mol('C'), mol('O')].join('\n');
    const rxnColumn = {template: false, units: 'rxn'};
    expectFlag('RXN-17', checkReactionEdit(block, 'CC(=O)O.CNC>>CC(=O)NC', rxnColumn, true), 'warn',
      ['RXN', 'reaction SMILES']);
    expectNone('RXN-17 unedited', checkReactionEdit(block, block, rxnColumn, false), ['warn']);
  }));

  test('rules: the enumerator\'s filters (TRF-05, TRF-06 block; TRF-07 warns)', async () => inCrux((m) => {
    const rdkit = getRdKitModule();
    expectFlag('TRF-05', enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: ['[N;H2']}],
      buildingBlocks: ['CC(=O)O', 'NCCN', 'CCN'], exclusionSmarts: []}), 'block', ['[N;H2', 'cannot be read']);
    expectNone('TRF-05 fixed', enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: ['[N;H2]']}],
      buildingBlocks: ['CC(=O)O', 'NCCN', 'CCN'], exclusionSmarts: []}), ['block', 'warn']);
    // the enumerator cuts a blocking cell at ';' and '|' (shared.ts): a SMARTS holding ';' comes in pieces, named whole
    expectFlag('TRF-05 cut, unreadable', enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: ['[N', 'H2']}],
      buildingBlocks: ['CC(=O)O', 'NCCN', 'CCN'], exclusionSmarts: []}), 'block', ['[N;H2', 'cannot be read']);
    expectFlag('TRF-05 cut, readable whole', enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: ['[N',
      'H2]']}],
    buildingBlocks: ['CC(=O)O', 'NCCN', 'CCN'], exclusionSmarts: []}), 'block', ['[N;H2]', 'cut at ";"', '&']);
    const exclusion = enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: []}],
      buildingBlocks: ['CC(=O)O', 'CCN'], exclusionSmarts: ['C(=O)Cl]', 'S(=O)(=O)[#6]']});
    expectFlag('TRF-06', exclusion, 'block', ['C(=O)Cl]', 'cannot be read']);
    expect(exclusion.some((f) => f.text.includes('S(=O)(=O)[#6]')), false, 'TRF-06: a readable SMARTS is not named');
    expectFlag('TRF-07', enumeratorFlags(m, rdkit, {templates: [{blockingSmartsList: ['[NX3](=O)=O']}],
      buildingBlocks: ['O=C(O)c1ccc([N+](=O)[O-])cc1', 'CCN'], exclusionSmarts: []}), 'warn',
    ['[NX3](=O)=O', 'matches none of the building blocks']);
  }));

  test('RDKit mode: nothing is flagged, no slot is added (12A)', async () => inRdkit(() => {
    expect(checksModule(), null, 'no checks module in RDKit mode');
    expect(checkSketcherWrite('C(C)(C)(C)(C)C', {units: 'smiles'}, null, '').length, 0, 'the sketcher write');
    expect(checkReaction('[C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]', T).length, 0, 'a template');
    expect(checkValue('C(C)(C)(C)(C)C', DG.SEMTYPE.MOLECULE).length, 0, 'a molecule value');
    expect(checksShown('C(C)(C)(C)(C)C', DG.SEMTYPE.MOLECULE), false, 'the Checks section');
    const dialog = okCancelDialog('RDKit mode');
    try {
      const col = DG.Column.fromStrings('molecule', ['CCO']);
      col.semType = DG.SEMTYPE.MOLECULE;
      const sketcher = fakeSketcher('C(C)(C)(C)(C)C', '');
      expect(watchMoleculeEditor(dialog, sketcher.host, col, 'CCO'), null, 'no slot in RDKit mode');
      expect(dialog.root.querySelector(`.${FOOTER_SLOT_CLASS}`), null, 'no slot in the footer');
      expect(dialog.getButton('OK').disabled, false, 'OK as it was');
    } finally {
      dialog.close();
    }
  }));

  test('footer slot: OK and Cancel stay where they are, a block disables OK, a long text is cut', async () => {
    const dialog = okCancelDialog('Footer slot');
    try {
      await frame();
      const before = buttonRects(dialog);
      const slot = footerSlot(dialog);
      await frame();
      const root = dialog.root.querySelector(`.${FOOTER_SLOT_CLASS}`) as HTMLElement;
      expect(root !== null, true, 'the slot is in the footer');
      const long: Flag[] = [
        {rule: 'moleculeUnreadable', level: 'block',
          text: 'Chem can\'t read this molecule back: an impossible valence ' +
          'on C (atom 1) and on N (atom 7); an aromatic ring with no Kekulé structure. It would be read as a query ' +
          'at best, and every calculation, search and canonical SMILES skips it'},
        {rule: 'mapsInMoleculeColumn', level: 'warn', text: 'Atom map numbers (1) in a molecule column: canonical ' +
          'SMILES keep them, so the molecule is no duplicate of its unmapped self'},
      ];
      const states: [string, Flag[], boolean][] = [['empty', [], false], ['a long block and a warn', long, true],
        ['a warn', [long[1]], false], ['empty again', [], false]];
      for (const [name, flags, blocked] of states) {
        slot.show(flags);
        await frame();
        const now = buttonRects(dialog);
        expect(samePlace(now.ok, before.ok), true, `OK stays put (${name})`);
        expect(samePlace(now.cancel, before.cancel), true, `Cancel stays put (${name})`);
        expect(dialog.getButton('OK').disabled, blocked, `OK ${blocked ? 'disabled' : 'enabled'} (${name})`);
        expect(slot.blocked, blocked, `the slot's block (${name})`);
        expect(root.dataset.level,
          flags.length ? flags.map((f) => f.level).sort()[0] === 'block' ? 'block' : 'warn' : '',
          `the slot's level (${name})`);
        const slotRect = root.getBoundingClientRect();
        expect(slotRect.right <= Math.min(now.ok.left, now.cancel.left) + 0.5 ||
          slotRect.left >= Math.max(now.ok.right, now.cancel.right) - 0.5, true,
        `the slot beside the buttons (${name})`);
      }
      slot.show(long);
      await frame();
      const label = root.lastElementChild as HTMLElement;
      expect(label.scrollWidth > label.clientWidth, true, 'a long text is cut');
      expect((root.getAttribute('aria-label') ?? '').includes('canonical SMILES keep them'), true,
        'the whole text in its label (and its tooltip)');
      expect(root.getBoundingClientRect().height <= 20.5, true, 'one line high');
    } finally {
      dialog.close();
    }
  });

  test('cell editor: a query whose SMILES is empty blocks OK in a SMILES column; a molecule unblocks it (MOL-03)',
    async () => inCrux(async (m) => {
      const dialog = okCancelDialog('Cell editor');
      try {
        const col = DG.Column.fromStrings('molecule', ['CCO']);
        col.semType = DG.SEMTYPE.MOLECULE;
        col.meta.units = DG.chem.Notation.Smiles;
        // another sketcher's way: an empty SMILES of a query drawing (Crux's own says so: CruxSketcher.smilesIsQuery)
        const sketcher = fakeSketcher('', molblockOf(m, '[C,N]1=CC=CC=C1', true));
        const slot = watchMoleculeEditor(dialog, sketcher.host, col, 'CCO');
        expect(slot !== null, true, 'a slot in Crux mode');
        await frame();
        expect(dialog.getButton('OK').disabled, true, `OK disabled: ${said([...slot!.shown])}`);
        expectFlag('the slot', [...slot!.shown], 'block', ['query', 'SMILES']);
        sketcher.draw('CCN', molblockOf(m, 'CCN'));
        expect(dialog.getButton('OK').disabled, false, 'OK enabled again');
        expect(slot!.shown.length, 0, 'the slot empty');
      } finally {
        dialog.close();
      }
    }));

  test('cell editor: an edit of only enhanced stereo, an attachment point or an S-group is written (MOL-07)',
    async () => {
      const pairs: [string, string][] = [
        ['C[C@H](O)[C@@H](C)F', 'C[C@H](O)[C@@H](C)F |o1:1,3|'],
        ['C[C@H](O)[C@@H](C)F', 'C[C@H](O)[C@@H](C)F |&1:1,3|'],
        ['*c1ccccc1', '*c1ccccc1 |$_AP1;;;;;;$|'],
        ['CCO', 'CCO |Sg:n:1,2::ht|'],
      ];
      await inCrux((m) => {
        // what the comparison rests on: Crux keeps an S-group in Chem's stereo setting (chem-wasm 732af7a03; before
        // it, none was read, and the S-groups' text was compared)
        const sg = m.get_mol('CCO |Sg:n:1,2::ht|');
        try {
          expect(sg.get_cxsmiles(JSON.stringify({CX_ALL_BUT_COORDS: true})).includes('Sg:n:'), true,
            'Crux keeps the S-group in its canonical CXSMILES');
        } finally {
          sg.delete();
        }
        for (const [before, after] of pairs)
          expect(unchangedByEdit(before, after), false, `Crux mode: ${before} -> ${after} is an edit`);
        expect(unchangedByEdit('OCC', 'CCO'), true, 'Crux mode: the same molecule redrawn is no edit');
        expect(unchangedByEdit('C[C@H](O)[C@@H](C)F |o1:1,3|', 'C[C@H](O)[C@@H](C)F |o1:1,3|'), true,
          'Crux mode: the same group kept is no edit');
      });
      await inRdkit(() => {
        // RDKit mode as it always was: canonical SMILES, which drop the stereo group
        expect(unchangedByEdit(pairs[0][0], pairs[0][1]), true, 'RDKit mode compares canonical SMILES, as before');
      });
    });

  // crux-sketch refuses a ChemDraw reaction in a molecule document with words, so MOL-17's silent loss does not
  // happen; this pins the refusal, so a crux-sketch that came to drop the arrow silently fails here. A BDD scenario
  // cannot deliver a paste into the sketcher's shadow root, so the same file goes through addMolecule, which reads and
  // places a text as a paste does (crux-sketch API-025).
  test('cell editor: a ChemDraw reaction is refused by the Crux sketcher with words, the cell\'s molecule kept ' +
    '(MOL-17)',
  async () => inCrux(async () => {
    // RDKit's test data, as Crux Sketch keeps it (its notice: NOTICE.chemdraw_template1 beside this file)
    const cdxml = await _package.files.readAsText('tests/chemdraw_template1.cdxml');
    expect(/<arrow\b/.test(cdxml), true, 'the file holds a reaction arrow');
    const sketcherType = grok.chem.currentSketcherType;
    grok.chem.currentSketcherType = 'Crux';
    const host = new DG.chem.Sketcher();
    host.setSmiles('CCO'); // the cell's value, as the cell editor opens on it
    const dialog = ui.dialog('Cell editor').add(host).onOK(() => {}).show();
    try {
      await host.sketcherReady();
      const impl = host.sketcher!;
      expect(impl.root.classList.contains('crux-sketcher'), true, 'the sketcher is Crux');
      const crux = impl.root.querySelector<CruxAdd>('crux-sketch');
      expect(crux !== null, true, 'Crux\'s element is there');
      await settled(impl);
      const cell = host.getSmiles();
      expect(cell, 'CCO', 'the cell\'s molecule drawn');
      let changes = 0;
      const sub = impl.onChanged.subscribe(() => changes++);
      try {
          crux!.addMolecule(cdxml);
          // CDXML is read once crux-sketch's CDXML module has loaded
          await waitUntil(() => crux!.lastError !== null || changes > 0,
            () => `neither refused nor added; the drawing: ${host.getSmiles()}`, 30000);
          await settled(impl);
          const error = JSON.stringify(crux!.lastError);
          expect(changes, 0,
            `no change event (a reaction read as one molecule is one); the drawing: ${host.getSmiles()}`);
          expect(crux!.lastError?.format, 'cdxml', `refused as CDXML: ${error}`);
          expect(/Reactions are not supported/.test(crux!.lastError?.message ?? ''), true,
            `the refusal's words: ${error}`);
          expect(host.getSmiles(), cell, 'what OK writes into the cell is still the cell\'s molecule');
          const losses = (impl as unknown as {readLosses: string[]}).readLosses;
          expect(losses.length, 0, `nothing read and left out: ${losses.join(', ')}`);
      } finally {
        sub.unsubscribe();
      }
    } finally {
      dialog.close();
      grok.chem.currentSketcherType = sketcherType;
    }
  }));

  // RXN-06: the reaction cell editor says in place of its sketcher why Crux Sketch did not open the cell's reaction,
  // and writes nothing on OK; the slot keeps OK disabled, saying nothing more (the refusal is said once, not twice).
  // The editor here is a stand-in with the editor's fields: the slot reads `refused` once the editor is ready and at
  // each change.
  test('reaction cell editor: a reaction Crux Sketch does not open keeps OK off, the slot blank (RXN-06)',
    async () => inCrux(async () => {
      const route = 'CC(=O)O.NC>>CC(=O)NC--**--c1ccccc1Br>>c1ccccc1B(O)O--**--CC(=O)NC.c1ccccc1B(O)O>>CC(=O)NCc1ccccc1';
      const col = DG.Column.fromStrings('reaction', [route]);
      col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
      col.meta.units = 'route';
      // what the slot said before it read `refused`: the value's text alone, unedited, which raises nothing
      expectNone('the text alone', checkReactionEdit(route, route, {template: false, units: 'route'}, false),
        ['block', 'warn']);
      const dialog = okCancelDialog('Reaction');
      try {
        const onChange = new Subject<string>();
        const editor = {dialog, cell: {cell: {column: col, value: route}}, original: route, edited: false, onChange,
          ready: Promise.resolve(), value: route, slot: null,
          refused: 'a route that is not a chain of steps opens in a later version'};
        const slot = watchReactionEditor(editor as unknown as ReactionCellEditor);
        expect(slot !== null, true, 'a slot in Crux mode');
        await editor.ready;
        await frame();
        // the editor says why once, in place of its sketcher; the slot does not say it again
        expectNone('RXN-06: the slot', [...slot!.shown], ['block', 'warn']);
        expect(slot!.root.textContent, '', 'the slot blank: the refusal is said once, by the editor');
        expect(dialog.getButton('OK').disabled, true, `OK disabled: ${said([...slot!.shown])}`);
        // something drawn after it changes nothing: OK writes nothing while the editor holds a refusal
        editor.edited = true;
        onChange.next('CC>>CO');
        expect(slot!.root.textContent, '', 'still blank after an edit');
        expect(dialog.getButton('OK').disabled, true, 'OK still disabled');
      } finally {
        dialog.close();
      }
    }));

  // A reagent name drawn as an atom label ("HATU") used to come back as "*" without a word: the slot names what OK's
  // text does not keep of the drawing, from Crux Sketch's warnings (reactions/reaction-losses.ts), as warns, OK
  // enabled; an RXN block keeps a label; molecule losses are left to their own rules. The editor is a stand-in with
  // the editor's fields, as RXN-06's.
  test('reaction cell editor: a label written as "*" and the drawing\'s texts and arrow styles warn in the slot; OK ' +
    'stays enabled (16A)', async () => inCrux(async () => {
    const row = (format: string, code: string, atoms: number[] = [], message = code) =>
      ({kind: 'format' as const, format, code, message, atoms});
    const warnings = [row('smiles', 'alias', [4]), row('smarts', 'alias', [4]), row('molV3000', 'alias', [4]),
      row('smiles', 'texts'), row('smarts', 'texts'), row('smiles', 'arrow-style'), row('smiles', 'stereo-groups', [1]),
      row('smiles', 'shapes', [], 'the shapes')];
    const smiles = 'CC(=O)O.CCN.*>>CCNC(C)=O';
    const losses = drawingLosses(warnings, smiles);
    expect(JSON.stringify(losses), JSON.stringify({format: 'smiles', labels: 1, objects: [
      'the texts (conditions and captions)', 'the arrows\' styles (each is written as a plain arrow)',
      'the shapes']}), 'reaction SMILES: the label, the objects in Chem\'s words, a later one in its own');
    const rxn = drawingLosses(warnings, '$RXN\n\n      Crux\n\n  1  1\n');
    expect(rxn.labels, 0, 'an RXN block keeps a label');
    expect(rxn.objects.length, 3, 'but no reaction text keeps the drawing\'s objects');
    expect(drawingLosses([row('smiles', 'stereo-groups', [1]), row('smiles', 'query', [2])], smiles).objects.length, 0,
      'molecule losses are not the drawing\'s');
    const col = DG.Column.fromStrings('reaction', ['CC(=O)O.CCN>>CCNC(C)=O']);
    col.semType = DG.SEMTYPE.CHEMICAL_REACTION;
    col.meta.units = 'smiles';
    const dialog = okCancelDialog('Reaction');
    try {
      const onChange = new Subject<string>();
      const editor = {dialog, cell: {cell: {column: col, value: col.get(0)}}, original: col.get(0), edited: false,
        onChange, ready: Promise.resolve(), value: col.get(0), slot: null, refused: null, warnings};
      const slot = watchReactionEditor(editor as unknown as ReactionCellEditor);
      await editor.ready;
      await frame();
      expectNone('unedited: OK writes nothing', [...slot!.shown], ['block', 'warn']);
      editor.edited = true;
      onChange.next(smiles);
      expectFlag('the label', [...slot!.shown], 'warn', ['label', '"*"', 'reaction SMILES']);
      expectFlag('the objects', [...slot!.shown], 'warn', ['not kept in reaction SMILES', 'texts', 'styles']);
      expect(dialog.getButton('OK').disabled, false, `OK left clickable, the footer saying: ${said([...slot!.shown])}`);
    } finally {
      dialog.close();
    }
  }));

  // Chem's cell editor as a grid opens it: through the platform's function, which runs in package.js, the bundle whose
  // functions make the sketchers the platform opens. These tests run in package-test.js: its own editMoleculeCell would
  // look for its own CruxSketcher class, which no such sketcher is (the test then saw no block, the product being
  // fine).
  async function editCellAsTheGridDoes(cell: DG.GridCell): Promise<void> {
    const editor = DG.Func.find({package: 'Chem', name: 'editMoleculeCell'})[0];
    expect(editor !== undefined, true, 'the cell editor function of Chem');
    await editor.apply({cell});
  }

  // The engine-poison foot-gun: a Crux getter that throws after an edit (the engine failed on that read, or the
  // drawing cannot be written in the column's notation) answered '', and OK wrote '' over the cell. Chem's own cell
  // editor opens on a real Crux whose handle's four getters are then made to throw (own properties over its class's:
  // no Proxy, the handle unchanged otherwise). The engine recovers by itself (crux-sketch API-069), and a user is never
  // asked to reload, so its failure is a warning and OK stays on, writing nothing; a drawing the notation cannot hold
  // blocks OK, pressed even then (enabled by hand, as Enter could). Either way the cell must keep its value.
  test('cell editor: a Crux getter that throws never empties the cell (an engine failure: a warning, no reload; a ' +
    'drawing not written as SMILES: a block)', async () => inCrux(async () => {
    const engineMessage = 'Crux Sketch: the engine could not do this (panicked at s' +
      'rc/lib.rs:1:1: test); the drawing is kept as it was';
    const cases: [string, () => unknown, 'warn' | 'block', string[]][] = [
      ['a raw RuntimeError', () => new WebAssembly.RuntimeError('unreachable'), 'warn', ['could not be read back']],
      ['crux-sketch\'s engine error', () => Object.freeze({format: 'engine', message: engineMessage}), 'warn',
        ['could not be read back']],
      // the sketcher's own message names the engine; the footer says it in plain words, naming no engine
      ['a drawing not written as SMILES', () => Object.freeze({format: 'smiles',
        message: 'Crux Sketch: the engine could not write this drawing as smiles (test)'}), 'block',
      ['cannot be written as SMILES', 'keeps its value']],
    ];
    const sketcherType = grok.chem.currentSketcherType;
    grok.chem.currentSketcherType = 'Crux';
    try {
      for (const [name, error, level, words] of cases) {
        const col = DG.Column.fromStrings('molecule', ['CCO']);
        col.semType = DG.SEMTYPE.MOLECULE;
        col.meta.units = DG.chem.Notation.Smiles;
        const df = DG.DataFrame.fromColumns([col]);
        const tv = grok.shell.addTableView(df);
        const dialogs = () => [...document.querySelectorAll<HTMLElement>('.d4-dialog')]
          .filter((d) => d.querySelector('.crux-sketcher'));
        const before = dialogs().length;
        let dialog: HTMLElement | null = null;
        let restore: (() => void) | null = null;
        try {
          await editCellAsTheGridDoes(tv.grid.cell('molecule', 0));
          await waitUntil(() => dialogs().length > before, `${name}: the cell editor opened`);
          dialog = dialogs()[dialogs().length - 1];
          const root = dialog.querySelector('[data-u2-name="cruxSketch"]') ?? dialog.querySelector('.crux-sketcher')!;
          const impl = () => DG.Widget.find(root) as unknown as
            {sketch: SketcherHandle | null, isRenderPending: boolean, onChanged: Subject<any>,
              explicitMol: unknown} | null;
          await waitUntil(() => !!impl()?.sketch && !impl()!.isRenderPending, `${name}: Crux ready`, 30000);
          const crux = impl()!;
          restore = throwingGetters(crux.sketch!, error);
          // a person's edit, as Crux's own 'user' change handler makes it: the host's explicit value (the cell's own
          // SMILES, given back as the editor opened) ends, then onChanged. The host relays Crux's changes once its own
          // init has resolved, a moment after Crux is ready (a change told at once is lost), so the edit is
          // told again every 200 ms until the slot answers, for at most 5 s
          const slot = dialog.querySelector<HTMLElement>(`.${FOOTER_SLOT_CLASS}`);
          for (let i = 0; i < 25 && slot?.dataset.level !== level; i++) {
            crux.explicitMol = null;
            crux.onChanged.next(null);
            await new Promise((r) => setTimeout(r, 200));
          }
          const said = slot?.getAttribute('aria-label') ?? '(no slot)';
          expect(slot?.dataset.level, level, `${name}: a ${level} in the footer; it says: ${said}`);
          expect(words.every((w) => said.toLowerCase().includes(w.toLowerCase())), true,
            `${name}: its words name ${JSON.stringify(words)}; it says: ${said}`);
          expect(said.includes('Crux'), false, `${name}: the footer names no engine; it says: ${said}`);
          expect(/reload|restart/i.test(said), false, `${name}: nothing asks for a reload; it says: ${said}`);
          const ok = [...dialog.querySelectorAll<HTMLButtonElement>('.d4-dialog-footer button')]
            .find((b) => /^\s*OK\s*$/i.test(b.textContent ?? ''))!;
          expect(ok.disabled || ok.classList.contains('d4-disabled'), level === 'block',
            `${name}: OK ${level === 'block' ? 'disabled' : 'enabled'}`);
          ok.disabled = false;
          ok.classList.remove('d4-disabled', 'disabled');
          ok.click();
          await frame();
          expect(df.get('molecule', 0), 'CCO', `${name}: the cell keeps its value after OK`);
        } finally {
          // the getters back before the editor closes, if it is still open (the handle itself never changed)
          restore?.();
          restore = null;
          if (dialog?.isConnected) {
            [...dialog.querySelectorAll<HTMLButtonElement>('.d4-dialog-footer button')]
              .find((b) => /cancel/i.test(b.textContent ?? ''))?.click();
          }
          tv.close();
        }
      }
    } finally {
      grok.chem.currentSketcherType = sketcherType;
    }
  }));

  // A sketcher opened on an engine that has failed says so as it starts (lastError {format: 'engine'}, no event a
  // later listener hears) and holds nothing, its getters answering '' without a throw: OK wrote '' after an edit, with
  // no getter throwing. Since crux-sketch API-069 the engine recovers from a trap by itself, and only an engine that
  // can make no fresh instance (its last resort) or a start it refused says so; the guard stays, without a block or a
  // balloon, as an engine failure never asks the user for anything. So the test gives the sketcher such a handle: the
  // real one, its lastError and its four getters shadowed (own properties, no Proxy) before the sketcher sees it. It
  // wraps the loader Chem's sketchers use (given to the hook): this bundle's own copy of Crux Sketch cannot make a
  // sketcher once package.js's has defined <crux-sketch> (the test timed out after earlier tests had opened sketchers).
  test('cell editor: a Crux sketcher started on a failed engine never empties the cell, though no getter throws, ' +
    'and asks for no reload', async () => inCrux(async () => {
    const stopped = Object.freeze({format: 'engine', message: 'Crux Sketch: the engine failed (RuntimeError: ' +
      'unreachable), and no fresh engine could be started'});
    const shadowed: (() => void)[] = [];
    setCruxSketchLoaderForTests(async (load) => {
      const real = await load();
      return {...real, createSketcher: async (...args: Parameters<typeof real.createSketcher>) => {
        const handle = await real.createSketcher(...args);
        const names = ['lastError', 'smiles', 'molfile', 'molV3000', 'smarts'];
        const saved = names.map((n) => [n, Object.getOwnPropertyDescriptor(handle, n)] as const);
        for (const n of names)
          Object.defineProperty(handle, n, {configurable: true, get: () => n === 'lastError' ? stopped : ''});
        shadowed.push(() => {
          for (const [n, d] of saved) {
            if (d)
              Object.defineProperty(handle, n, d);
            else
              delete (handle as unknown as {[k: string]: unknown})[n];
          }
        });
        return handle;
      }} as typeof real;
    });
    const balloons: string[] = [];
    const shell = grok.shell as unknown as {error: (s: unknown) => void};
    const shellError = shell.error;
    shell.error = (s: unknown) => {
      balloons.push(String(s));
    };
    const sketcherType = grok.chem.currentSketcherType;
    grok.chem.currentSketcherType = 'Crux';
    const col = DG.Column.fromStrings('molecule', ['CCO']);
    col.semType = DG.SEMTYPE.MOLECULE;
    col.meta.units = DG.chem.Notation.Smiles;
    const df = DG.DataFrame.fromColumns([col]);
    const tv = grok.shell.addTableView(df);
    const dialogs = () => [...document.querySelectorAll<HTMLElement>('.d4-dialog')]
      .filter((d) => d.querySelector('.crux-sketcher'));
    const before = dialogs().length;
    let dialog: HTMLElement | null = null;
    try {
      await editCellAsTheGridDoes(tv.grid.cell('molecule', 0));
      await waitUntil(() => dialogs().length > before, 'the cell editor opened');
      dialog = dialogs()[dialogs().length - 1];
      const root = dialog.querySelector('[data-u2-name="cruxSketch"]') ?? dialog.querySelector('.crux-sketcher')!;
      const impl = () => DG.Widget.find(root) as unknown as
        {sketch: SketcherHandle | null, isRenderPending: boolean, onChanged: Subject<any>, explicitMol: unknown} | null;
      await waitUntil(() => !!impl()?.sketch && !impl()!.isRenderPending, 'Crux ready', 30000);
      expect(shadowed.length, 1, 'the sketcher was given the stopped engine\'s handle');
      const slot = () => dialog!.querySelector<HTMLElement>(`.${FOOTER_SLOT_CLASS}`);
      // as the editor opens, before any edit: the warning (the host checks the slot once the sketcher is ready)
      await waitUntil(() => slot()?.dataset.level === 'warn',
        () => `a warning as it opens; the slot says: ${slot()?.getAttribute('aria-label') ?? '(no slot)'}`, 5000);
      expect(/could not be read back/i.test(slot()!.getAttribute('aria-label') ?? ''), true,
        'the warning says the drawing could not be read back');
      expect(/reload|restart/i.test(slot()!.getAttribute('aria-label') ?? ''), false, 'it asks for no reload');
      // an edit (as Crux's own 'user' change handler makes it): the getters now answer the empty drawing's ''
      const crux = impl()!;
      crux.explicitMol = null;
      crux.onChanged.next(null);
      await frame();
      expect(slot()?.dataset.level, 'warn', 'still a warning after the edit');
      const ok = [...dialog.querySelectorAll<HTMLButtonElement>('.d4-dialog-footer button')]
        .find((b) => /^\s*OK\s*$/i.test(b.textContent ?? ''))!;
      expect(ok.disabled || ok.classList.contains('d4-disabled'), false, 'OK enabled');
      ok.click();
      await frame();
      expect(df.get('molecule', 0), 'CCO', 'the cell keeps its value after OK');
      expect(balloons.length, 0, `no balloon; shown: ${JSON.stringify(balloons)}`);
    } finally {
      shadowed.forEach((put) => put());
      setCruxSketchLoaderForTests(null);
      shell.error = shellError;
      grok.chem.currentSketcherType = sketcherType;
      if (dialog?.isConnected) {
        [...dialog.querySelectorAll<HTMLButtonElement>('.d4-dialog-footer button')]
          .find((b) => /cancel/i.test(b.textContent ?? ''))?.click();
      }
      tv.close();
    }
  }));

  // An engine failure is no stop (crux-sketch API-069 recovers by itself), so a getter the engine could not answer
  // shows no balloon and asks for no reload; the failure is kept for the hosts that write until the getters answer
  // again
  test('a Crux getter\'s engine failure shows no balloon, and is forgotten once the getters answer', async () =>
    inCrux(async () => {
      const shown: string[] = [];
      const shell = grok.shell as unknown as {error: (s: unknown) => void, warning: (s: unknown) => void};
      const said = {error: shell.error, warning: shell.warning};
      shell.error = (s: unknown) => {
        shown.push(String(s));
      };
      shell.warning = (s: unknown) => {
        shown.push(String(s));
      };
      const refused = Object.freeze({format: 'engine', message: 'Crux Sketch: the engine could not do thi' +
        's (panicked at src/lib.rs:1:1: test); the drawing is kept as it was'});
      try {
        const impl = new CruxSketcher();
        // a stand-in handle (no engine behind it): its four getters throw as a refused read's do, then answer
        const handle = {warnings: []};
        const answer = throwingGetters(handle, () => refused);
        (impl as unknown as {sketch: object}).sketch = handle;
        expect(impl.smiles, '', 'the platform\'s getter answers \'\'');
        expect(impl.molFile, '', 'and again');
        expect(impl.readFailure(['smiles'])?.engine, true, 'the failure is kept for the hosts that write');
        expect(shown.length, 0, `no balloon for two failing reads; shown: ${JSON.stringify(shown)}`);
        answer();
        Object.assign(handle, {smiles: 'CCO', molfile: DG.chem.convert('CCO', DG.chem.Notation.Smiles,
          DG.chem.Notation.MolBlock)});
        expect(impl.smiles, 'CCO', 'the getter answers again');
        expect(impl.molFile.includes('M  END'), true, 'and the molfile');
        expect(impl.readFailure(['smiles']), null, 'nothing kept once the getters answered');
      } finally {
        shell.error = said.error;
        shell.warning = said.warning;
      }
    }));

  test('Checks section: the worst first, each line its level and rule; none for a clean value', async () =>
    inCrux(() => {
      const value = DG.SemanticValue.fromValueType('C(C)(C)(C)(C)C.[CH3:1][OH:2]', DG.SEMTYPE.MOLECULE);
      expect(checksShown(value.value, DG.SEMTYPE.MOLECULE), true, 'shown for a flagged molecule');
      const list = checksPanel(value).root.querySelector(`.${CHECKS_LIST_CLASS}`) ?? checksPanel(value).root;
      const lines = [...list.querySelectorAll('[data-level]')] as HTMLElement[];
      expect(lines.length >= 2, true, `two lines at least: ${lines.map((l) => l.textContent).join(' | ')}`);
      expect(lines[0].dataset.level, 'block', 'the block first');
      expect(lines.some((l) => l.dataset.level === 'warn' && l.dataset.rule === 'mapsInMoleculeColumn'), true,
        'the maps warn');
      expect(checksShown('CCO', DG.SEMTYPE.MOLECULE), false, 'not shown for a clean molecule');
      // (a mapped value whose bracket atoms carry no H count reads as a template: TPL-13; this one is a molecule each)
      expect(checksShown('CC(=O)O.NC>>CC(=O)NC', DG.SEMTYPE.CHEMICAL_REACTION), false,
        'a clean reaction: nothing to say');
    }));

  test('urea formation (dr014) takes two amines and runs in Two-Component', async () => inCrux(async () => {
    const reactions: NamedReaction[] = JSON.parse(await _package.files.readAsText('reactions/default-reactions.json'));
    const urea = reactions.find((r) => r.id === 'urea-formation')!;
    expect(urea !== undefined, true, 'the template is there');
    expectNone('its reactant count', reactantCountFlags(urea.reactionSmarts, 2), ['warn']);
    const result = await runTwoComponentReaction(getRdKitModule(), urea.reactionSmarts, ['CCN'], ['CN'],
      {mode: 'pairwise', removeSaltsAndWater: false});
    expect(result.products.length, 1, 'one pair');
    const rdkit = getRdKitModule();
    const canon = (s: string): string => {
      const mol = rdkit.get_mol(s);
      try {
        return mol.get_smiles();
      } finally {
        mol.delete();
      }
    };
    expect(result.products[0] ? canon(result.products[0]) : '', canon('CCNC(=O)NC'),
      `a urea from ethylamine and methylamine; got ${result.products[0]}`);
  }));

  test('a variable used twice is filled in twice in Crux mode, once in RDKit mode (TPL-08)', async () => {
    const smarts = '[c:1][Br]>>[c:1][*${n}].[*${n}]';
    await inCrux(() => expect(resolveReactionVariables(smarts, {n: 1}), '[c:1][Br]>>[c:1][*1].[*1]', 'Crux mode'));
    await inRdkit(() => expect(resolveReactionVariables(smarts, {n: 1}), '[c:1][Br]>>[c:1][*1].[*${n}]',
      'RDKit mode as it was'));
  });

  test('Reaction panel: a route\'s steps are its own text; a restated intermediate is no step (RXN-13)', async () =>
    inCrux(() => {
      const steps = (reaction: string): string[] => {
        const widget = reactionPanel(reaction);
        try {
          return ([...widget.root.querySelectorAll('canvas[data-reaction]')] as HTMLCanvasElement[]).slice(1)
            .map((c) => c.dataset.reaction!);
        } finally {
          widget.detach();
        }
      };
      const a = '[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]';
      const b = '[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C';
      expect(steps(`${a}--**--${b}`).join(' | '), `${a} | ${b}`, 'a template route: its own steps, H counts kept');
      const s1 = 'Oc1ccccc1.O=C=O>[Na+].[OH-].O=S(=O)(O)O>O=C(O)c1ccccc1O';
      const s2 = 'O=C(O)c1ccccc1O.CC(=O)OC(C)=O>O=P(O)(O)O>CC(=O)Oc1ccccc1C(=O)O.CC(=O)O';
      expect(steps(`${s1}--**--${s2}`).join(' | '), `${s1} | ${s2}`, 'a specific route: agents and all');
      expect(steps('CCO.CC(=O)O>>CCOC(C)=O>>CCOC(C)=O.O>>CCO').join(' | '),
        'CCO.CC(=O)O>>CCOC(C)=O | CCOC(C)=O.O>>CCO',
        'a restated intermediate (A.B>>C>>C.D>>E) is two steps, as Crux reads it, each as the value writes it');
    }));

  // The panel maps on its own: a template, or a reaction with a component no molecule reader reads, is refused as the
  // panel opens, with the checker's flags, and the mapper is never asked. A stand-in mapper
  // counts the asks (with it set the panel has a Mapping section whether or not Chem's Crux file has the mapper).
  test('Reaction panel: the mapping refuses a template, saying it is one, and names an unreadable component, never ' +
    'asking the mapper (MAP-02, MAP-03)', async () => inCrux(async () => {
    let asked = 0;
    setReactionMapperForTests(async () => {
      asked++;
      throw new Error('the mapper was asked');
    });
    const run = async (reaction: string): Promise<{state: string, flags: string[], said: string}> => {
      const widget = reactionPanel(reaction);
      try {
        const status = widget.root.querySelector<HTMLElement>(`.${MAPPING_CLASS}`);
        const flags = ([...widget.root.querySelectorAll('[data-level]')] as HTMLElement[])
          .map((l) => `${l.dataset.level}: ${l.textContent}`);
        // past the panel's delay before it maps (300 ms), its widget still attached: a mapping it started would ask
        await new Promise((r) => setTimeout(r, 600));
        return {state: status?.dataset.state ?? '(no Mapping section)', flags, said: widget.root.textContent ?? ''};
      } finally {
        widget.detach();
      }
    };
    try {
      const t = await run('[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]');
      expect(t.state, 'refused', `MAP-02: the template refused; said: ${t.said}`);
      expect(t.flags.some((f) => f.startsWith('warn:') && /template/i.test(f) && /map numbers/i.test(f)), true,
        `MAP-02: a warn as the panel opens; flags: ${t.flags.join(' | ')}`);
      expect(t.said.includes('is not a valid SMILES string'), false, `MAP-02: no SMILES error; said: ${t.said}`);
      expect(/not mapped/i.test(t.said), true, `MAP-02: nothing mapped; said: ${t.said}`);
      const u = await run('C(C)(C)(C)(C)C>>CO');
      expect(u.state, 'refused', `MAP-03: refused; said: ${u.said}`);
      expect(/reactant 1/i.test(u.said), true, `MAP-03: reactant 1 named; said: ${u.said}`);
      expect(u.said.includes('is not a valid SMILES string'), false, `MAP-03: no SMILES error; said: ${u.said}`);
      expect(asked, 0, 'the mapper was never asked');
    } finally {
      setReactionMapperForTests(null);
    }
  }));

  test('enumerator: a route or agents over the arrow is skipped saying so (TPL-06, TPL-11); Chem\'s demo exclusion ' +
    'file excludes a nitro product, and an exclusion SMARTS that matched nothing is named (TRF-07)', async () =>
    inCrux(async () => {
      const config = cloneConfig(DEFAULT_CONFIG);
      config.enumeration.num_rounds = 1;
      config.products_specs.min_num_carbon_atoms = -1;
      const rdkit = getRdKitModule();
      const canon = (s: string): string => {
        const mol = rdkit.get_mol(s);
        try {
          return mol.get_smiles();
        } finally {
          mol.delete();
        }
      };
      const template = (smarts: string, name: string) => ({smarts, blockingSmartsList: [], reactionName: name});
      const route = await enumerate({rdkit, config, buildingBlocks: ['CC(=O)O', 'CCN'], exclusionSmarts: [],
        templates: [template('[C:1](=[O:2])[OH].[NH2:3]>>[C:1](=[O:2])[NH:3]>>[C:1](=[O:2])[N:3]C', 'route'),
          template('[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]', 'amide')]});
      expect(route.rows.some((r) => r.product === canon('CCNC(C)=O')), true,
        `TPL-06: the amide made; rows: ${route.rows.map((r) => r.product).join(', ')}`);
      expect(route.warnings.some((w) => w.includes('Template 1') && /route/i.test(w)), true,
        `TPL-06: template 1 named as a route; warnings: ${route.warnings.join(' | ')}`);
      const agents = await enumerate({rdkit, config, buildingBlocks: ['CC(=O)O', 'CCN'], exclusionSmarts: [],
        templates: [template('[C:1](=[O:2])[OH].[NH2:3]>CN(C)C=O>[C:1](=[O:2])[NH:3]', 'agents')]});
      expect(agents.warnings.some((w) => w.includes('Template 1') && /agents over the arrow/i.test(w)), true,
        `TPL-11: template 1 named for its agents; warnings: ${agents.warnings.join(' | ')}`);
      // TRF-07, a data fix in both modes: the demo file's nitro line now matches a nitro as both engines hold it
      // ([N+](=O)[O-]); its old line [NX3](=O)=O matched none, so it excluded nothing, and the run says so
      const exclusions = DG.DataFrame.fromCsv(await _package.files.readAsText('enumerations/ex_smarts.csv'))
        .col('SMARTS')!.toList().map((s) => String(s ?? '').trim()).filter((s) => s);
      // "Reject charged species" (on by default) drops a nitro product whatever the exclusion list says (on, the run
      // gives no rows at all): off here, so the list decides
      const nitroConfig = cloneConfig(config);
      nitroConfig.products_specs.remove_charged_species = false;
      const nitroRun = (exclusionSmarts: string[]) => enumerate({rdkit, config: nitroConfig, exclusionSmarts,
        buildingBlocks: ['O=C(O)c1ccc([N+](=O)[O-])cc1', 'CCN'],
        templates: [template('[C:1](=[O:2])[OH].[N;H2:3]>>[C:1](=[O:2])[N:3]', 'amide')]});
      const nitroProduct = canon('CCNC(=O)c1ccc([N+](=O)[O-])cc1');
      const fixed = await nitroRun(exclusions);
      expect(fixed.rows.some((r) => r.product === nitroProduct), false,
        `TRF-07: the demo file (${exclusions.join(' ; ')}) excludes the nitro product; rows: ` +
      fixed.rows.map((r) => r.product).join(', '));
      expect(fixed.warnings.some((w) => w.includes(exclusions[0]) && /matched nothing/i.test(w)), false,
        `TRF-07: the demo file's nitro line matched; warnings: ${fixed.warnings.join(' | ')}`);
      const old = await nitroRun(['[NX3](=O)=O', 'S(=O)(=O)[#6]']);
      expect(old.rows.some((r) => r.product === nitroProduct), true,
        `TRF-07: the old line excluded no nitro product; rows: ${old.rows.map((r) => r.product).join(', ')}`);
      expect(old.warnings.some((w) => w.includes('[NX3](=O)=O') && /matched nothing/i.test(w)), true,
        `TRF-07: the run names the old line as matching nothing; warnings: ${old.warnings.join(' | ')}`);
    }));
});

/** A stand-in for the cell editor's sketcher host: what it would write (`smiles`), its drawing (`molfile`), and a
 * change to tell the editor of (`draw`). Its implementation is no Crux sketcher: the empty-SMILES path is checked. */
function fakeSketcher(smiles: string, molfile: string): {host: DG.chem.Sketcher, draw: (s: string,
  mb: string) => void} {
  const onChanged = new Subject<any>();
  const state = {smiles, molfile};
  const host = {
    sketcher: null,
    onChanged,
    getSmiles: () => state.smiles,
    getMolFile: () => state.molfile,
    sketcherReady: () => new Promise(() => {}), // never: the editor updates on its own and on each change
  } as unknown as DG.chem.Sketcher;
  return {host, draw: (s: string, mb: string) => {
    state.smiles = s;
    state.molfile = mb;
    onChanged.next(null);
  }};
}
