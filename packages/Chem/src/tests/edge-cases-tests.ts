/* eslint-disable max-len */
/**
 * Edge inputs (night 2026-11-02): blank and broken cells, odd chemistry, queries no one reads, analyses with nothing to
 * find, in both chemical engines. Each test is a table of inputs; RDKit mode's answer is the reference unless a
 * comment says why Crux's differs (RDKit 2026.03 vs Chem's RDKit.js 2024.09).
 */
import * as DG from 'datagrok-api/dg';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemicalEngine} from '../engine/chem-engine';
import {convertMolNotationByEngine, fingerprintsByEngine, inchisByEngine} from '../engine/worker-jobs';
import {mcsByEngine, rGroupsByEngine} from '../engine/analysis-jobs';
import {chemSubstructureSearchLibrary} from '../chem-searches';
import {Fingerprint} from '../utils/chem-common';
import {RDKitCellRenderer} from '../rendering/rdkit-cell-renderer';
import {PackageFunctions} from '../package';
import {inBothEngines, inEngine} from './engine-utils';
import {getCruxService} from '../crux/crux-searches';
import {cruxModuleSync} from '../engine/crux-module';
import {_package} from '../package-test';

const molblock = (lines: string[]) => lines.join('\n');
const BENZENE_V2000 = molblock(['', '  edge', '', '  6  6  0  0  0  0  0  0  0  0999 V2000',
  '    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000   -1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '   -1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '    0.0000    1.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
  '  1  2  2  0', '  2  3  1  0', '  3  4  2  0', '  4  5  1  0', '  5  6  2  0', '  6  1  1  0', 'M  END']);

/** The edge column: name → cell (what crux-core's spike n1102-chem-edge changes is tested apart, below). */
const EDGE: [string, string][] = [
  ['empty', ''],
  ['space', ' '],
  ['tab and newline', '\t\n'],
  ['methane', 'C'],
  ['benzene, Kekulé', 'C1=CC=CC=C1'],
  ['salt', 'CC(=O)Oc1ccccc1C(=O)O.[Na+]'],
  ['isotopes', '[2H]C([2H])([2H])[13CH3]'],
  ['radical', '[CH3]'],
  ['charges', 'C[N+](C)(C)CC(=O)[O-]'],
  ['dative', 'C[NH2]->[Pt](Cl)(Cl)<-[NH2]C'],
  ['ferrocene', '[cH-]1cccc1.[cH-]1cccc1.[Fe+2]'],
  ['azulene', 'c1ccc2cccc2cc1'],
  ['pyrrole N-H', 'c1cc[nH]c1'],
  ['tropylium', 'c1cc[cH+]ccc1'],
  ['hypervalent C', 'C(C)(C)(C)(C)C'],
  ['stereo groups', 'C[C@H](F)[C@@H](C)O |&1:1,o2:3|'],
  ['R groups', '[*:1]c1ccccc1[*:2]'],
  ['SMARTS', '[#6]~[#7]'],
  ['invalid paren', 'CC(C'],
  ['invalid element', 'C[Xx]'],
  ['quaternary N, uncharged', 'CN(C)(C)CC(=O)O'],
  ['name', 'aspirin'],
  ['NaN', 'NaN'],
  ['molblock', BENZENE_V2000],
  ['molblock CRLF', BENZENE_V2000.replace(/\n/g, '\r\n')],
];

const molecules = () => EDGE.map(([, m]) => m);
/** Per row whose two values differ: `name: rdkit | crux`. */
function differences<T>(rdkit: T[], crux: T[], show: (v: T) => string = (v) => JSON.stringify(v)): string[] {
  return EDGE.flatMap(([name], i) => show(rdkit[i]) === show(crux[i]) ? [] :
    [`${name}: RDKit ${show(rdkit[i])} | Crux ${show(crux[i])}`]);
}
/** A conversion's result as compared: a MOL block's counts line and atom symbols (no coordinates: the layouts of RDKit
 * 2024.09 and 2026.03 differ), any other string as it is. */
const molblockShape = (s: string | null) => !s?.includes('M  END') ? s :
  s.split('\n').slice(3).filter((l) => /^ *\d+ *\d+ .*V2000$|^ +-?\d+\.\d+ +-?\d+\.\d+ +-?\d+\.\d+ /.test(l))
    .map((l) => l.includes('V2000') ? l.trim() : l.slice(31, 34).trim()).join(' ');
const short = (s: string | null | undefined) => s === null || s === undefined ? String(s) :
  JSON.stringify(s.length > 60 ? s.slice(0, 60) + '…' : s);

/**
 * What the census finds RDKit mode and Crux mode doing differently, each for a reason (stand runs of 2026-11-03):
 * - Crux writes RDKit 2026.03's strings where Chem's RDKit.js 2024.09 wrote others: the dative SMILES' bracket atoms,
 *   and SMARTS of molecules without their atoms' H counts (`[#7]` for `[nH]`; RDKit 2026.03.6's MolToSmarts, Python's
 *   too);
 * - the search's RDKit side for the rows Crux does not index (crux-searches.ts `rdkitMatcher`) matches a cell RDKit
 *   reads only as a query (an uncharged quaternary N, a pentavalent C, a SMARTS such as `NaN` or `C[Xx]`) as the query
 *   it is; RDKit mode prescreens by that query's pattern fingerprint, which keeps some of them and drops others with no
 *   rule a fallback could follow (it finds `NaN` for `N` but not for `*`, the quaternary N for `N` but not for `CN`).
 *   Crux's answer is the one kept: right for the molecules, arbitrary like RDKit mode's for the non-molecules.
 */
const CENSUS_DIFFERENCES = [
  'convert smiles: dative: RDKit "CN->[Pt](Cl)(Cl)<-NC" | Crux "C[NH2]->[Pt]([Cl])([Cl])<-[NH2]C"',
  'convert smarts: isotopes: RDKit "[2H]-[#6](-[2H])(-[2H])-[13#6H3]" | Crux "[2H]-[#6](-[2H])(-[2H])-[13#6]"',
  'convert smarts: radical: RDKit "[#6H3]" | Crux "[#6]"',
  'convert smarts: dative: RDKit "[#6]-[#7H2]->[Pt](-[#17])(-[#17])<-[#7H2]-[#6]" | Crux "[#6]-[#7]->[Pt](-[#17])(-[#17])<-[#7]-[#6]"',
  'convert smarts: ferrocene: RDKit "[#6H-]1:[#6]:[#6]:[#6]:[#6]:1.[#6H-]1:[#6]:[#6]:[#6]:[#6]:1.…" | Crux "[#6-]1:[#6]:[#6]:[#6]:[#6]:1.[#6-]1:[#6]:[#6]:[#6]:[#6]:1.[F…"',
  'convert smarts: pyrrole N-H: RDKit "[#6]1:[#6]:[#6]:[#7H]:[#6]:1" | Crux "[#6]1:[#6]:[#6]:[#7]:[#6]:1"',
  'convert smarts: tropylium: RDKit "[#6]1:[#6]:[#6]:[#6H+]:[#6]:[#6]:[#6]:1" | Crux "[#6]1:[#6]:[#6]:[#6+]:[#6]:[#6]:[#6]:1"',
  'contains *: RDKit rows R groups | Crux rows R groups, invalid element, NaN',
  'contains CC: RDKit rows salt, isotopes, charges, azulene, stereo groups | Crux rows salt, isotopes, charges, azulene, ' +
    'hypervalent C, stereo groups, quaternary N, uncharged',
  'contains C(=O)O: RDKit rows salt, charges | Crux rows salt, charges, quaternary N, uncharged',
];

category('chemical engine: edge cases', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  test('census: every edge cell through each function, both engines', async () => {
    const found: string[] = [];
    const col = DG.Column.fromStrings('molecules', molecules());
    col.semType = DG.SEMTYPE.MOLECULE;
    const runs = await inBothEngines(async () => {
      const out: {[what: string]: any[]} = {};
      for (const notation of [DG.chem.Notation.Smiles, DG.chem.Notation.MolBlock, DG.chem.Notation.Smarts])
        out[`convert ${notation}`] = (await convertMolNotationByEngine(molecules(), notation, false)).map(molblockShape);
      out['InChI'] = await inchisByEngine(molecules(), false);
      out['Morgan null'] = (await fingerprintsByEngine(Fingerprint.Morgan, molecules(), false)).fps.map((f) => f === null);
      out['Pattern null'] = (await fingerprintsByEngine(Fingerprint.Pattern, molecules(), false)).fps.map((f) => f === null);
      const renderer = new RDKitCellRenderer(chemCommonRdKit.getRdKitModule());
      out['drawn'] = molecules().map((m) => m === '' ? null : renderer._fetchMol(m, [], false, false, {}, false).molCtx.mol !== null);
      out['isSmiles'] = molecules().map((m) => PackageFunctions.isSmiles(m));
      out['isSmarts'] = molecules().map((m) => PackageFunctions.isSmarts(m));
      for (const query of ['C', 'c1ccccc1', '*', '[#7]', 'CC', 'C(=O)O']) {
        const hits = await chemSubstructureSearchLibrary(col, query, '');
        out[`contains ${query}`] = EDGE.flatMap((_, i) => hits.getBit(i) ? [i] : []);
      }
      return out;
    });
    for (const what of Object.keys(runs.rdkit)) {
      const [r, c] = [runs.rdkit[what], runs.crux[what]];
      if (what.startsWith('contains')) {
        if (JSON.stringify(r) !== JSON.stringify(c))
          found.push(`${what}: RDKit rows ${r.map((i: number) => EDGE[i][0]).join(', ')} | Crux rows ${c.map((i: number) => EDGE[i][0]).join(', ')}`);
      } else
        found.push(...differences(r, c, (v: any) => typeof v === 'string' ? short(v) : JSON.stringify(v)).map((d) => `${what}: ${d}`));
    }
    // the analyses with nothing to find, Crux mode alone (RDKit mode's MCS and R-groups throw or return nothing alike)
    const analyses = await inBothEngines(async () => ({
      mcsNone: await mcsByEngine(['C', 'O', 'N'], true, true).catch((e) => `throws ${e?.message ?? e}`),
      mcsBlank: await mcsByEngine(['c1ccccc1C', ' ', 'c1ccccc1O'], true, true).catch((e) => `throws ${e?.message ?? e}`),
      rgroupsNone: await rGroupsByEngine(['CCO', 'CCN'], 'c1ccccc1', false).then((r) => JSON.stringify(r).slice(0, 200))
        .catch((e) => `throws ${e?.message ?? e}`),
    }));
    for (const k of Object.keys(analyses.rdkit)) {
      const [r, c] = [(analyses.rdkit as any)[k], (analyses.crux as any)[k]];
      if (JSON.stringify(r) !== JSON.stringify(c))
        found.push(`${k}: RDKit ${short(String(r))} | Crux ${short(String(c))}`);
    }
    // (grok test prints a failure's message after its last ': ', so none is in it)
    const shown = (list: string[]) => list.join(' ## ').replace(/: /g, ' = ');
    expect(shown(found), shown(CENSUS_DIFFERENCES));
  }, {timeout: 300000});

  test('Crux reads ring closures, blanks and truncated MOL blocks as RDKit does', async () => {
    // RDKit 2026.03.6 (and Chem's RDKit.js) reads none of these as a molecule or a query; Crux read `C11` as methane,
    // `C1C1` / `C12CC12` as a molecule with two bonds between the same atoms (written `C1`, `C1C2`), and a blank string
    // or a MOL block without `M  END` as an empty query (crux-core spike n1102-chem-edge: smi-compiler's ring closures,
    // get_qmol's blanks); get_qmol drops the blanks at the ends (`  c1ccccc1` is benzene, it was an empty query)
    await inEngine(ChemicalEngine.Crux, async () => {
      const crux = cruxModuleSync()!;
      const reads = (s: string) => {
        const mol = crux.get_mol(s, '{}');
        const query = mol ? null : crux.get_qmol(s);
        try {
          return mol ? `mol ${mol.get_smiles()}` : query ? `query ${query.get_smarts()}` : 'none';
        } finally {
          mol?.delete();
          query?.delete();
        }
      };
      const cells: [string, string][] = [['C11', 'none'], ['C1C1', 'none'], ['C12CC12', 'none'], ['C%101CC%101', 'none'],
        [' ', 'none'], ['\t\n', 'none'], [BENZENE_V2000.split('\n').slice(0, 8).join('\n'), 'none'],
        [BENZENE_V2000.replace('M  END', ''), 'none'], ['C1CC1', 'mol C1CC1'], ['C12CC1C2', 'mol C1C2CC12']];
      expect(cells.map(([s]) => `${JSON.stringify(s.slice(0, 12))}: ${reads(s)}`).join('; '),
        cells.map(([s, r]) => `${JSON.stringify(s.slice(0, 12))}: ${r}`).join('; '));
      const query = crux.get_qmol('  c1ccccc1');
      expect(query?.get_num_atoms() ?? 0, 6, 'a query with blanks at its start');
      query?.delete();
      // get_aromatic_form writes the aromatic bonds (MinimalLib's kekulize = false): azulene's ten ring bonds
      const azulene = crux.get_mol('c1ccc2cccc2cc1');
      const block = azulene!.get_aromatic_form();
      azulene!.delete();
      expect(block.split('\n').slice(14, 25).map((l) => l.trim().split(/\s+/)[2]).join(''), '44444444441');
    });
  });

  test('a blank cell is no molecule in either engine: crossed out, not converted', async () => {
    const blanks = [' ', '\t\n', '\r\n', '   '];
    const runs = await inBothEngines(async () => {
      const renderer = new RDKitCellRenderer(chemCommonRdKit.getRdKitModule());
      return {
        drawn: blanks.map((m) => renderer._fetchMol(m, [], false, false, {}, false).molCtx.mol !== null).join(','),
        smiles: (await convertMolNotationByEngine(blanks, DG.chem.Notation.Smiles, false)).join(','),
        isSmiles: blanks.map((m) => PackageFunctions.isSmiles(m)).join(','),
      };
    });
    expect(JSON.stringify(runs.crux), JSON.stringify(runs.rdkit));
    expect(runs.crux.drawn, 'false,false,false,false');
  });

  test('a trap the RDKit worker code swallows fails its job, and the workers are replaced', async () => {
    // a ring of 20,000 atoms: writing its SMILES overflows the stack (V8's RangeError inside the wasm call), which the
    // RDKit worker code's per-molecule `catch` swallowed, leaving a worker whose stack pointer the overflow moved
    const n = 20000;
    const ring = ['', '  ring', '', '  0  0  0     0  0            999 V3000', 'M  V30 BEGIN CTAB', `M  V30 COUNTS ${n} ${n} 0 0 0`,
      'M  V30 BEGIN ATOM', ...Array.from({length: n}, (_, i) => `M  V30 ${i + 1} C ${(10 * Math.cos(2 * Math.PI * i / n)).toFixed(4)} ${(10 * Math.sin(2 * Math.PI * i / n)).toFixed(4)} 0 0`),
      'M  V30 END ATOM', 'M  V30 BEGIN BOND', ...Array.from({length: n}, (_, i) => `M  V30 ${i + 1} 1 ${i + 1} ${(i + 1) % n + 1}`),
      'M  V30 END BOND', 'M  V30 END CTAB', 'M  END'].join('\n');
    await inEngine(ChemicalEngine.Crux, async () => {
      const service = await getCruxService();
      const generation = service.generation;
      const outcome = await service.convertMolNotation([ring, 'CCO'], DG.chem.Notation.Smiles, false)
        .then((res) => ({error: null, ring: res.results[0] ?? ''}), (e) => ({error: String(e?.message ?? e), ring: ''}));
      if (outcome.error === null && outcome.ring.length > n) {
        // Crux's SMILES writer no longer recurses (crux-core spike n1102-b3 vendored): the ring is written, nothing
        // overflows, so this input has no trap to watch for (a swallowed one leaves the ring unconverted, below)
        expect(outcome.ring.startsWith('C1'), true, 'the ring\'s SMILES');
        return;
      }
      const error = outcome.error;
      expect(error !== null && /RuntimeError/.test(error), true, `the job failed as a trap: ${error}`);
      expect(service.generation > generation, true, 'the workers were replaced');
      const after = await service.convertMolNotation(['C1CC1', 'CCO'], DG.chem.Notation.Smiles, false);
      expect(after.results.join(' '), 'C1CC1 CCO', 'the new workers convert');
    });
  }, {timeout: 120000});

  test('known differences, pinned: dative bonds, cells RDKit reads only as a query', async () => {
    // a dative bond is a single bond to Crux's descriptor graph (crux-core engine gaps EG-15, EG-30): the N->Pt bonds of
    // this cisplatin count as rotatable, and the NH2 donors / acceptors are not seen (RDKit: 0 rotatable, 2 and 2)
    const dative = 'C[NH2]->[Pt](Cl)(Cl)<-[NH2]C';
    const counts = (module: any) => {
      const mol = module.get_mol(dative);
      try {
        const d = JSON.parse(mol.get_descriptors());
        return `rotatable ${d.NumRotatableBonds}, donors ${d.NumHBD}, acceptors ${d.NumHBA}`;
      } finally {
        mol.delete();
      }
    };
    await inEngine(ChemicalEngine.Crux, async () =>
      expect(counts(cruxModuleSync()), 'rotatable 2, donors 0, acceptors 0', 'Crux (a known gap)'));
    expect(counts(chemCommonRdKit.getRdKitModule()), 'rotatable 0, donors 2, acceptors 2', 'RDKit');

    // a cell RDKit reads only as a query (an uncharged quaternary N, a pentavalent C): Crux mode matches the query it is,
    // RDKit mode's prescreen by that query's pattern fingerprint misses it for most queries (it finds the quaternary N
    // for N, but not for CN): Crux's answer is the one kept
    const cells = ['CN(C)(C)CC(=O)O', 'C(C)(C)(C)(C)C', 'CC(=O)O', 'CCN'];
    const col = DG.Column.fromStrings('molecules', cells);
    col.semType = DG.SEMTYPE.MOLECULE;
    const hits = await inBothEngines(async () => {
      const out: string[] = [];
      for (const query of ['N', 'CN', 'CC', 'C(=O)O']) {
        const found = await chemSubstructureSearchLibrary(col, query, '');
        out.push(`${query}: ${cells.filter((_, i) => found.getBit(i)).join(' ')}`);
      }
      return out.join(' ## ');
    });
    expect(hits.crux, 'N: CN(C)(C)CC(=O)O CCN ## CN: CN(C)(C)CC(=O)O CCN ## CC: CN(C)(C)CC(=O)O C(C)(C)(C)(C)C CC(=O)O CCN ## ' +
      'C(=O)O: CN(C)(C)CC(=O)O CC(=O)O', 'Crux mode');
    expect(hits.rdkit, 'N: CN(C)(C)CC(=O)O CCN ## CN: CCN ## CC: CC(=O)O CCN ## C(=O)O: CC(=O)O', 'RDKit mode');
  });
});
