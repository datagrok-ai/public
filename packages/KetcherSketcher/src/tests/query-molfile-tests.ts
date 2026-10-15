/* The query molfile Ketcher writes for Datagrok's substructure filter (src/query-molfile.ts; crux-sketch spike
 * query-roundtrip, K1): every SMARTS-only atom property and custom bond query Ketcher offers, alone and together, written
 * so that RDKit, as Chem's substructure search reads a query molblock, finds what the query's intended SMARTS finds; and
 * read back onto Ketcher as the query it was. */
import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {KetSerializer, MolSerializer} from 'ketcher-core';
import {asKetcherQuery, customQueryOf, mdlBondQuery, rdkitSmarts, smartsOnly, smartsqGroups, withQueries,
  withoutSmartsq} from '../query-molfile';
import * as DG from 'datagrok-api/dg';

/** Molecules the queries are matched over: carbons in chains and rings, aromatic and not, with 0 to 3 H, a charge. */
const PROBES = ['CO', 'CCO', 'OC1CCCCC1', 'Oc1ccccc1', 'CC(C)(C)O', 'C=CO', 'CC=O', 'CN', 'c1ccncc1', 'C1CC1O',
  'NCC(=O)O', 'CC(C)CO', 'C1CCCC1C', 'OCc1ccccc1', 'C[N+](C)(C)C', 'CC1CCCCC1', 'c1ccccc1C', 'C=C', 'OC(=O)C1CCCCC1'];

let rdkit: any = null;
function module(): any {
  if (rdkit === null) {
    const call = DG.Func.find({package: 'Chem', name: 'getRdKitModule'})[0].prepare();
    call.callSync();
    rdkit = call.getOutputParamValue();
  }
  return rdkit;
}

/** Which probes a query (SMARTS, or a molblock read as Chem's search reads a query molblock) matches: '1' and '0'. */
function matched(query: string): string {
  const rd = module();
  const q = rd.get_qmol(query);
  if (!q)
    throw new Error(`RDKit cannot read ${query}`);
  try {
    if (query.includes('M  END')) {
      try {
        q.convert_to_aromatic_form();
      } catch {
        // as Chem's search: a query that cannot be aromatized is matched as it is
      }
    }
    return PROBES.map((p) => {
      const m = rd.get_mol(p);
      try {
        return m.get_substruct_match(q) !== '{}' ? '1' : '0';
      } finally {
        m.delete();
      }
    }).join('');
  } finally {
    q.delete();
  }
}

type Ket = {atoms: any[], bonds: any[]};
const two = (a: any, b: any, bond: any = {}): Ket => ({atoms: [{location: [0, 0, 0], ...a}, {location: [1.5, 0, 0], ...b}],
  bonds: [{type: 1, atoms: [0, 1], ...bond}]});
const ketText = (k: Ket) => JSON.stringify({root: {nodes: [{$ref: 'mol0'}]}, mol0: {type: 'molecule', ...k}});

/** Ketcher's structure of a KET, its V2000 by ketcher-core, and the query molfile written for it. */
function written(k: Ket): {v2000: string, text: string} {
  const struct = new KetSerializer().deserialize(ketText(k));
  const v2000 = new MolSerializer().serialize(struct);
  return {v2000, text: withQueries(struct, v2000)};
}

// [what, the structure as Ketcher's tools leave it, its intended SMARTS]
const FEATURES: [string, Ket, string][] = [
  ['aromatic', two({label: 'C', queryProperties: {aromaticity: 'aromatic'}}, {label: 'O'}), '[#6;a]-[#8]'],
  ['aliphatic', two({label: 'C', queryProperties: {aromaticity: 'aliphatic'}}, {label: 'O'}), '[#6;A]-[#8]'],
  ['implicit H count 3', two({label: 'C', implicitHCount: 3}, {label: 'O'}), '[#6;h3]-[#8]'],
  ['implicit H count 0', two({label: 'C', implicitHCount: 0}, {label: 'O'}), '[#6;h0]-[#8]'],
  ['ring membership 1', two({label: 'C', queryProperties: {ringMembership: 1}}, {label: 'O'}), '[#6;R1]-[#8]'],
  ['ring membership 0', two({label: 'C', queryProperties: {ringMembership: 0}}, {label: 'O'}), '[#6;R0]-[#8]'],
  ['ring size 6', two({label: 'C', queryProperties: {ringSize: 6}}, {label: 'O'}), '[#6;r6]-[#8]'],
  ['ring size 3', two({label: 'C', queryProperties: {ringSize: 3}}, {label: 'O'}), '[#6;r3]-[#8]'],
  ['connectivity 4', two({label: 'C', queryProperties: {connectivity: 4}}, {label: 'O'}), '[#6;X4]-[#8]'],
  ['custom atom query', two({label: 'C', queryProperties: {customQuery: '#6;$([#6]=[#8])'}}, {label: 'O'}), '[#6;$([#6]=[#8])]-[#8]'],
  ['custom bond query, single or double in a ring', two({label: 'C'}, {label: 'C'}, {type: undefined, customQuery: '-,=;@'}), '[#6]-,=;@[#6]'],
  ['custom bond query, any in a chain', two({label: 'C'}, {label: 'O'}, {type: undefined, customQuery: '~;!@'}), '[#6]~;!@[#8]'],
  ['custom bond query, aromatic', two({label: 'C'}, {label: 'C'}, {type: undefined, customQuery: ':'}), '[#6]:[#6]'],
  // together
  ['a list with a ring size', two({type: 'atom-list', elements: ['C', 'N'], notList: false, queryProperties: {ringSize: 6}},
    {label: 'O'}), '[#6,#7;r6]-[#8]'],
  ['a NOT list, aliphatic', two({type: 'atom-list', elements: ['N', 'O'], notList: true, queryProperties: {aromaticity: 'aliphatic'}},
    {label: 'O'}), '[!#7&!#8;A]-[#8]'],
  ['ring bond count as drawn with connectivity 4', two({label: 'C', ringBondCount: -2, queryProperties: {connectivity: 4}},
    {label: 'O'}), '[#6;x0;X4]-[#8]'],
  ['ring bond count 2 with ring size 6', two({label: 'C', ringBondCount: 2, queryProperties: {ringSize: 6}}, {label: 'O'}),
    '[#6;x2;r6]-[#8]'],
  ['substitution count 3, aromatic', two({label: 'C', substitutionCount: 3, queryProperties: {aromaticity: 'aromatic'}},
    {label: 'O'}), '[#6;D3;a]-[#8]'],
  ['unsaturated, aliphatic', two({label: 'C', unsaturatedAtom: true, queryProperties: {aromaticity: 'aliphatic'}}, {label: 'O'}),
    '[#6;$(*=,:,#*);A]-[#8]'],
  ['H count 1 or more, ring membership 1', two({label: 'C', hCount: 2, queryProperties: {ringMembership: 1}}, {label: 'O'}),
    '[#6;h{1-};R1]-[#8]'],
  ['a charge, aliphatic', two({label: 'C'}, {label: 'N', charge: 1, queryProperties: {aromaticity: 'aliphatic'}}),
    '[#6]-[#7+;A]'],
  ['generic A with ring membership 1', two({label: 'A', queryProperties: {ringMembership: 1}}, {label: 'O'}), '[!#1;R1]-[#8]'],
  ['two atoms with properties and a custom bond', two({label: 'C', queryProperties: {ringSize: 6}},
    {label: 'C', queryProperties: {connectivity: 3}}, {type: undefined, customQuery: '-,:;@'}), '[#6;r6]-,:;@[#6;X3]'],
];

category('ketcher: query molfile', () => {
  test('the SMARTS-only atom properties as RDKit writes them, each and together', async () => {
    expect(smartsOnly({}).join(' '), '', 'no property');
    expect(smartsOnly({queryProperties: {aromaticity: 'aromatic'}}).join(' '), 'a');
    expect(smartsOnly({queryProperties: {aromaticity: 'aliphatic'}}).join(' '), 'A');
    expect(smartsOnly({implicitHCount: 0}).join(' '), 'h0');
    expect(smartsOnly({implicitHCount: 3}).join(' '), 'h3');
    expect(smartsOnly({queryProperties: {ringMembership: 0}}).join(' '), 'R0');
    expect(smartsOnly({queryProperties: {ringSize: 5}}).join(' '), 'r5');
    expect(smartsOnly({queryProperties: {connectivity: 2}}).join(' '), 'X2');
    expect(smartsOnly({implicitHCount: 1, queryProperties: {aromaticity: 'aromatic', ringMembership: 2, ringSize: 6,
      connectivity: 3}}).join(' '), 'a h1 R2 r6 X3');
    expect(smartsOnly({queryProperties: {aromaticity: null, ringMembership: null, ringSize: null, connectivity: null,
      customQuery: null}, implicitHCount: null}).join(' '), '', 'every field Any');
  });

  test('a custom bond query as MDL\'s type and topology, where MDL can say it', async () => {
    const cases: [string, string][] = [
      ['-', '1 0'], ['=', '2 0'], ['#', '3 0'], [':', '4 0'], ['-,=', '5 0'], ['=,-', '5 0'], ['-,:', '6 0'], ['=,:', '7 0'],
      ['~', '8 0'], ['@', '8 1'], ['!@', '8 2'], ['-,=;@', '5 1'], ['-,:;!@', '6 2'], [':;@', '4 1'], ['-&@', '1 1'],
      ['~;!@', '8 2'], [' -, = ; @ ', '5 1'],
      // what MDL cannot say, or a SMARTS parser reads otherwise
      ['-,=&@', 'none'], ['-,#', 'none'], ['-,=,:', 'none'], ['$([#6])', 'none'], ['/', 'none'], ['', 'none'], ['-;@;!@', 'none'],
    ];
    for (const [custom, want] of cases) {
      const got = mdlBondQuery(custom);
      expect(got === null ? 'none' : `${got.type} ${got.topology}`, want, `custom bond query "${custom}"`);
    }
  });

  test('each query feature Ketcher offers is written so that RDKit finds what its intended SMARTS finds', async () => {
    const problems: string[] = [];
    for (const [what, ket, intended] of FEATURES) {
      const want = matched(intended);
      if (!want.includes('1') || !want.includes('0'))
        problems.push(`${what}: the probes do not tell ${intended} apart (${want})`);
      const {text} = written(ket);
      const got = matched(text);
      if (got !== want)
        problems.push(`${what}: written ${text.includes('V3000') ? 'V3000' : 'V2000'}, RDKit finds ${got}, ${intended} finds ${want}`);
    }
    expect(problems.length, 0, problems.join('\n'));
  });

  test('a drawing with no SMARTS-only property and no custom bond query is written as ketcher-core writes it', async () => {
    for (const ket of [two({label: 'C'}, {label: 'O'}), two({label: 'C', ringBondCount: 2}, {label: 'O'}),
      two({type: 'atom-list', elements: ['N', 'O'], notList: false}, {label: 'C'}), two({label: 'C'}, {label: 'O'}, {type: 5, topology: 1})]) {
      const {v2000, text} = written(ket);
      expect(text, v2000);
    }
  });

  test('a query molfile written here shows Ketcher its structure again; another V3000 with SMARTSQ, its SMARTS as custom queries', async () => {
    const ket = two({label: 'C', queryProperties: {ringSize: 6}}, {label: 'O'});
    const struct = new KetSerializer().deserialize(ketText(ket));
    const text = withQueries(struct, new MolSerializer().serialize(struct));
    expect(text.includes('V3000') && smartsqGroups(text).size === 1, true, 'written V3000 with one SMARTSQ group');
    const back = JSON.parse(asKetcherQuery(text));
    expect(back.mol0.atoms[0].queryProperties?.ringSize, 6, 'the ring size, as drawn');
    // a V3000 another sketcher wrote (Crux's: RDKit's SMARTSQ groups)
    const crux = ['', '     RDKit          2D', '', '  0  0  0  0  0  0  0  0  0  0999 V3000', 'M  V30 BEGIN CTAB',
      'M  V30 COUNTS 2 1 1 0 0', 'M  V30 BEGIN ATOM', 'M  V30 1 C 0.000000 0.000000 0.000000 0',
      'M  V30 2 O 1.299038 0.750000 0.000000 0', 'M  V30 END ATOM', 'M  V30 BEGIN BOND', 'M  V30 1 1 1 2', 'M  V30 END BOND',
      'M  V30 BEGIN SGROUP', 'M  V30 1 DAT 0 ATOMS=(1 1) QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&r6]"',
      'M  V30 END SGROUP', 'M  V30 END CTAB', 'M  END', ''].join('\n');
    expect(customQueryOf('[#6&r6:3]'), '#6&r6');
    expect(withoutSmartsq(crux)!.includes('SGROUP'), false, 'the SMARTSQ groups gone');
    const shown = JSON.parse(asKetcherQuery(crux));
    expect(shown.mol0.atoms[0].queryProperties?.customQuery, '#6&r6', 'Crux\'s ring size as Ketcher\'s custom query');
    expect(matched(rdkitSmarts(crux)), matched('[#6;r6]-[#8]'), 'the SMARTS read as RDKit reads it');
    // anything else is shown as it is
    expect(asKetcherQuery('CCO'), 'CCO');
  });
});
