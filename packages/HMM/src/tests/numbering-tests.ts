import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';
import {HmmerPool} from '../pool';
import {allowedChains, numberingFrame, SCHEMES} from '../numbering';
import ANARCI_FIXTURE from './anarci-fixture.json';
import OPTIONS_FIXTURE from './anarci-options-fixture.json';
import type {Scheme} from '../hmmer/anarci/types.ts';

const BIO_COLUMNS = ['position_names', 'chain_type', 'annotations_json', 'numbering_detail', 'numbering_map'];

function fixtureColumn(): DG.Column<string> {
  const col = DG.Column.fromStrings('sequence', ANARCI_FIXTURE.map((f) => f.sequence));
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = 'fasta';
  col.setTag('aligned', 'SEQ');
  col.setTag('alphabet', 'PT');
  return col as DG.Column<string>;
}

category('ANARCI numbering', () => {
  for (const scheme of SCHEMES) {
    test(`identical to ANARCI: ${scheme}`, async () => {
      const df = DG.DataFrame.fromColumns([fixtureColumn()]);
      const result: DG.DataFrame = await grok.functions.call('HMM:anarciNumbering',
        {df, seqCol: df.col('sequence'), scheme});
      for (const name of BIO_COLUMNS) expect(result.col(name) !== null, true, `missing column ${name}`);
      ANARCI_FIXTURE.forEach((f, i) => {
        const expected = f.expected[scheme as Scheme];
        const positions = result.get('position_names', i) ?? '';
        if (expected === null) {
          expect(positions, '', `${f.name}: ANARCI numbers nothing`);
          return;
        }
        expect(positions, expected.positions, `${f.name} positions`);
        const detail: {position: string; aa: string}[] = JSON.parse(result.get('numbering_detail', i) || '[]');
        expect(detail.map((d) => d.aa).join(''), expected.residues, `${f.name} residues`);
        expect(result.get('species', i), expected.species, `${f.name} species`);
        expect(result.get('evalue', i), expected.evalue, `${f.name} E-value`);
        expect(result.get('bitscore', i), expected.bitscore, `${f.name} bit score`);
        expect(result.get('v_gene', i) || null, expected.vGene, `${f.name} V gene`);
        expect(result.get('j_gene', i) || null, expected.jGene, `${f.name} J gene`);
        expect(result.get('domains', i), expected.domains, `${f.name} domains`);
        // numbering_map points into the gap-free sequence at the numbered residues.
        const map: Record<string, number> = JSON.parse(result.get('numbering_map', i) || '{}');
        for (const d of detail) expect(f.sequence[map[d.position]], d.aa, `${f.name} ${d.position}`);
      });
    });
  }

  test('options (species, bit score threshold) identical to ANARCI', async () => {
    for (const setting of [...new Set(OPTIONS_FIXTURE.map((f) => f.setting))]) {
      const cases = OPTIONS_FIXTURE.filter((f) => f.setting === setting);
      const [label, threshold] = setting.split('@');
      const col = DG.Column.fromStrings('sequence', cases.map((f) => f.sequence));
      col.semType = DG.SEMTYPE.MACROMOLECULE;
      col.meta.units = 'fasta';
      const df = DG.DataFrame.fromColumns([col]);
      const result: DG.DataFrame = await grok.functions.call('HMM:anarciNumbering', {df, seqCol: col, scheme: 'imgt',
        species: label === 'human_mouse' ? 'human and mouse' : label, bitScoreThreshold: Number(threshold)});
      cases.forEach((f, i) => {
        expect(result.get('position_names', i) ?? '', f.expected?.positions ?? '', `${setting} ${f.name} positions`);
        expect(result.get('species', i) || null, f.expected?.species ?? null, `${setting} ${f.name} species`);
        expect(result.get('v_gene', i) || null, f.expected?.vGene ?? null, `${setting} ${f.name} V gene`);
        expect(result.get('j_gene', i) || null, f.expected?.jGene ?? null, `${setting} ${f.name} J gene`);
      });
    }
  });

  test('separator notation numbers like FASTA', async () => {
    const fasta = fixtureColumn();
    const separator = DG.Column.fromStrings('sequence', ANARCI_FIXTURE.map((f) => f.sequence.split('').join('-')));
    separator.semType = DG.SEMTYPE.MACROMOLECULE;
    separator.meta.units = 'separator';
    separator.setTag('separator', '-');
    separator.setTag('alphabet', 'PT');
    const run = async (col: DG.Column) => (await grok.functions.call('HMM:anarciNumbering',
      {df: DG.DataFrame.fromColumns([col]), seqCol: col, scheme: 'imgt'}) as DG.DataFrame);
    const [a, b] = [await run(fasta), await run(separator)];
    for (const name of BIO_COLUMNS)
      expect(JSON.stringify(b.col(name)!.toList()), JSON.stringify(a.col(name)!.toList()), name);
  });

  test('worker count does not change results', async () => {
    const sequences = ANARCI_FIXTURE.map((f, i): [string, string] => [String(i), f.sequence]);
    const many = [...sequences, ...sequences, ...sequences, ...sequences]
      .map(([, s], i): [string, string] => [String(i), s]);
    const options = {scheme: 'kabat' as Scheme, allow: allowedChains('kabat'), assignGermline: true};
    const one = await HmmerPool.get(_package).anarci(many, options, 1);
    const all = await HmmerPool.get(_package).anarci(many, options, 8);
    expect(JSON.stringify(all), JSON.stringify(one));
  });

  test('region annotations resolve in every row', async () => {
    const sequences = ANARCI_FIXTURE.map((f, i): [string, string] => [String(i), f.sequence]);
    const results = await HmmerPool.get(_package).anarci(sequences,
      {scheme: 'wolfguy', allow: allowedChains('wolfguy'), assignGermline: false});
    const frame = numberingFrame('wolfguy', results);
    for (let i = 0; i < frame.rowCount; i++) {
      if (!frame.get('position_names', i)) continue;
      const map = JSON.parse(frame.get('numbering_map', i));
      const regions: {start: string; end: string}[] = JSON.parse(frame.get('annotations_json', i));
      expect(regions.length, 7, `row ${i} regions`);
      for (const r of regions) expect(map[r.start] !== undefined && map[r.end] !== undefined, true);
    }
    // TCR chains get IMGT regions from Bio's tables too.
    const imgt = numberingFrame('imgt', await HmmerPool.get(_package).anarci(sequences,
      {scheme: 'imgt', allow: allowedChains('imgt'), assignGermline: false}));
    let tcr = 0;
    for (let i = 0; i < imgt.rowCount; i++) {
      if (!['Alpha', 'Beta', 'Gamma', 'Delta'].includes(imgt.get('chain', i))) continue;
      tcr++;
      const regions: {name: string}[] = JSON.parse(imgt.get('annotations_json', i));
      expect(regions.map((r) => r.name).join(','), 'FR1,CDR1,FR2,CDR2,FR3,CDR3,FR4', `TCR row ${i} regions`);
    }
    expect(tcr > 0, true, 'fixture has TCR chains');
  });
});
