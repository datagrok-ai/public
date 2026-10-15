import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package-test';

/** Fixtures (synced from Rusty-HMMER web/align): 60 real antibody chains and
 * C hmmalign's aligned rows to Pfam V-set, plain and --trim. */
async function alignFixture(): Promise<{col: DG.Column<string>; expected: Record<string, string>[]}> {
  const df = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/align-sequences.csv'));
  const col = df.col('sequence') as DG.Column<string>;
  col.semType = DG.SEMTYPE.MACROMOLECULE;
  col.meta.units = 'fasta';
  col.setTag('alphabet', 'PT');
  const [header, ...lines] = (await _package.files.readAsText('tests/align-expected.csv')).trim().split('\n');
  const names = header.split(',');
  return {col, expected: lines.map((l) => Object.fromEntries(l.split(',').map((v, i) => [names[i], v])))};
}

/** hmmalign's afa row as the Bio column holds it (upper case, `-` gaps). */
const asColumn = (row: string) => row.toUpperCase().replace(/\./g, '-');

category('HMMER alignment', () => {
  test('align to profile matches C hmmalign (V-set, plain and trimmed)', async () => {
    const {col, expected} = await alignFixture();
    for (const trim of [false, true]) {
      const result: DG.Column<string> = await grok.functions.call('HMM:alignToProfile',
        {sequences: col, profile: 'Pfam family or accession', family: 'V-set', trim});
      expect(result.length, expected.length, 'rows');
      expected.forEach((e, i) =>
        expect(result.get(i), asColumn(trim ? e.trimmed : e.aligned), `row ${i} trim=${trim}`));
      expect(result.getTag('aligned'), 'SEQ.MSA');
      const names = result.getTag('.positionNames').split(', ');
      expect(names.length, result.get(0)!.length, 'a position name per column');
      expect(names.filter((n) => !n.includes('.')).length, 109, 'V-set has 109 match positions');
    }
  });

  test('align to profile picks a germline model for antibodies', async () => {
    const {col} = await alignFixture();
    const result: DG.Column<string> = await grok.functions.call('HMM:alignToProfile',
      {sequences: col, profile: 'Antibody/TCR germline (auto)'});
    const lengths = new Set(result.toList().map((r: string) => r.length));
    expect(lengths.size, 1, 'equal row lengths');
    expect(/^(human|mouse|rat|rabbit|rhesus|pig|alpaca|cow)_[HKLABGD]$/.test(result.getTag('.hmmProfile')), true,
      `profile ${result.getTag('.hmmProfile')}`);
    for (let i = 0; i < col.length; i++)
      expect(result.get(i)!.replace(/-/g, ''), col.get(i)!.toUpperCase(), `row ${i} keeps its residues`);
  });
});

/** The WebLogo profile a package function supplies (Bio's `webLogoProfile` role). */
interface Profile {
  title: string;
  positions: {name: string; heights: Record<string, number>; info: Record<string, string>}[];
  columnPositions?: number[];
}

async function webLogoProfile(table: DG.DataFrame, column: string, profile: string): Promise<Profile> {
  return grok.functions.call('HMM:hmmWebLogoProfile', {table, column, profile});
}

category('HMMER logo', () => {
  test('logo profile values match C hmmlogo (V-set)', async () => {
    const expected = (await _package.files.readAsText('tests/logo-vset.txt')).trim().split('\n');
    const heights = expected.slice(2, expected.indexOf('Indel values'));
    const indels = expected.slice(expected.indexOf('Indel values') + 1);
    const profile = await webLogoProfile(DG.DataFrame.create(1), '', 'V-set');
    expect(profile.positions.length, heights.length, 'positions');
    // %6.3f of exact float values; ties at three decimals do not occur in this model.
    heights.forEach((line, i) => {
      const printed = line.split(':')[1].split('(')[0].trim().split(/\s+/);
      const h = profile.positions[i].heights;
      expect(Array.from('ACDEFGHIKLMNPQRSTVWY', (a) => h[a].toFixed(3)).join(' '), printed.join(' '),
        `position ${i + 1}`);
      const [insertP, insertLength, occupancy] = indels[i].split(':')[1].trim().split(/\s+/);
      const info = profile.positions[i].info;
      expect([info['Insert probability'], info['Expected insert length'], info['Occupancy']].join(' '),
        [insertP, insertLength, occupancy].join(' '), `indels ${i + 1}`);
    });
  });

  test('WebLogo draws the profile; letters select rows of the aligned column', async () => {
    const {col} = await alignFixture();
    const aligned: DG.Column<string> = await grok.functions.call('HMM:alignToProfile',
      {sequences: col, profile: 'Pfam family or accession', family: 'V-set', trim: false});
    const df = DG.DataFrame.fromColumns([aligned]);
    // Positions of an hmmalign column of the profile are its match positions.
    const profile = await webLogoProfile(df, aligned.name, 'V-set');
    const names = aligned.getTag('.positionNames').split(', ');
    expect(profile.columnPositions?.length, 109, 'a column position per match position');
    expect(profile.columnPositions!.every((at, i) => names[at] === String(i + 1)), true, 'match columns');
    expect((await webLogoProfile(df, aligned.name, 'C1-set')).columnPositions === undefined, true, 'other profiles');

    const tv = grok.shell.addTableView(df);
    try {
      const viewer = await df.plot.fromType('WebLogo',
        {profile: 'V-set', sequenceColumnName: aligned.name}) as DG.Viewer;
      tv.dockManager.dock(viewer, DG.DOCK_TYPE.DOWN);
      for (let i = 0; i < 150 && ((viewer as any).positions?.length ?? 0) !== 109; i++)
        await new Promise((r) => setTimeout(r, 100));
      expect((viewer as any).positions.length, 109, 'WebLogo positions');
    } finally {
      tv.close();
    }
  });
});

category('HMMER profile building', () => {
  test('build profile equals C hmmbuild, search with it equals C hmmsearch', async () => {
    const alignment = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/build-alignment.csv'));
    const aligned = alignment.col('aligned') as DG.Column<string>;
    aligned.semType = DG.SEMTYPE.MACROMOLECULE;
    aligned.meta.units = 'fasta';
    aligned.setTag('aligned', 'SEQ.MSA');
    aligned.setTag('alphabet', 'PT');
    const text: string = await grok.functions.call('HMM:buildProfile',
      {table: alignment, alignment: aligned, name: 'c1set-small'});
    const undated = (t: string) => t.replace(/^DATE .*$/m, 'DATE');
    expect(undated(text), undated(await _package.files.readAsText('tests/build-model.hmm')), 'model text');

    const {mapColumns, storeBuiltProfile} = await import('../logo');
    storeBuiltProfile(alignment, 'c1set-small', text, aligned);
    const profile = await webLogoProfile(alignment, aligned.name, 'c1set-small');
    const map = mapColumns(text);
    expect(map.length, Number(/^LENG\s+(\d+)/m.exec(text)![1]), 'a MAP column per match position');
    expect(map.every((c, i) => c > (map[i - 1] ?? 0) && c <= aligned.get(0)!.length), true, 'increasing columns');
    expect(profile.columnPositions?.join(' '), map.map((c) => c - 1).join(' '), 'positions of the alignment');

    const targets = DG.DataFrame.fromCsv(await _package.files.readAsText('tests/phmmer-targets.csv'));
    const col = targets.col('sequence') as DG.Column<string>;
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    col.meta.units = 'fasta';
    col.setTag('alphabet', 'PT');
    const result: DG.DataFrame = await grok.functions.call('HMM:searchWithHmm',
      {table: targets, sequence: col, model: text, evalue: 10, gathering: false});
    const [header, ...lines] = (await _package.files.readAsText('tests/build-expected.csv')).trim().split('\n');
    const names = header.split(',');
    const expected = lines.map((l) => Object.fromEntries(l.split(',').map((v, i) => [names[i], v])));
    const rowNames = targets.col('name')!.toList() as string[];
    for (const e of expected) {
      const row = rowNames.indexOf(e.target);
      expect(Number(result.columns.byIndex(0).get(row).toFixed(1)), Number(e.score), `${e.target} score`);
      expect(Number(result.columns.byIndex(1).get(row).toPrecision(2)), Number(e.evalue), `${e.target} E-value`);
      expect(result.columns.byIndex(2).get(row), Number(e.reported_domains), `${e.target} domains`);
    }
    expect(result.columns.byIndex(0).stats.valueCount, expected.length, 'reported sequences');
  });
});
