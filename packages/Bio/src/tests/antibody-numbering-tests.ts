import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {_package} from '../package';
import {numberSequencesWithImmunum} from '../utils/antibody-numbering/immunum-client';
import {numberAntibodyColumn} from '../utils/antibody-numbering/number-antibody';

/** Canonical test sequences picked from samples/antibodies.csv.
 *  - heavyChain1/2 are IGH variable regions starting with the classic EVQL/QVQL motifs
 *  - lightChain1/2 are IGK/IGL variable regions (DIQM/DIVM/DIVL...)
 *  These are stable inputs for immunum so unit tests can assert exact chain type
 *  and region coverage without fetching the CSV from the server. */
const HEAVY_1 = 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARVAPGALDYWGQGTLVTVSS';
const HEAVY_2 = 'EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDHYSGSGSYYYYFDYWGQGTLVTVSS';
const LIGHT_KAPPA = 'DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPRTFGQGTKVEIK';
const LIGHT_LAMBDA = 'QSALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTLVFGGGTKLTVL';

/** FR/CDR counts we expect in the immunum annotation JSON for IMGT/Kabat.
 *  The engine only accepts IMGT and Kabat — those are the choices declared in
 *  package.ts and surfaced in the dialog's scheme dropdown. */
const EXPECTED_REGION_COUNT = 7; // FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4
const EXPECTED_REGION_NAMES = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4'];

/** Sanity range for alignment confidence on canonical antibody sequences. */
const MIN_CONFIDENCE = 0.5;

category('antibody numbering (immunum)', () => {
  // Each numberSequencesWithImmunum call spawns a fresh worker and terminates
  // it before returning — no shared setup / teardown needed.

  test('worker: heavy chain (IMGT)', async () => {
    const [row] = await numberSequencesWithImmunum([HEAVY_1], 'imgt');
    expect(row.chainType, 'Heavy');
    expect(row.chainCode, 'H');
    expect(row.confidence >= MIN_CONFIDENCE, true);
    expect(row.positionNames.length > 0, true);
    expect(row.numberingDetail.length > 0, true);
    // numbering_map indices must fall inside the input sequence
    for (const idx of Object.values(row.numberingMap))
      expect(idx >= 0 && idx < HEAVY_1.length, true);
    expect(row.numberingDetail.length, Object.keys(row.numberingMap).length);
  });

  test('worker: light kappa chain (IMGT)', async () => {
    const [row] = await numberSequencesWithImmunum([LIGHT_KAPPA], 'imgt');
    expect(row.chainType, 'Light');
    expect(row.chainCode === 'K' || row.chainCode === 'L', true);
    expect(row.confidence >= MIN_CONFIDENCE, true);
    expect(row.numberingDetail.length > 0, true);
  });

  test('worker: light lambda chain (IMGT)', async () => {
    const [row] = await numberSequencesWithImmunum([LIGHT_LAMBDA], 'imgt');
    expect(row.chainType, 'Light');
    expect(row.confidence >= MIN_CONFIDENCE, true);
  });

  test('worker: batch numbering', async () => {
    const rows = await numberSequencesWithImmunum(
      [HEAVY_1, LIGHT_KAPPA, HEAVY_2, LIGHT_LAMBDA], 'imgt');
    expect(rows.length, 4);
    expect(rows[0].chainType, 'Heavy');
    expect(rows[1].chainType, 'Light');
    expect(rows[2].chainType, 'Heavy');
    expect(rows[3].chainType, 'Light');
  });

  test('worker: empty / short sequences fail gracefully', async () => {
    const rows = await numberSequencesWithImmunum(['', 'AAAA', '   '], 'imgt');
    expect(rows.length, 3);
    for (const r of rows) {
      expect(r.positionNames, '');
      expect(r.numberingDetail.length, 0);
      expect(r.error.length > 0, true);
    }
  });

  test('worker: kabat scheme returns kabat-style position codes', async () => {
    const [imgt] = await numberSequencesWithImmunum([HEAVY_1], 'imgt');
    const [kabat] = await numberSequencesWithImmunum([HEAVY_1], 'kabat');
    expect(imgt.chainType, 'Heavy');
    expect(kabat.chainType, 'Heavy');
    // Kabat numbering keys should not match IMGT one-for-one — the schemes
    // number the same residues differently. A weak but robust check: the set
    // of keys differs.
    const imgtKeys = new Set(Object.keys(imgt.numberingMap));
    const kabatKeys = new Set(Object.keys(kabat.numberingMap));
    let differ = false;
    for (const k of kabatKeys) if (!imgtKeys.has(k)) {differ = true; break;}
    expect(differ, true);
  });

  test('numberAntibodyColumn: DataFrame shape matches antpack script', async () => {
    const col = DG.Column.fromStrings('seq', [HEAVY_1, LIGHT_KAPPA, '']);
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    const result = await numberAntibodyColumn(col, 'imgt');

    // Required columns, same names as the Python script
    for (const name of ['position_names', 'chain_type', 'annotations_json',
      'numbering_detail', 'numbering_map']) {
      expect(result.col(name) !== null, true);
    }
    expect(result.rowCount, 3);

    // Row 0 — heavy chain: all 5 fields populated
    expect(result.get('position_names', 0).length > 0, true);
    expect(result.get('chain_type', 0), 'Heavy');
    const annot0 = JSON.parse(result.get('annotations_json', 0));
    expect(annot0.length, EXPECTED_REGION_COUNT);
    expect(annot0.map((a: any) => a.name).join(','), EXPECTED_REGION_NAMES.join(','));
    for (const a of annot0) {
      expect(a.visualType, 'region');
      expect(a.category, 'structure');
      expect(a.sourceScheme, 'IMGT');
      expect(a.autoGenerated, true);
    }

    // Row 1 — light chain: region JSON has same structure
    expect(result.get('chain_type', 1), 'Light');
    const annot1 = JSON.parse(result.get('annotations_json', 1));
    expect(annot1.length, EXPECTED_REGION_COUNT);

    // Row 2 — empty input: all fields blank / '[]'
    expect(result.get('position_names', 2), '');
    expect(result.get('chain_type', 2), '');
    expect(result.get('annotations_json', 2), '[]');
    expect(result.get('numbering_detail', 2), '');
    expect(result.get('numbering_map', 2), '');
  });

  test('numberAntibodyColumn: numbering_map indices line up with sequence', async () => {
    const col = DG.Column.fromStrings('seq', [HEAVY_1]);
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    const result = await numberAntibodyColumn(col, 'imgt');

    const detail = JSON.parse(result.get('numbering_detail', 0));
    const map = JSON.parse(result.get('numbering_map', 0));

    // For each numbered position: sequence[charIdx] must equal the recorded aa
    for (const entry of detail) {
      const idx = map[entry.position];
      expect(typeof idx === 'number', true);
      expect(HEAVY_1[idx], entry.aa);
    }
  });

  test('numberAntibodyColumn: annotations_json start/end resolve via numbering_map', async () => {
    const col = DG.Column.fromStrings('seq', [HEAVY_1]);
    col.semType = DG.SEMTYPE.MACROMOLECULE;
    const result = await numberAntibodyColumn(col, 'imgt');

    const annotations = JSON.parse(result.get('annotations_json', 0));
    const map = JSON.parse(result.get('numbering_map', 0));
    // FR1 start (position "1") must be present and resolve to a valid char index.
    const fr1 = annotations.find((a: any) => a.name === 'FR1');
    expect(fr1 !== undefined, true);
    const startIdx = map[fr1.start];
    expect(typeof startIdx === 'number', true);
    expect(startIdx >= 0 && startIdx < HEAVY_1.length, true);
  });

  test('numberAntibodyColumn: loads antibodies.csv sample subset', async () => {
    let df: DG.DataFrame;
    try {
      df = await _package.files.readCsv('samples/antibodies.csv');
    } catch (err) {
      // Sample may not be deployed on every server; skip instead of failing.
      console.warn('antibodies.csv not available — skipping', err);
      return;
    }
    const hcCol = df.col('AntibodyHC') ?? df.col('HeavyChain') ?? df.columns.byName('AntibodyHC');
    if (!hcCol) return;

    // Subset to the first 10 rows so the test finishes in seconds.
    const subset = DG.Column.fromStrings('seq',
      Array.from({length: Math.min(10, hcCol.length)}, (_, i) => hcCol.get(i) ?? ''));
    subset.semType = DG.SEMTYPE.MACROMOLECULE;

    const result = await numberAntibodyColumn(subset, 'imgt');
    expect(result.rowCount, subset.length);

    let heavyCount = 0;
    for (let i = 0; i < result.rowCount; i++)
      if (result.get('chain_type', i) === 'Heavy') heavyCount++;
    // Expect the majority of the HC column to be classified as heavy.
    expect(heavyCount >= Math.ceil(subset.length * 0.6), true);
  });
});
