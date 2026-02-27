/* eslint-disable max-lines-per-function */
/* eslint-disable max-len */
/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {smi} from './randsmiles';
// Placeholder imports for glyph PNGs (ensure your bundler loads them as base64 strings)
import {glyphPool} from './glyphs/glyphs';
import {DBExplorer} from '@datagrok-libraries/db-explorer/src/db-explorer';
import {moleculeRenderer, imageRenderer, rawImageRenderer, helmRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
import {biologicsConfig} from './config';
import {helms} from './helms';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


//name: createDemoBiologicsData
export async function createDemoBiologicsData(): Promise<any> {
  const cn = 'Biologics:biologics'; // connection name
  const AA = 'ACDEFGHIKLMNPQRSTVWYGGGGGGAAAADEE';
  const randInt = (min: number, max: number) => Math.floor(min + Math.random() * (max - min + 1));
  const randPick = <T>(arr: T[]) => arr[randInt(0, arr.length - 1)];
  const escape = (s: string) => s.replace(/'/g, '\'\'').replaceAll('@', '');

  const exec = async (sql: string) => await grok.data.db.query(cn, sql);

  // 0. Clear existing demo data (keep assay types & target organisms)
  await exec(`
    DELETE FROM biologics.assay_results;
    DELETE FROM biologics.sequence_liabilities;
    DELETE FROM biologics.sequence_properties;
    DELETE FROM biologics.sequence_regions;
    DELETE FROM biologics.adc;
    DELETE FROM biologics.expression_batches;
    DELETE FROM biologics.purification_batches;
    DELETE FROM biologics.linkers;
    DELETE FROM biologics.drugs;
    DELETE FROM biologics.sequences;
    DELETE FROM biologics.peptides;
  `); // cascade should handle the rest

  const getConservedRegion = () => {
    const conservedRegionSeq = 'CYSASNITLYCVHQRGGGGAAAAGGGGGGGSSSVTVS';
    let actCons = conservedRegionSeq.substring(0, Math.max(Math.floor(Math.random() * (conservedRegionSeq.length - 1)), Math.floor(conservedRegionSeq.length / 2)));
    const substitutions = randInt(1, 4);
    for (let i = 0; i < substitutions; i++) {
      const pos = randInt(0, actCons.length - 1);
      const newAA = randPick(AA as unknown as string[]);
      actCons = actCons.substring(0, pos) + newAA + actCons.substring(pos + 1);
    }
    return actCons;
  };


  // 0. insert helms into peptides table
  const insertedHelms: {id: number, helm: string, name: string}[] = [];
  for (let i = 0; i < helms.length; i++) {
    const helm = helms[i];
    const name = `Peptide_${i + 1}`;
    const df = await exec(`
      INSERT INTO biologics.peptides(name, helm)
      VALUES ('${escape(name)}', '${escape(helm)}')
      RETURNING id
    `);
    insertedHelms.push({id: df.get('id', 0), helm, name});
  }


  // 1. Generate 500 protein sequences (~1000 length)
  const seqCount = 200;
  const sequences: { name: string, seq: string }[] = [];
  for (let i = 1; i <= seqCount; i++) {
    const len = 950 + randInt(0, 100); // 950-1050
    let s = '';
    for (let j = 0; j < len; j++) {
      s += AA[randInt(0, AA.length - 1)];
      if (j % 150 == 1) {
        s += getConservedRegion();
        j = s.length + 1;
      }
    }
    sequences.push({name: `Sequence_${i}`, seq: s});
  }

  async function insertSequences(list: { name: string, seq: string }[]) {
    const chunkSize = 20;
    const ids: number[] = [];
    for (let i = 0; i < list.length; i += chunkSize) {
      const chunk = list.slice(i, i + chunkSize);
      const values = chunk
        .map((c) => `('PROTEIN','${escape(c.seq)}','${escape(c.name)}')`)
        .join(',');
      const sql = `
        INSERT INTO biologics.sequences(sequence_type, sequence, name)
        VALUES ${values}
        RETURNING id, name
      `;
      const df = await exec(sql);
      for (let r = 0; r < df.rowCount; r++)
        ids.push(df.get('id', r));
    }
    return ids;
  }

  // 1b. Insert annotated regions for each sequence
  async function insertSequenceRegions(seqIds: number[], seqs: { name: string, seq: string }[]) {
    const regionChunks: string[] = [];
    for (let i = 0; i < seqIds.length; i++) {
      const seqId = seqIds[i];
      const seqLen = seqs[i].seq.length;
      // Antibody domain boundaries (with small random jitter)
      const jit = () => randInt(-3, 3);
      const hcEnd = Math.min(450 + jit(), seqLen);
      const lcStart = hcEnd + 1;
      const lcEnd = Math.min(lcStart + 219 + jit(), seqLen);
      // Heavy chain domains
      const vhEnd = Math.min(120 + jit(), hcEnd);
      const ch1Start = vhEnd + 1; const ch1End = Math.min(ch1Start + 109 + jit(), hcEnd);
      const hingeStart = ch1End + 1; const hingeEnd = Math.min(hingeStart + 14 + jit(), hcEnd);
      const ch2Start = hingeEnd + 1; const ch2End = Math.min(ch2Start + 109 + jit(), hcEnd);
      const ch3Start = ch2End + 1; const ch3End = Math.min(ch3Start + 109 + jit(), hcEnd);
      // Light chain domains
      const vlEnd = Math.min(lcStart + 109 + jit(), lcEnd);
      const clStart = vlEnd + 1; const clEnd = lcEnd;
      // CDRs in VH (Kabat-like numbering with jitter)
      const cdrH1s = 26 + jit(); const cdrH1e = Math.min(cdrH1s + 9 + randInt(0, 3), vhEnd);
      const cdrH2s = 50 + jit(); const cdrH2e = Math.min(cdrH2s + 15 + randInt(0, 4), vhEnd);
      const cdrH3s = 95 + jit(); const cdrH3e = Math.min(cdrH3s + 8 + randInt(0, 15), vhEnd);
      // CDRs in VL
      const cdrL1s = lcStart + 23 + jit(); const cdrL1e = Math.min(cdrL1s + 12 + randInt(0, 5), vlEnd);
      const cdrL2s = lcStart + 49 + jit(); const cdrL2e = Math.min(cdrL2s + 6 + randInt(0, 2), vlEnd);
      const cdrL3s = lcStart + 89 + jit(); const cdrL3e = Math.min(cdrL3s + 8 + randInt(0, 3), vlEnd);

      const regions: [string, number, number][] = [
        ['HC', 1, hcEnd], ['LC', lcStart, lcEnd],
        ['VH', 1, vhEnd], ['VL', lcStart, vlEnd],
        ['CH1', ch1Start, ch1End], ['Hinge', hingeStart, hingeEnd],
        ['CH2', ch2Start, ch2End], ['CH3', ch3Start, ch3End],
        ['CL', clStart, clEnd],
        ['CDR-H1', cdrH1s, cdrH1e], ['CDR-H2', cdrH2s, cdrH2e], ['CDR-H3', cdrH3s, cdrH3e],
        ['CDR-L1', cdrL1s, cdrL1e], ['CDR-L2', cdrL2s, cdrL2e], ['CDR-L3', cdrL3s, cdrL3e],
      ];
      for (const [type, s, e] of regions) {
        if (s >= 1 && e >= s && e <= seqLen)
          regionChunks.push(`(${seqId}, '${type}', ${s}, ${e})`);
      }
    }
    const chunkSize = 200;
    for (let i = 0; i < regionChunks.length; i += chunkSize) {
      const part = regionChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_regions(sequence_id, region_type, start_pos, end_pos) VALUES ${part.join(',')}`);
    }
    return regionChunks.length;
  }

  // 1c. Insert computed biophysical properties per sequence
  async function insertSequenceProperties(seqIds: number[], seqs: { name: string, seq: string }[]) {
    const propChunks: string[] = [];
    for (let i = 0; i < seqIds.length; i++) {
      const seqLen = seqs[i].seq.length;
      const mw = +(seqLen * 110.0 + randInt(-500, 500)).toFixed(1); // Da, ~110 per AA
      const pI = +(6.0 + Math.random() * 3.0).toFixed(2); // 6.0-9.0
      const gravy = +(-0.5 + Math.random() * 1.0).toFixed(3); // -0.5 to 0.5
      const charge = +(-10 + Math.random() * 20).toFixed(1); // -10 to +10
      const aggProp = +(Math.random()).toFixed(3); // 0 to 1
      propChunks.push(`(${seqIds[i]}, ${mw}, ${pI}, ${gravy}, ${charge}, ${aggProp})`);
    }
    const chunkSize = 100;
    for (let i = 0; i < propChunks.length; i += chunkSize) {
      const part = propChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_properties(sequence_id, molecular_weight, isoelectric_point, hydrophobicity_index, charge_at_ph7, aggregation_propensity) VALUES ${part.join(',')}`);
    }
  }

  // 1d. Insert developability liabilities per sequence
  async function insertSequenceLiabilities(seqIds: number[], seqs: { name: string, seq: string }[]) {
    const liabChunks: string[] = [];
    const risks = ['Low', 'Medium', 'High'];
    for (let i = 0; i < seqIds.length; i++) {
      const seq = seqs[i].seq;
      // Scan for common liability motifs
      for (let p = 0; p < seq.length - 1; p++) {
        const di = seq.substring(p, p + 2);
        if ((di === 'NG' || di === 'NS') && Math.random() < 0.5)
          liabChunks.push(`(${seqIds[i]}, 'Deamidation', ${p + 1}, '${di}', '${randPick(risks)}')`);
        else if ((di === 'DG' || di === 'DS') && Math.random() < 0.3)
          liabChunks.push(`(${seqIds[i]}, 'Isomerization', ${p + 1}, '${di}', '${randPick(risks)}')`);
        if (seq[p] === 'M' && Math.random() < 0.2)
          liabChunks.push(`(${seqIds[i]}, 'Oxidation', ${p + 1}, 'Met', '${randPick(risks)}')`);
      }
    }
    const chunkSize = 200;
    for (let i = 0; i < liabChunks.length; i += chunkSize) {
      const part = liabChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_liabilities(sequence_id, liability_type, position, motif, risk_level) VALUES ${part.join(',')}`);
    }
    return liabChunks.length;
  }

  // 2. Drugs (empty smiles list for user to fill later)
  const drugSmiles: { name: string, smiles: string }[] = [
    // { name: 'Drug_1', smiles: 'CCO...' },
  ];

  for (let i = 1; i <= 200; i++)
    drugSmiles.push({name: `Drug_${i}`, smiles: smi[i]});


  async function insertDrugs(list: { name: string, smiles: string }[]) {
    if (list.length === 0) return [];
    const values = list.map((d) => `('${escape(d.name)}','${escape(d.smiles)}')`).join(',');
    const sql = `
      INSERT INTO biologics.drugs(name, smiles)
      VALUES ${values}
      RETURNING id
    `;
    const df = await exec(sql);
    const ids: number[] = [];
    for (let r = 0; r < df.rowCount; r++)
      ids.push(df.get('id', r));
    return ids;
  }

  // 3. Linkers (mix of SMALL & PROTEIN)
  const linkerProteinSeqs = ['GGGGS', 'GGSGGS', 'GGGSGGGS', 'GSGSG', 'GPGPG'];
  const linkerSmallSmiles = ['CCOCCOCCOCCOCCOCCOCCOCCO', 'CCNCCNCCNCCNCCNCCNCCN', 'CNC(=O)CCNC(=O)CCNC(=O)CCNC(=O)C', 'NCC(=O)NCC(=O)NCC(=O)NCC=O', 'COCOCOCOCOCOCOCO'];
  const linkers: {
    linker_type: 'SMALL' | 'PROTEIN',
    linker_sequence?: string,
    linker_molecule_smiles?: string
  }[] = [];

  for (let i = 0; i < 5; i++)
    linkers.push({linker_type: 'PROTEIN', linker_sequence: linkerProteinSeqs[i]});
  for (let i = 0; i < 5; i++)
    linkers.push({linker_type: 'SMALL', linker_molecule_smiles: linkerSmallSmiles[i]});

  async function insertLinkers(list: typeof linkers) {
    if (!list.length) return [];
    const values = list.map((l) => {
      if (l.linker_type === 'PROTEIN')
        return `('PROTEIN', NULL, '${escape(l.linker_sequence!)}')`;
      else
        return `('SMALL', '${escape(l.linker_molecule_smiles!)}', NULL)`;
    }).join(',');
    const sql = `
      INSERT INTO biologics.linkers(linker_type, linker_molecule_smiles, linker_sequence)
      VALUES ${values}
      RETURNING id
    `;
    const df = await exec(sql);
    const ids: number[] = [];
    for (let r = 0; r < df.rowCount; r++)
      ids.push(df.get('id', r));
    return ids;
  }

  // 4. Fetch existing assay_types, target organisms, and targets
  const assayTypesDf = await exec('SELECT id, name, category FROM biologics.assay_types');
  const assayTypes: { id: number, name: string, category: string }[] = [];
  for (let r = 0; r < assayTypesDf.rowCount; r++)
    assayTypes.push({id: assayTypesDf.get('id', r), name: assayTypesDf.get('name', r), category: assayTypesDf.get('category', r) || 'Legacy'});

  const newPanelTypes = assayTypes.filter((at) => at.category !== 'Legacy');
  const legacyTypes = assayTypes.filter((at) => at.category === 'Legacy');

  const orgDf = await exec('SELECT id, name FROM biologics.target_organisms');
  const organisms: number[] = [];
  for (let r = 0; r < orgDf.rowCount; r++)
    organisms.push(orgDf.get('id', r));

  const targetsDf = await exec('SELECT id, name FROM biologics.targets');
  const targets: number[] = [];
  for (let r = 0; r < targetsDf.rowCount; r++)
    targets.push(targetsDf.get('id', r));

  // 5. Insert everything
  const sequenceIds = await insertSequences(sequences);
  const drugIds = await insertDrugs(drugSmiles);
  const linkerIds = await insertLinkers(linkers);

  // 5b. Insert sequence annotations
  const regionCount = await insertSequenceRegions(sequenceIds, sequences);
  await insertSequenceProperties(sequenceIds, sequences);
  const liabilityCount = await insertSequenceLiabilities(sequenceIds, sequences);

  // 6. Purification batches (random subset)
  const purBatchesValues: string[] = [];
  new Array(sequenceIds.length * 3).fill(null).map((_, i) => randPick(sequenceIds)).forEach((id) => {
    purBatchesValues.push(`(${id}, 'PurBatch_${id}_${Math.floor(Math.random() * 1000)}', 'Auto-generated purification batch')`);
  });
  if (purBatchesValues.length) {
    await exec(`
      INSERT INTO biologics.purification_batches(sequence_id, name, notes)
      VALUES ${purBatchesValues.join(',')}
    `);
  }

  // 7. Expression batches (random subset)
  const exprSystems = ['CHO', 'HEK293', 'E.coli'];
  const exprValues: string[] = [];
  new Array(sequenceIds.length * 3).fill(null).map((_, i) => randPick(sequenceIds)).forEach((id) => {
    exprValues.push(`(${id}, '${escape(randPick(exprSystems))}', ${ (Math.random() * 150).toFixed(2) }, 'ExprBatch_${id}_${(Math.floor(Math.random() * 1000))}', 'Auto-generated expression batch')`);
  });
  if (exprValues.length) {
    await exec(`
      INSERT INTO biologics.expression_batches(sequence_id, expression_system, yield_mg, name, notes)
      VALUES ${exprValues.join(',')}
    `);
  }

  // 8. ADCs (only if we have drugs)
  let adcCount = 0;
  const adcIds: number[] = [];
  if (drugIds.length && linkerIds.length) {
    const adcValues: string[] = [];
    const adcToMake = Math.min(drugIds.length * linkerIds.length, sequenceIds.length * 3);
    for (let i = 0; i < adcToMake; i++) {
      const antibodyId = randPick(sequenceIds);
      const drugId = randPick(drugIds);
      const linkerId = randPick(linkerIds);
      const name = `ADC_${i + 1}`;
      const glyph = ''; // placeholder (base64 PNG string)
      adcValues.push(`('${escape(name)}', ${antibodyId}, ${linkerId}, ${drugId}, ${glyph === '' ? 'NULL' : `'${escape(glyph)}'`})`);
    }
    if (adcValues.length) {
      const df = await exec(`
        INSERT INTO biologics.adc(name, antibody_id, linker_id, drug_id, glyph)
        VALUES ${adcValues.join(',')}
        RETURNING id
      `);
      for (let r = 0; r < df.rowCount; r++)
        adcIds.push(df.get('id', r));
      adcCount = adcValues.length;
    }
  }

  // 9. Assay results (random)
  function randomAssayValue(name: string): { val: number, units: string } {
    // New panel: Primary Screening
    if (/ELISA.*OD/i.test(name)) return {val: +(0.1 + Math.random() * 2.9).toFixed(3), units: 'OD'};
    if (/ELISA.*EC50/i.test(name)) return {val: +(0.1 + Math.random() * 99.9).toFixed(2), units: 'nM'};
    if (/Flow.*Positive/i.test(name)) return {val: +(10 + Math.random() * 89).toFixed(1), units: '%'};
    if (/Flow.*EC50/i.test(name)) return {val: +(0.5 + Math.random() * 499.5).toFixed(2), units: 'nM'};
    // Epitope Binning
    if (/SPR.*Inhibition/i.test(name)) return {val: +(5 + Math.random() * 90).toFixed(1), units: '%'};
    // Affinity Characterization (exact match to avoid cross-matching)
    if (name === 'SPR ka') return {val: +(1e4 + Math.random() * 9.9e5).toFixed(0), units: '1/Ms'};
    if (name === 'SPR kd') return {val: +(1e-4 + Math.random() * 9.9e-3).toFixed(5), units: '1/s'};
    if (name === 'SPR KD') return {val: +(0.1 + Math.random() * 99.9).toFixed(2), units: 'nM'};
    // Functional Assays
    if (/Reporter.*EC50/i.test(name)) return {val: +(0.1 + Math.random() * 499.9).toFixed(2), units: 'nM'};
    if (/Reporter.*Emax/i.test(name)) return {val: +(50 + Math.random() * 100).toFixed(1), units: '%'};
    // Developability
    if (/Thermal.*Tm/i.test(name)) return {val: +(55 + Math.random() * 30).toFixed(1), units: 'C'};
    if (/Monomer/i.test(name)) return {val: +(85 + Math.random() * 14.9).toFixed(1), units: '%'};
    if (/Aggregate Size/i.test(name)) return {val: +(5 + Math.random() * 45).toFixed(1), units: 'nm'};
    if (/Viscosity/i.test(name)) return {val: +(1 + Math.random() * 49).toFixed(1), units: 'cP'};
    // Legacy types
    if (/IC50/i.test(name)) return {val: +(Math.random() * 500).toFixed(2), units: 'nM'};
    if (/Caspase/i.test(name)) return {val: +(Math.random() * 10).toFixed(2), units: 'RFU'};
    if (/half-life/i.test(name)) return {val: +(Math.random() * 240).toFixed(1), units: 'min'};
    if (/Binding affinity/i.test(name)) return {val: +(Math.random() * 50).toFixed(2), units: 'nM'};
    if (/Cell binding/i.test(name)) return {val: +(Math.random() * 100).toFixed(2), units: 'nM'};
    if (/DAR/i.test(name)) return {val: +(2 + Math.random() * 6).toFixed(2), units: 'ratio'};
    if (/Cmax/i.test(name)) return {val: +(Math.random() * 100).toFixed(2), units: 'ug/mL'};
    if (/Tmax/i.test(name)) return {val: +(Math.random() * 48).toFixed(2), units: 'h'};
    if (/AUC/i.test(name)) return {val: +(Math.random() * 5000).toFixed(1), units: 'ug*h/mL'};
    return {val: +(Math.random() * 100).toFixed(2), units: ''};
  }

  const assayResultChunks: string[] = [];
  const maxResults = 5000; // cap
  let inserted = 0;
  for (const seqId of sequenceIds) {
    const perSeqAssays = randInt(3, 26);
    for (let i = 0; i < perSeqAssays; i++) {
      const at = randPick(legacyTypes.length ? legacyTypes : assayTypes);
      const org = randPick(organisms);
      const tgt = randPick(targets);
      const av = randomAssayValue(at.name);
      assayResultChunks.push(`(${at.id}, ${randPick(adcIds)}, ${org}, ${av.val}, '${escape(av.units)}', ${tgt})`);
      inserted++;
      if (inserted >= maxResults) break;
    }
    if (inserted >= maxResults) break;
  }
  if (assayResultChunks.length) {
    const chunkSize = 200;
    for (let i = 0; i < assayResultChunks.length; i += chunkSize) {
      const part = assayResultChunks.slice(i, i + chunkSize);
      await exec(`
        INSERT INTO biologics.assay_results(assay_id, adc_id, target_organism_id, result_value, units, target_id)
        VALUES ${part.join(',')}
      `);
    }
  }
  // 10. add peptide assay results
  const peptideAssayChunks: string[] = [];
  const maxPeptideResults = 2000;
  let peptideAssayInserted = 0;
  for (const peptide of insertedHelms) {
    const perPeptideAssays = randInt(3, 10);
    for (let i = 0; i < perPeptideAssays; i++) {
      const at = randPick(assayTypes);
      const org = randPick(organisms);
      const tgt = randPick(targets);
      const av = randomAssayValue(at.name);
      peptideAssayChunks.push(`(${at.id}, ${peptide.id}, ${org}, ${av.val}, '${escape(av.units)}', ${tgt})`);
      peptideAssayInserted++;
      if (peptideAssayInserted >= maxPeptideResults) break;
    }
    if (peptideAssayInserted >= maxPeptideResults) break;
  }
  if (peptideAssayChunks.length) {
    const chunkSize = 200;
    for (let i = 0; i < peptideAssayChunks.length; i += chunkSize) {
      const part = peptideAssayChunks.slice(i, i + chunkSize);
      await exec(`
        INSERT INTO biologics.assay_results(assay_id, peptide_id, target_organism_id, result_value, units, target_id)
        VALUES ${part.join(',')}
      `);
    }
  }

  // 11. Comprehensive assay panel: every ADC gets one result per new-panel assay type
  let comprehensiveInserted = 0;
  if (newPanelTypes.length && adcIds.length) {
    const compChunks: string[] = [];
    for (const adcId of adcIds) {
      const tgt = randPick(targets);
      const org = randPick(organisms);
      for (const at of newPanelTypes) {
        const av = randomAssayValue(at.name);
        compChunks.push(`(${at.id}, ${adcId}, ${org}, ${av.val}, '${escape(av.units)}', ${tgt})`);
      }
    }
    comprehensiveInserted = compChunks.length;
    const chunkSize = 200;
    for (let i = 0; i < compChunks.length; i += chunkSize) {
      const part = compChunks.slice(i, i + chunkSize);
      await exec(`
        INSERT INTO biologics.assay_results(assay_id, adc_id, target_organism_id, result_value, units, target_id)
        VALUES ${part.join(',')}
      `);
    }
  }

  // 12. Backfill sequence_id on assay_results from ADC antibody_id
  await exec(`
    UPDATE biologics.assay_results ar
    SET sequence_id = adc.antibody_id
    FROM biologics.adc adc
    WHERE ar.adc_id = adc.id AND ar.sequence_id IS NULL
  `);

  await populateAdcGlyphs(adcCount + 1);

  return {
    sequences: sequenceIds.length,
    drugs: drugIds.length,
    linkers: linkerIds.length,
    assay_results_legacy: inserted,
    assay_results_comprehensive: comprehensiveInserted,
    purification_batches: purBatchesValues.length,
    expression_batches: exprValues.length,
    adcs: adcCount,
    peptides: insertedHelms.length,
    sequence_regions: regionCount,
    sequence_liabilities: liabilityCount,
  };
}

//name: populateAdcGlyphs
//description: Populates missing ADC glyphs with random PNG (base64) strings
export async function populateAdcGlyphs(limit: number = 50): Promise<{updated: number}> {
  const cn = 'Biologics:biologics';
  const exec = async (sql: string) => await grok.data.db.query(cn, sql);
  const escape = (s: string) => s.replace(/'/g, '\'\'');

  // Use real imported glyph images if available, else fall back
  if (!glyphPool.length)
    throw new Error('No glyph images available');

  // Fetch ADC ids without glyphs (NULL or empty)
  const df = await exec(`SELECT id FROM biologics.adc WHERE (glyph IS NULL OR glyph='') LIMIT ${limit}`);
  if (df.rowCount === 0)
    return {updated: 0};

  const updates: string[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    const id = df.get('id', i);
    const glyph = glyphPool[Math.floor(Math.random() * glyphPool.length)];
    updates.push(`UPDATE biologics.adc SET glyph='${escape(glyph)}' WHERE id=${id}`);
  }
  const batchSize = 20;
  for (let i = 0; i < updates.length; i += batchSize)
    await exec(updates.slice(i, (i + batchSize)).join(';') + ';');
  return {updated: updates.length};
}

//name: autostartbiologics
//meta.role: autostart
export async function autostartBiologics() {
  const exp = DBExplorer.initFromConfig(biologicsConfig);
  if (!exp) {
    grok.shell.error('Failed to initialize Biologics DB Explorer');
    return;
  }
  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'structure' || lc.includes('smiles') || lc.includes('compound_structure')) && typeof value === 'string' && grok.chem.checkSmiles(value);
  }, (value) => moleculeRenderer(value as string));

  exp.addCustomRenderer((_, colName, value) => colName?.toLowerCase()?.includes('helm') && typeof value === 'string' && value?.toLowerCase()?.startsWith('peptide'),
    (value) => helmRenderer(value as string)
  );

  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'glyph' || lc === 'image' || lc === 'png' || lc === 'thumbnail') && typeof value === 'string' && value.startsWith('iVBORw0KGgo');
  }, (value) => rawImageRenderer(value as string));

  // exp.addDefaultHeaderReplacerColumns(['units']);
  console.log('Biologics object handlers registered');
}
