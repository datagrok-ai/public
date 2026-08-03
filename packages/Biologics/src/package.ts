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
  const randInt = (min: number, max: number) => Math.floor(min + Math.random() * (max - min + 1));
  const randPick = <T>(arr: T[]) => arr[randInt(0, arr.length - 1)];
  const escape = (s: string) => s.replace(/'/g, '\'\'').replaceAll('@', '');

  const exec = async (sql: string) => await grok.data.db.query(cn, sql);

  // 0. Clear existing demo data (keep assay types & target organisms)
  await exec(`
    DELETE FROM biologics.assay_curves;
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

  // Load real antibody chain data (heavy + light) from the bundled CSV.
  const antibodiesDf = await _package.files.readCsv('antibodies.csv');
  const hcCol = antibodiesDf.col('AntibodyHC')!;
  const lcCol = antibodiesDf.col('AntibodyLC')!;
  const idCol = antibodiesDf.col('Antibody_ID');

  // Load curve pool used to populate assay_curves for dose-response style assays.
  const curvesDf = await _package.files.readCsv('random_curves.csv');
  const curveCol = curvesDf.col('curve')!;
  const curvePool: string[] = [];
  for (let i = 0; i < curvesDf.rowCount; i++) {
    const c = curveCol.get(i);
    if (typeof c === 'string' && c.length > 0)
      curvePool.push(c);
  }


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


  // 1. Pick up to 400 real antibody rows (heavy + light chain) from the bundled CSV.
  const seqCount = Math.min(400, antibodiesDf.rowCount);
  const sequences: { name: string, heavy: string, light: string }[] = [];
  for (let i = 0; i < seqCount; i++) {
    const heavy = (hcCol.get(i) as string) ?? '';
    const light = (lcCol.get(i) as string) ?? '';
    if (!heavy && !light) continue;
    const baseName = idCol ? (idCol.get(i) as string) : '';
    const name = baseName && baseName.length ? `Ab_${baseName}` : `Antibody_${i + 1}`;
    sequences.push({name, heavy, light});
  }

  async function insertSequences(list: { name: string, heavy: string, light: string }[]) {
    const chunkSize = 20;
    const ids: number[] = [];
    for (let i = 0; i < list.length; i += chunkSize) {
      const chunk = list.slice(i, i + chunkSize);
      const values = chunk
        .map((c) => `('PROTEIN','${escape(c.heavy)}','${escape(c.light)}','${escape(c.name)}')`)
        .join(',');
      const sql = `
        INSERT INTO biologics.sequences(sequence_type, heavy_chain, light_chain, name)
        VALUES ${values}
        RETURNING id, name
      `;
      const df = await exec(sql);
      for (let r = 0; r < df.rowCount; r++)
        ids.push(df.get('id', r));
    }
    return ids;
  }

  // 1b. Insert annotated regions per chain. Positions are within the heavy_chain or light_chain.
  async function insertSequenceRegions(seqIds: number[], seqs: { name: string, heavy: string, light: string }[]) {
    const regionChunks: string[] = [];
    const jit = () => randInt(-3, 3);

    for (let i = 0; i < seqIds.length; i++) {
      const seqId = seqIds[i];
      const hcLen = seqs[i].heavy.length;
      const lcLen = seqs[i].light.length;

      // Heavy chain regions (positions relative to heavy_chain)
      if (hcLen > 0) {
        const vhEnd = Math.min(120 + jit(), hcLen);
        const ch1Start = vhEnd + 1; const ch1End = Math.min(ch1Start + 109 + jit(), hcLen);
        const hingeStart = ch1End + 1; const hingeEnd = Math.min(hingeStart + 14 + jit(), hcLen);
        const ch2Start = hingeEnd + 1; const ch2End = Math.min(ch2Start + 109 + jit(), hcLen);
        const ch3Start = ch2End + 1; const ch3End = Math.min(ch3Start + 109 + jit(), hcLen);
        const cdrH1s = 26 + jit(); const cdrH1e = Math.min(cdrH1s + 9 + randInt(0, 3), vhEnd);
        const cdrH2s = 50 + jit(); const cdrH2e = Math.min(cdrH2s + 15 + randInt(0, 4), vhEnd);
        const cdrH3s = 95 + jit(); const cdrH3e = Math.min(cdrH3s + 8 + randInt(0, 15), vhEnd);

        const hcRegions: [string, number, number][] = [['VH', 1, vhEnd]];
        if (ch1Start <= hcLen) hcRegions.push(['CH1', ch1Start, ch1End]);
        if (hingeStart <= hcLen) hcRegions.push(['Hinge', hingeStart, hingeEnd]);
        if (ch2Start <= hcLen) hcRegions.push(['CH2', ch2Start, ch2End]);
        if (ch3Start <= hcLen) hcRegions.push(['CH3', ch3Start, ch3End]);
        hcRegions.push(['CDR-H1', cdrH1s, cdrH1e]);
        hcRegions.push(['CDR-H2', cdrH2s, cdrH2e]);
        hcRegions.push(['CDR-H3', cdrH3s, cdrH3e]);

        for (const [type, s, e] of hcRegions) {
          if (s >= 1 && e >= s && e <= hcLen)
            regionChunks.push(`(${seqId}, '${type}', 'HC', ${s}, ${e})`);
        }
      }

      // Light chain regions (positions relative to light_chain)
      if (lcLen > 0) {
        const vlEnd = Math.min(109 + jit(), lcLen);
        const clStart = vlEnd + 1; const clEnd = lcLen;
        const cdrL1s = 23 + jit(); const cdrL1e = Math.min(cdrL1s + 12 + randInt(0, 5), vlEnd);
        const cdrL2s = 49 + jit(); const cdrL2e = Math.min(cdrL2s + 6 + randInt(0, 2), vlEnd);
        const cdrL3s = 89 + jit(); const cdrL3e = Math.min(cdrL3s + 8 + randInt(0, 3), vlEnd);

        const lcRegions: [string, number, number][] = [['VL', 1, vlEnd]];
        if (clStart <= lcLen) lcRegions.push(['CL', clStart, clEnd]);
        lcRegions.push(['CDR-L1', cdrL1s, cdrL1e]);
        lcRegions.push(['CDR-L2', cdrL2s, cdrL2e]);
        lcRegions.push(['CDR-L3', cdrL3s, cdrL3e]);

        for (const [type, s, e] of lcRegions) {
          if (s >= 1 && e >= s && e <= lcLen)
            regionChunks.push(`(${seqId}, '${type}', 'LC', ${s}, ${e})`);
        }
      }
    }

    const chunkSize = 200;
    for (let i = 0; i < regionChunks.length; i += chunkSize) {
      const part = regionChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_regions(sequence_id, region_type, chain, start_pos, end_pos) VALUES ${part.join(',')}`);
    }
    return regionChunks.length;
  }

  // 1c. Insert computed biophysical properties per chain (HC and LC each get their own row)
  async function insertSequenceProperties(seqIds: number[], seqs: { name: string, heavy: string, light: string }[]) {
    const propChunks: string[] = [];
    for (let i = 0; i < seqIds.length; i++) {
      for (const ch of [{tag: 'HC', seq: seqs[i].heavy}, {tag: 'LC', seq: seqs[i].light}]) {
        if (!ch.seq.length) continue;
        const seqLen = ch.seq.length;
        const mw = +(seqLen * 110.0 + randInt(-500, 500)).toFixed(1); // Da, ~110 per AA
        const pI = +(6.0 + Math.random() * 3.0).toFixed(2); // 6.0-9.0
        const gravy = +(-0.5 + Math.random() * 1.0).toFixed(3); // -0.5 to 0.5
        const charge = +(-10 + Math.random() * 20).toFixed(1); // -10 to +10
        const aggProp = +(Math.random()).toFixed(3); // 0 to 1
        propChunks.push(`(${seqIds[i]}, '${ch.tag}', ${mw}, ${pI}, ${gravy}, ${charge}, ${aggProp})`);
      }
    }
    const chunkSize = 100;
    for (let i = 0; i < propChunks.length; i += chunkSize) {
      const part = propChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_properties(sequence_id, chain, molecular_weight, isoelectric_point, hydrophobicity_index, charge_at_ph7, aggregation_propensity) VALUES ${part.join(',')}`);
    }
  }

  // 1d. Insert developability liabilities; record the chain so the position is unambiguous
  async function insertSequenceLiabilities(seqIds: number[], seqs: { name: string, heavy: string, light: string }[]) {
    const liabChunks: string[] = [];
    const risks = ['Low', 'Medium', 'High'];
    for (let i = 0; i < seqIds.length; i++) {
      for (const ch of [{tag: 'HC', seq: seqs[i].heavy}, {tag: 'LC', seq: seqs[i].light}]) {
        const seq = ch.seq;
        for (let p = 0; p < seq.length - 1; p++) {
          const di = seq.substring(p, p + 2);
          if ((di === 'NG' || di === 'NS') && Math.random() < 0.5)
            liabChunks.push(`(${seqIds[i]}, '${ch.tag}', 'Deamidation', ${p + 1}, '${di}', '${randPick(risks)}')`);
          else if ((di === 'DG' || di === 'DS') && Math.random() < 0.3)
            liabChunks.push(`(${seqIds[i]}, '${ch.tag}', 'Isomerization', ${p + 1}, '${di}', '${randPick(risks)}')`);
          if (seq[p] === 'M' && Math.random() < 0.2)
            liabChunks.push(`(${seqIds[i]}, '${ch.tag}', 'Oxidation', ${p + 1}, 'Met', '${randPick(risks)}')`);
        }
      }
    }
    const chunkSize = 200;
    for (let i = 0; i < liabChunks.length; i += chunkSize) {
      const part = liabChunks.slice(i, i + chunkSize);
      await exec(`INSERT INTO biologics.sequence_liabilities(sequence_id, chain, liability_type, position, motif, risk_level) VALUES ${part.join(',')}`);
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
  const assayTypesDf = await exec('SELECT id, name, category, fit_model FROM biologics.assay_types');
  const assayTypes: { id: number, name: string, category: string, fit_model: string | null }[] = [];
  for (let r = 0; r < assayTypesDf.rowCount; r++) {
    assayTypes.push({
      id: assayTypesDf.get('id', r),
      name: assayTypesDf.get('name', r),
      category: assayTypesDf.get('category', r) || 'Legacy',
      fit_model: assayTypesDf.get('fit_model', r) || null,
    });
  }

  const newPanelTypes = assayTypes.filter((at) => at.category !== 'Legacy');
  const legacyTypes = assayTypes.filter((at) => at.category === 'Legacy');
  const curveEligibleAssayIds = new Set(
    assayTypes.filter((at) => !!at.fit_model && curvePool.length > 0).map((at) => at.id)
  );

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

  // 6. Purification batches — heavy and light chains are purified separately, then assembled.
  const chains = ['HC', 'LC'] as const;
  const purBatchesValues: string[] = [];
  new Array(sequenceIds.length * 3).fill(null).map(() => randPick(sequenceIds)).forEach((id) => {
    const ch = randPick(chains as unknown as string[]);
    purBatchesValues.push(`(${id}, '${ch}', 'PurBatch_${id}_${ch}_${Math.floor(Math.random() * 1000)}', 'Auto-generated purification batch (${ch})')`);
  });
  if (purBatchesValues.length) {
    await exec(`
      INSERT INTO biologics.purification_batches(sequence_id, chain, name, notes)
      VALUES ${purBatchesValues.join(',')}
    `);
  }

  // 7. Expression batches — same per-chain split.
  const exprSystems = ['CHO', 'HEK293', 'E.coli'];
  const exprValues: string[] = [];
  new Array(sequenceIds.length * 3).fill(null).map(() => randPick(sequenceIds)).forEach((id) => {
    const ch = randPick(chains as unknown as string[]);
    exprValues.push(`(${id}, '${ch}', '${escape(randPick(exprSystems))}', ${ (Math.random() * 150).toFixed(2) }, 'ExprBatch_${id}_${ch}_${(Math.floor(Math.random() * 1000))}', 'Auto-generated expression batch (${ch})')`);
  });
  if (exprValues.length) {
    await exec(`
      INSERT INTO biologics.expression_batches(sequence_id, chain, expression_system, yield_mg, name, notes)
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

  // Helper: bulk-insert assay results (with optional peptide_id), return newly-inserted ids
  // in the same order, and queue curves for rows whose assay type has a fit_model.
  const curveQueue: { resultId: number, curve: string }[] = [];

  async function insertAssayResultsBatch(
    rows: { assayId: number, refId: number, refKind: 'adc' | 'peptide', org: number, value: number, units: string, target: number }[]
  ): Promise<number[]> {
    const allIds: number[] = [];
    if (!rows.length) return allIds;
    const chunkSize = 200;
    for (let i = 0; i < rows.length; i += chunkSize) {
      const part = rows.slice(i, i + chunkSize);
      const adcPart = part.filter((r) => r.refKind === 'adc');
      const pepPart = part.filter((r) => r.refKind === 'peptide');
      // Group by ref kind because the column list differs.
      for (const grp of [
        {kind: 'adc' as const, list: adcPart, col: 'adc_id'},
        {kind: 'peptide' as const, list: pepPart, col: 'peptide_id'},
      ]) {
        if (!grp.list.length) continue;
        const values = grp.list
          .map((r) => `(${r.assayId}, ${r.refId}, ${r.org}, ${r.value}, '${escape(r.units)}', ${r.target})`)
          .join(',');
        const df = await exec(`
          INSERT INTO biologics.assay_results(assay_id, ${grp.col}, target_organism_id, result_value, units, target_id)
          VALUES ${values}
          RETURNING id, assay_id
        `);
        for (let r = 0; r < df.rowCount; r++) {
          const id = df.get('id', r);
          const assayId = df.get('assay_id', r);
          allIds.push(id);
          if (curveEligibleAssayIds.has(assayId))
            curveQueue.push({resultId: id, curve: randPick(curvePool)});
        }
      }
    }
    return allIds;
  }

  // 9. Per-sequence (ADC) assay results, mostly legacy types
  const legacyAssayRows: Parameters<typeof insertAssayResultsBatch>[0] = [];
  const maxResults = 5000; // cap
  let inserted = 0;
  outer: for (const _seqId of sequenceIds) {
    const perSeqAssays = randInt(3, 26);
    for (let i = 0; i < perSeqAssays; i++) {
      const at = randPick(legacyTypes.length ? legacyTypes : assayTypes);
      const av = randomAssayValue(at.name);
      legacyAssayRows.push({
        assayId: at.id, refId: randPick(adcIds), refKind: 'adc',
        org: randPick(organisms), value: av.val, units: av.units, target: randPick(targets),
      });
      inserted++;
      if (inserted >= maxResults) break outer;
    }
  }
  await insertAssayResultsBatch(legacyAssayRows);

  // 10. Per-peptide assay results
  const peptideAssayRows: Parameters<typeof insertAssayResultsBatch>[0] = [];
  const maxPeptideResults = 2000;
  let peptideAssayInserted = 0;
  outer2: for (const peptide of insertedHelms) {
    const perPeptideAssays = randInt(3, 10);
    for (let i = 0; i < perPeptideAssays; i++) {
      const at = randPick(assayTypes);
      const av = randomAssayValue(at.name);
      peptideAssayRows.push({
        assayId: at.id, refId: peptide.id, refKind: 'peptide',
        org: randPick(organisms), value: av.val, units: av.units, target: randPick(targets),
      });
      peptideAssayInserted++;
      if (peptideAssayInserted >= maxPeptideResults) break outer2;
    }
  }
  await insertAssayResultsBatch(peptideAssayRows);

  // 11. Comprehensive assay panel: every ADC gets one result per new-panel assay type
  let comprehensiveInserted = 0;
  if (newPanelTypes.length && adcIds.length) {
    const compRows: Parameters<typeof insertAssayResultsBatch>[0] = [];
    for (const adcId of adcIds) {
      const tgt = randPick(targets);
      const org = randPick(organisms);
      for (const at of newPanelTypes) {
        const av = randomAssayValue(at.name);
        compRows.push({
          assayId: at.id, refId: adcId, refKind: 'adc',
          org, value: av.val, units: av.units, target: tgt,
        });
      }
    }
    comprehensiveInserted = compRows.length;
    await insertAssayResultsBatch(compRows);
  }

  // 11b. Persist queued curves
  let curvesInserted = 0;
  if (curveQueue.length) {
    const chunkSize = 50; // curves are large JSON blobs; keep statements small
    for (let i = 0; i < curveQueue.length; i += chunkSize) {
      const part = curveQueue.slice(i, i + chunkSize);
      const values = part.map((c) => `(${c.resultId}, '${escape(c.curve)}')`).join(',');
      await exec(`INSERT INTO biologics.assay_curves(assay_result_id, curve) VALUES ${values}`);
      curvesInserted += part.length;
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
    assay_curves: curvesInserted,
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
  const curvesGcdFunc = DG.Func.find({name: '_FitChartCellRenderer'})[0];
  const curveGCRenderer = curvesGcdFunc ? (await curvesGcdFunc.apply({})) as DG.GridCellRenderer : null;
  const curvesRenderer = (curve: string) => {
    const host = ui.div();
    if (!curveGCRenderer || !curve) return host;
    try {
      const gc = DG.GridCell.fromValue(curve);
      const canvas = ui.canvas(300, 200);
      curveGCRenderer.render(canvas.getContext('2d')!, 0, 0, 300, 200, gc, null);
      host.appendChild(canvas);
    } catch (e) {
      console.error('Error rendering curve:', e);
      host.textContent = 'Error rendering curve';
    }
    return host;
  };


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

  exp.addCustomRenderer((_, colName, value) => colName?.toLowerCase()?.includes('curve') && typeof value === 'string',
    (value) => curvesRenderer(value as string)
  );
  // exp.addDefaultHeaderReplacerColumns(['units']);
  console.log('Biologics object handlers registered');
}
