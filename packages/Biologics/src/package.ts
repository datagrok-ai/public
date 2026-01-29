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
import {moleculeRenderer, imageRenderer, rawImageRenderer} from '@datagrok-libraries/db-explorer/src/renderer';
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

  // 4. Fetch existing assay_types & target organisms
  const assayTypesDf = await exec('SELECT id, name FROM biologics.assay_types');
  const assayTypes: { id: number, name: string }[] = [];
  for (let r = 0; r < assayTypesDf.rowCount; r++)
    assayTypes.push({id: assayTypesDf.get('id', r), name: assayTypesDf.get('name', r)});

  const orgDf = await exec('SELECT id, name FROM biologics.target_organisms');
  const organisms: number[] = [];
  for (let r = 0; r < orgDf.rowCount; r++)
    organisms.push(orgDf.get('id', r));

  // 5. Insert everything
  const sequenceIds = await insertSequences(sequences);
  const drugIds = await insertDrugs(drugSmiles);
  const linkerIds = await insertLinkers(linkers);

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
    if (/IC50/i.test(name)) return {val: +(Math.random() * 500).toFixed(2), units: 'nM'};
    if (/Caspase/i.test(name)) return {val: +(Math.random() * 10).toFixed(2), units: 'RFU'};
    if (/half-life/i.test(name)) return {val: +(Math.random() * 240).toFixed(1), units: 'min'};
    if (/Binding affinity/i.test(name)) return {val: +(Math.random() * 50).toFixed(2), units: 'nM'};
    if (/Cell binding/i.test(name)) return {val: +(Math.random() * 100).toFixed(2), units: 'nM'};
    if (/DAR/i.test(name)) return {val: +(2 + Math.random() * 6).toFixed(2), units: 'ratio'};
    if (/Cmax/i.test(name)) return {val: +(Math.random() * 100).toFixed(2), units: 'µg/mL'};
    if (/Tmax/i.test(name)) return {val: +(Math.random() * 48).toFixed(2), units: 'h'};
    if (/AUC/i.test(name)) return {val: +(Math.random() * 5000).toFixed(1), units: 'µg·h/mL'};
    return {val: +(Math.random() * 100).toFixed(2), units: ''};
  }

  const assayResultChunks: string[] = [];
  const maxResults = 5000; // cap
  let inserted = 0;
  for (const seqId of sequenceIds) {
    const perSeqAssays = randInt(3, 26);
    for (let i = 0; i < perSeqAssays; i++) {
      const at = randPick(assayTypes);
      const org = randPick(organisms);
      const av = randomAssayValue(at.name);
      assayResultChunks.push(`(${at.id}, ${randPick(adcIds)}, ${org}, ${av.val}, '${escape(av.units)}')`);
      inserted++;
      if (inserted >= maxResults) break;
    }
    if (inserted >= maxResults) break;
  }
  if (assayResultChunks.length) {
    // Insert in chunks to avoid overly large statements
    const chunkSize = 200;
    for (let i = 0; i < assayResultChunks.length; i += chunkSize) {
      const part = assayResultChunks.slice(i, i + chunkSize);
      await exec(`
        INSERT INTO biologics.assay_results(assay_id, adc_id, target_organism_id, result_value, units)
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
      const av = randomAssayValue(at.name);
      peptideAssayChunks.push(`(${at.id}, ${peptide.id}, ${org}, ${av.val}, '${escape(av.units)}')`);
      peptideAssayInserted++;
      if (peptideAssayInserted >= maxPeptideResults) break;
    }
    if (peptideAssayInserted >= maxPeptideResults) break;
  }
  if (peptideAssayChunks.length) {
    // Insert in chunks to avoid overly large statements
    const chunkSize = 200;
    for (let i = 0; i < peptideAssayChunks.length; i += chunkSize) {
      const part = peptideAssayChunks.slice(i, i + chunkSize);
      await exec(`
        INSERT INTO biologics.assay_results(assay_id, peptide_id, target_organism_id, result_value, units)
        VALUES ${part.join(',')}
      `);
    }
  }


  await populateAdcGlyphs(adcCount + 1);

  return {
    sequences: sequenceIds.length,
    drugs: drugIds.length,
    linkers: linkerIds.length,
    assay_results: inserted,
    purification_batches: purBatchesValues.length,
    expression_batches: exprValues.length,
    adcs: adcCount,
    peptides: insertedHelms.length,
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
  function renderHelm(value: string): HTMLElement {
    return ui.wait(async () => {
    //@ts-ignore
      const helmInput = await ui.input.helmAsync('helm', {
        editable: false,
      });
      helmInput.setStringValue(value);
      await DG.delay(200); // wait for proper sizing
      helmInput.getInput().addEventListener('click', () => {
        grok.shell.o = helmInput.getValue();
      });
      helmInput.getInput().addEventListener('dblclick', () => {
        helmInput.showEditorDialog();
      });

      helmInput.getInput().style.width = '400px';
      helmInput.getInput().style.setProperty('height', '300px', 'important');
      return helmInput.getInput();
    });
  }

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
    (value) => renderHelm(value as string)
  );

  exp.addCustomRenderer((_, colName, value) => {
    const lc = colName?.toLowerCase() || '';
    return (lc === 'glyph' || lc === 'image' || lc === 'png' || lc === 'thumbnail') && typeof value === 'string' && value.startsWith('iVBORw0KGgo');
  }, (value) => rawImageRenderer(value as string));

  // exp.addDefaultHeaderReplacerColumns(['units']);
  console.log('Biologics object handlers registered');
}
