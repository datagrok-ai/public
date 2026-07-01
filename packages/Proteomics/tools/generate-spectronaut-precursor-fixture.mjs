#!/usr/bin/env node
/*
 * Generates files/demo/spectronaut-hye-precursor.tsv — a small SYNTHETIC
 * precursor-level Spectronaut long-format report used to exercise the streaming
 * import path (Plan 02) and the duckdb-equivalence golden test (Plan 03).
 *
 * The fixture carries the D-01 precursor signature columns
 * (EG.ModifiedPeptide / FG.Charge / PEP.StrippedSequence) so the header-sniff
 * routes it to the streaming aggregator, uses CondA/CondB (so the committed
 * tools/spectronaut-aggregate.sql DMD↔WT flip is a structural no-op), keeps the
 * quantity value CONSTANT per (protein × condition × replicate) so duckdb max()
 * equals the parser's first-encountered value, and carries NO
 * PG.Genes/PG.ProteinAccessions columns (the real hye-mix demo lacks them and
 * the committed SQL drops the matching any_value SELECT terms for this reason).
 *
 * It deliberately includes one CON__ row, one REV__ row, one protein whose
 * precursors mix a >0.01 sub-threshold q-value with a passing one, one
 * non-numeric-q-value protein, and one empty-q-value protein so every R2 filter
 * branch is exercised by a single fixture.
 *
 * Pure standalone Node (no datagrok-api import); the row shape mirrors the
 * extended makeLongFormatTsv in src/tests/spectronaut-parser.ts.
 *
 * Run:  node tools/generate-spectronaut-precursor-fixture.mjs
 */

import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const PKG_ROOT = path.resolve(__dirname, '..');
const OUTPUT = path.join(PKG_ROOT, 'files/demo/spectronaut-hye-precursor.tsv');

const CONDITIONS = ['CondA', 'CondB'];
const REPLICATES = [1, 2, 3];

// Column order mirrors the extended makeLongFormatTsv header. NO
// PG.Genes/PG.ProteinAccessions — see doc block.
const HEADERS = [
  'R.FileName', 'R.Condition', 'R.Replicate', 'PG.ProteinGroups',
  'PG.Organisms', 'PG.Quantity', 'EG.Qvalue', 'PEP.StrippedSequence',
  'EG.ModifiedPeptide', 'FG.Charge', 'FG.Id',
];

// Two distinct precursors per (protein × condition × replicate). Quantity is
// held constant within each group; only precursor identity varies, so duckdb
// max() collapses to the same value the parser's first-encountered-wins logic
// produces.
const PRECURSORS = [
  {strip: 'PEPTIDER', mod: '_PEPTIDER_', charge: '2', fg: 'FG1'},
  {strip: 'ANOTHERPEP', mod: '_ANOTHERPEP_', charge: '3', fg: 'FG2'},
];

const ORGANISMS = ['Homo sapiens', 'Saccharomyces cerevisiae', 'Escherichia coli'];

const rows = [HEADERS.join('\t')];
let nDataRows = 0;
const proteinSet = new Set();
const sampleSet = new Set();

function emit(fileName, cond, rep, id, organism, quantity, qVal) {
  for (const p of PRECURSORS) {
    rows.push([fileName, cond, String(rep), id, organism,
      String(quantity), String(qVal), p.strip, p.mod, p.charge, p.fg].join('\t'));
    nDataRows++;
  }
}

function emitProtein(id, organism, qVal, quantity) {
  proteinSet.add(id);
  for (const cond of CONDITIONS) {
    for (const rep of REPLICATES) {
      const fileName = `run_${cond}_${rep}`;
      sampleSet.add(`${cond}_${rep}`);
      emit(fileName, cond, rep, id, organism, quantity, qVal);
    }
  }
}

// ── Regular passing proteins: a few dozen, three-species mix, all q<=0.01 ──
const N_REGULAR = 36;
let quant = 1000;
for (let i = 0; i < N_REGULAR; i++) {
  const id = `P${String(i).padStart(5, '0')}`;
  const organism = ORGANISMS[i % ORGANISMS.length];
  emitProtein(id, organism, 0.001, quant++);
}

// ── Decoy / contaminant rows (dropped by the WHERE CON__/REV__ filter) ──
emitProtein('CON__P99999', 'Homo sapiens', 0.001, 9000);
emitProtein('REV__Q88888', 'Homo sapiens', 0.001, 9001);

// ── Mixed sub-threshold protein: one passing precursor, one >0.01 dropped ──
// Emit the rows by hand so the two precursors carry DIFFERENT q-values.
{
  const id = 'P70001';
  proteinSet.add(id);
  for (const cond of CONDITIONS) {
    for (const rep of REPLICATES) {
      const fileName = `run_${cond}_${rep}`;
      sampleSet.add(`${cond}_${rep}`);
      // passing precursor
      rows.push([fileName, cond, String(rep), id, 'Homo sapiens',
        '7000', '0.002', PRECURSORS[0].strip, PRECURSORS[0].mod,
        PRECURSORS[0].charge, PRECURSORS[0].fg].join('\t'));
      nDataRows++;
      // sub-threshold precursor (q-value 0.05 > 0.01 → dropped by WHERE)
      rows.push([fileName, cond, String(rep), id, 'Homo sapiens',
        '7000', '0.05', PRECURSORS[1].strip, PRECURSORS[1].mod,
        PRECURSORS[1].charge, PRECURSORS[1].fg].join('\t'));
      nDataRows++;
    }
  }
}

// ── Non-numeric q-value protein ('Profiled' → TRY_CAST NULL → passes) ──
emitProtein('P70002', 'Saccharomyces cerevisiae', 'Profiled', 6500);

// ── Empty-string q-value protein (nullstr → NULL → passes) ──
emitProtein('P70003', 'Escherichia coli', '', 6600);

fs.writeFileSync(OUTPUT, rows.join('\n') + '\n');
console.log(`Wrote ${OUTPUT}`);
console.log(`  ${nDataRows} data rows, ${proteinSet.size} protein groups ` +
  `(incl. CON__/REV__/edge), ${sampleSet.size} samples ` +
  `(${CONDITIONS.length} conditions × ${REPLICATES.length} replicates).`);
