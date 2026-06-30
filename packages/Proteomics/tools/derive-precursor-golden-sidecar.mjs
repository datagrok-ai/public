#!/usr/bin/env node
/*
 * Derives files/demo/spectronaut-hye-precursor-golden.json from the committed
 * duckdb golden files/demo/spectronaut-hye-precursor-golden.tsv.
 *
 * This script does NO aggregation — it is a pure, lossless transcription of the
 * already-committed duckdb output into a JSON map the in-browser `grok test`
 * runner can read with certainty when committed-file reads are unavailable in
 * the test runner. Because it never re-computes anything, the equivalence test
 * stays pinned to REAL duckdb output (the D-04 oracle), not to a hand-derived
 * approximation. Re-running it against an unchanged golden produces a
 * byte-identical JSON.
 *
 * Regen chain (run in order, or the artifacts silently drift):
 *   1. node tools/generate-spectronaut-precursor-fixture.mjs
 *   2. tools/spectronaut-aggregate.sh files/demo/spectronaut-hye-precursor.tsv \
 *        files/demo/spectronaut-hye-precursor-golden.tsv
 *   3. node tools/derive-precursor-golden-sidecar.mjs   (this script)
 *
 * Run:  node tools/derive-precursor-golden-sidecar.mjs
 */

import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const PKG_ROOT = path.resolve(__dirname, '..');
const GOLDEN = path.join(PKG_ROOT, 'files/demo/spectronaut-hye-precursor-golden.tsv');
const OUTPUT = path.join(PKG_ROOT, 'files/demo/spectronaut-hye-precursor-golden.json');

// Quantity column preference order — same as the in-package Spectronaut parser
// (src/parsers/spectronaut-parser.ts QUANTITY_COLUMNS).
const QUANTITY_COLUMNS = ['PG.IBAQ', 'PG.Quantity'];

const text = fs.readFileSync(GOLDEN, 'utf8');
const lines = text.split(/\r?\n/).filter((l) => l.length > 0);
if (lines.length < 2)
  throw new Error(`Golden ${GOLDEN} has no data rows — run the aggregate wrapper first.`);

const headers = lines[0].split('\t');
const idx = Object.fromEntries(headers.map((h, i) => [h, i]));

const protIdx = idx['PG.ProteinGroups'];
const condIdx = idx['R.Condition'];
const replIdx = idx['R.Replicate'];
const qvalIdx = idx['EG.Qvalue'];
const quantCol = QUANTITY_COLUMNS.find((c) => idx[c] !== undefined);

if (protIdx === undefined || condIdx === undefined || replIdx === undefined ||
    qvalIdx === undefined || quantCol === undefined) {
  throw new Error(
    `Golden ${GOLDEN} is missing an expected column. ` +
    `Need PG.ProteinGroups / R.Condition / R.Replicate / EG.Qvalue / ` +
    `${QUANTITY_COLUMNS.join('|')}. Got: ${headers.join(', ')}`);
}
const quantIdx = idx[quantCol];

const out = {};
for (let i = 1; i < lines.length; i++) {
  const f = lines[i].split('\t');
  const protein = f[protIdx];
  const condition = f[condIdx];
  const replicate = f[replIdx];
  // Numbers exactly as duckdb wrote them. Empty q-value (TRY_CAST→NULL min,
  // e.g. the non-numeric / empty-string-q-value proteins) → Number('') is 0;
  // keep the verbatim transcription rather than inventing a sentinel.
  const key = `${protein}${condition}_${replicate}`;
  out[key] = {
    quantity: Number(f[quantIdx]),
    qvalue: Number(f[qvalIdx] ?? ''),
  };
}

fs.writeFileSync(OUTPUT, JSON.stringify(out, null, 2) + '\n');
const keys = Object.keys(out);
console.log(`Wrote ${OUTPUT}`);
console.log(`  ${keys.length} (protein × condition × replicate) entries ` +
  `transcribed verbatim from ${path.basename(GOLDEN)} ` +
  `(quantity column: ${quantCol}).`);
