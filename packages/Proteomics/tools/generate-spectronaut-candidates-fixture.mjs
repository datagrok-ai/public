#!/usr/bin/env node
/*
 * Generates files/demo/spectronaut-hye-candidates.tsv from the already-shipped
 * files/demo/spectronaut-hye-mix.tsv (HYE three-species mix, SpectroPipeR origin).
 *
 * Source data is real Spectronaut long-format PG quantification. The Candidates
 * report it ships alongside is DERIVED — we run our own Welch's t-test with BH
 * FDR correction in this script and emit the result formatted as a Spectronaut
 * Candidates report. This is not a verbatim Spectronaut export; it exists to
 * exercise the Candidates import path with values that line up with the source
 * PG file so the two demo fixtures tell the same story.
 *
 * Run:  node tools/generate-spectronaut-candidates-fixture.mjs
 */

import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
import jStat from 'jstat';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const PKG_ROOT = path.resolve(__dirname, '..');
const SOURCE = path.join(PKG_ROOT, 'files/demo/spectronaut-hye-mix.tsv');
const OUTPUT = path.join(PKG_ROOT, 'files/demo/spectronaut-hye-candidates.tsv');

// Numerator (group1) and denominator (group2) per Spectronaut convention.
// AVG Log2 Ratio = log2(numerator / denominator). The HYE mix is symmetric;
// labelling B as the numerator just keeps the convention legible.
const NUMERATOR = 'HYE mix B';
const DENOMINATOR = 'HYE mix A';

/** Read TSV file into header + rows. */
function readTsv(filePath) {
  const text = fs.readFileSync(filePath, 'utf8');
  const lines = text.split(/\r?\n/);
  const headers = lines[0].split('\t');
  const idx = Object.fromEntries(headers.map((h, i) => [h, i]));
  const rows = [];
  for (let i = 1; i < lines.length; i++) {
    if (!lines[i]) continue;
    const fields = lines[i].split('\t');
    rows.push(fields);
  }
  return {headers, idx, rows};
}

/** Pivot long PG report: protein -> {sampleKey -> IBAQ}. First-encountered
 * non-empty IBAQ wins per protein × sample (constant per protein × sample in
 * SpectroPipeR's report; mirrors the in-package Spectronaut parser). */
function pivotProteinIbaq(headers, idx, rows) {
  const protein = {};
  for (const r of rows) {
    const pg = r[idx['PG.ProteinGroups']];
    const cond = r[idx['R.Condition']];
    const repl = r[idx['R.Replicate']];
    const ibaqRaw = r[idx['PG.IBAQ']];
    const org = r[idx['PG.Organisms']];
    if (!pg || !cond || !repl) continue;
    const ibaq = Number(ibaqRaw);
    if (!Number.isFinite(ibaq) || ibaq <= 0) continue;
    const sampleKey = `${cond}|${repl}`;
    if (!protein[pg]) protein[pg] = {samples: {}, organism: org || ''};
    if (protein[pg].samples[sampleKey] === undefined)
      protein[pg].samples[sampleKey] = ibaq;
    if (!protein[pg].organism && org) protein[pg].organism = org;
  }
  return protein;
}

/** Welch's two-sample t-test using jStat. Returns null when either sample has
 * fewer than two finite values. */
function welchTTest(a, b) {
  const af = a.filter(Number.isFinite);
  const bf = b.filter(Number.isFinite);
  if (af.length < 2 || bf.length < 2) return null;
  const meanA = jStat.mean(af);
  const meanB = jStat.mean(bf);
  const varA = jStat.variance(af, true);
  const varB = jStat.variance(bf, true);
  const seA = varA / af.length;
  const seB = varB / bf.length;
  const se = Math.sqrt(seA + seB);
  if (se === 0) return {t: 0, df: af.length + bf.length - 2, p: 1, meanA, meanB};
  const t = (meanB - meanA) / se;
  // Welch–Satterthwaite degrees of freedom
  const df = Math.pow(seA + seB, 2) /
    (Math.pow(seA, 2) / (af.length - 1) + Math.pow(seB, 2) / (bf.length - 1));
  const pTwoSided = 2 * (1 - jStat.studentt.cdf(Math.abs(t), df));
  return {t, df, p: pTwoSided, meanA, meanB};
}

/** Benjamini–Hochberg FDR correction. Returns array of adjusted p-values aligned
 * to input order; entries left as null where the input is null. */
function bhAdjust(pvals) {
  const indexed = pvals.map((p, i) => ({p, i})).filter((x) => x.p !== null);
  indexed.sort((a, b) => a.p - b.p);
  const n = indexed.length;
  const adj = new Array(pvals.length).fill(null);
  let cumMin = 1;
  for (let rank = n; rank >= 1; rank--) {
    const e = indexed[rank - 1];
    const q = e.p * n / rank;
    cumMin = Math.min(cumMin, q);
    adj[e.i] = Math.min(cumMin, 1);
  }
  return adj;
}

function main() {
  const {headers, idx, rows} = readTsv(SOURCE);
  console.log(`Read ${rows.length} long-format rows, ${headers.length} columns.`);

  const protein = pivotProteinIbaq(headers, idx, rows);
  const proteins = Object.keys(protein).sort();
  console.log(`Pivoted to ${proteins.length} protein groups.`);

  // Sample keys per group
  const sampleKey = (c, r) => `${c}|${r}`;
  const groupSamples = (cond) => [1, 2, 3, 4].map((r) => sampleKey(cond, String(r)));
  const numKeys = groupSamples(NUMERATOR);
  const denKeys = groupSamples(DENOMINATOR);

  // Per-protein log2FC + raw p-value
  const proteinStats = proteins.map((pg) => {
    const samples = protein[pg].samples;
    const num = numKeys.map((k) => samples[k]).filter(Number.isFinite);
    const den = denKeys.map((k) => samples[k]).filter(Number.isFinite);
    const numLog2 = num.map((x) => Math.log2(x));
    const denLog2 = den.map((x) => Math.log2(x));
    const t = welchTTest(denLog2, numLog2);  // a=denominator, b=numerator
    if (!t) {
      return {
        pg,
        organism: protein[pg].organism,
        log2fc: null,
        p: null,
        meanNum: num.length > 0 ? jStat.mean(num) : null,
        meanDen: den.length > 0 ? jStat.mean(den) : null,
        nRatios: Math.min(num.length, den.length),
        valid: false,
      };
    }
    return {
      pg,
      organism: protein[pg].organism,
      log2fc: t.meanB - t.meanA,  // log2(num) - log2(den)
      p: t.p,
      meanNum: jStat.mean(num),
      meanDen: jStat.mean(den),
      nRatios: Math.min(num.length, den.length),
      valid: true,
    };
  });

  const adj = bhAdjust(proteinStats.map((s) => s.p));

  // Write Candidates TSV
  const outHeaders = [
    'Valid', 'Comparison (group1/group2)',
    'Condition Numerator', 'Condition Denominator',
    '# of Ratios', 'Group',
    'AVG Group Quantity Numerator', 'AVG Group Quantity Denominator',
    'AVG Log2 Ratio', 'Absolute AVG Log2 Ratio',
    '% Change', 'Ratio',
    'Pvalue', 'Qvalue',
    'ProteinGroups', 'Genes', 'UniProtIds', 'Organisms',
  ];
  const comparison = `${NUMERATOR} / ${DENOMINATOR}`;

  const lines = [outHeaders.join('\t')];
  let nSig = 0, nUp = 0, nDown = 0;
  for (let i = 0; i < proteinStats.length; i++) {
    const s = proteinStats[i];
    const q = adj[i];
    const ratio = (s.meanNum != null && s.meanDen != null && s.meanDen > 0)
      ? (s.meanNum / s.meanDen) : null;
    const pctChange = ratio != null ? (ratio - 1) * 100 : null;
    const row = [
      s.valid ? 'True' : 'False',
      comparison,
      NUMERATOR,
      DENOMINATOR,
      String(s.nRatios),
      s.pg,
      s.meanNum != null ? s.meanNum.toFixed(4) : '',
      s.meanDen != null ? s.meanDen.toFixed(4) : '',
      s.log2fc != null ? s.log2fc.toFixed(6) : '',
      s.log2fc != null ? Math.abs(s.log2fc).toFixed(6) : '',
      pctChange != null ? pctChange.toFixed(4) : '',
      ratio != null ? ratio.toFixed(6) : '',
      s.p != null ? s.p.toExponential(6) : '',
      q != null ? q.toExponential(6) : '',
      s.pg,                                        // ProteinGroups (same as Group here)
      '',                                          // Genes (HYE PG report has no symbols)
      s.pg,                                        // UniProtIds (PG.ProteinGroups in this dataset are UniProt accs)
      s.organism,
    ];
    lines.push(row.join('\t'));

    if (s.log2fc != null && q != null && Math.abs(s.log2fc) >= 1 && q <= 0.05) {
      nSig++;
      if (s.log2fc > 0) nUp++; else nDown++;
    }
  }

  fs.writeFileSync(OUTPUT, lines.join('\n') + '\n');
  console.log(`Wrote ${OUTPUT}`);
  console.log(`  ${proteinStats.length} rows; ${nSig} significant (${nUp} up, ${nDown} down) at |log2FC|>=1, q<=0.05`);
}

main();
