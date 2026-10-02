//name: SolubilityEstimate
//friendlyName: Solubility estimate (example)
//description: Synthetic solubility estimate for the example small-molecule portfolio (EX-SMD): a made-up function of potency, so a gap in the SAR table has a stand-in. Invented for the example; not a QSPR model and not fitted to anything.
//language: javascript
//meta.role: model
//meta.domain: Physchem
//input: dataframe data {caption: SAR table; description: One row per compound, with compound and ic50_um columns}
//output: dataframe predictions {caption: Estimated solubility}

// Made up for the example: round numbers chosen to give plausible values, nothing more.
const n = data.rowCount;
const compound = data.columns.byName('compound');
const ic50 = data.columns.byName('ic50_um');
const measured = data.columns.byName('solubility_um');
const ids = [];
const est = new Float32Array(n);
const meas = new Float32Array(n);
for (let i = 0; i < n; i++) {
  ids.push(compound ? String(compound.get(i)) : `row ${i + 1}`);
  const x = ic50 && !ic50.isNone(i) ? ic50.get(i) : 1;
  est[i] = Math.round(10 * 90 / (1 + 3 * Math.sqrt(Math.max(x, 0)))) / 10;
  meas[i] = measured && !measured.isNone(i) ? measured.get(i) : DG.FLOAT_NULL;
}
predictions = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('compound', ids),
  DG.Column.fromFloat32Array('solubility_um', meas),
  DG.Column.fromFloat32Array('estimated_solubility_um', est),
]);
predictions.name = 'Estimated solubility (example)';
