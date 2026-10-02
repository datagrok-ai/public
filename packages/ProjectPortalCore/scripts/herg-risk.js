//name: HergRisk
//friendlyName: hERG risk (example)
//description: Synthetic hERG IC50 estimate for the example small-molecule portfolio (EX-SMD): a made-up function of potency, so a compound without a hERG measurement has a stand-in. Invented for the example; not a safety model and not fitted to anything.
//language: javascript
//meta.role: model
//meta.domain: ADMET
//input: dataframe data {caption: SAR table; description: One row per compound, with compound and ic50_um columns}
//input: double threshold = 10 {caption: Flag below; units: uM; min: 1; max: 100}
//output: dataframe predictions {caption: Estimated hERG IC50}

// Made up for the example: round numbers chosen to give plausible values, nothing more.
const n = data.rowCount;
const compound = data.columns.byName('compound');
const ic50 = data.columns.byName('ic50_um');
const measured = data.columns.byName('herg_ic50_um');
const ids = [];
const flags = [];
const est = new Float32Array(n);
const meas = new Float32Array(n);
for (let i = 0; i < n; i++) {
  ids.push(compound ? String(compound.get(i)) : `row ${i + 1}`);
  const x = ic50 && !ic50.isNone(i) ? ic50.get(i) : 1;
  est[i] = Math.round(10 * 30 * Math.pow(Math.max(x, 0.001), 0.25)) / 10;
  meas[i] = measured && !measured.isNone(i) ? measured.get(i) : DG.FLOAT_NULL;
  flags.push(est[i] < threshold ? 'risk' : '');
}
predictions = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('compound', ids),
  DG.Column.fromFloat32Array('herg_ic50_um', meas),
  DG.Column.fromFloat32Array('estimated_herg_ic50_um', est),
  DG.Column.fromStrings('flag', flags),
]);
predictions.name = 'Estimated hERG IC50 (example)';
