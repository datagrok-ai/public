<!-- Full file is too large to inline; the targeted edit was applied directly.
     The change adds a "Programmatic usage (JS API)" <details> block
     immediately after the existing "How to use" block in the
     "### Elemental analysis" section, before "## Structure relationship analysis".

     Inserted content:

<details>
<summary>Programmatic usage (JS API)</summary>

Call `Chem:elementalAnalysis` via `grok.functions.call`. Pass the source
DataFrame, the molecular column, and flags to control which visualization(s)
to produce:

```javascript
await grok.functions.call('Chem:elementalAnalysis', {
  table: t,
  molecules: t.col('smiles'), // DG.Column containing SMILES / molblock
  radarViewer: true,          // add a standalone Radar viewer
  radarGrid: false            // embed sparklines in grid cells
});
```

After the call completes, atom-count columns (one per element present, in
periodic-table order) and a **Molecule Charge** column are appended to `t`,
and the chosen radar visualization is displayed.

</details>
-->