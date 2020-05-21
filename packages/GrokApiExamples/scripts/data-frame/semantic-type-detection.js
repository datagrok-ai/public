// Custom handling of the semantic type detection
// See https://datagrok.ai/help/discover/semantic-types

let demog = grok.parseCsv('smiles\nFc1cc(Cl)ccc1Br');

// Normally, the smiles column will be detected as 'Molecule' when we create a table view.
// However, this time we will override it

demog.onSemanticTypeDetecting.subscribe((args) => {
    demog.col('smiles').semType = 'bananas';
    args.preventDefault();
});
demog.onSemanticTypeDetected.subscribe((_) => { grok.info(demog.col('smiles').semType); });

grok.addTableView(demog);