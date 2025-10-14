// The default Form viewer shows 20 columns ranked by relevance.
// The list of columns can be modified to match other heuristics.
const df = await grok.data.getDemoTable('chem/smiles.csv');
const view = grok.shell.addTableView(df);

grok.events.onFormCreating.subscribe((data) => { 
  data.columns = data.filter((col) => col.semType === DG.SEMTYPE.MOLECULE);
});
