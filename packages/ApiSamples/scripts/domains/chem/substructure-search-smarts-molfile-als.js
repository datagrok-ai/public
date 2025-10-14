// An example of using substructure search with SMARTS pattern,
// encoded through an ALS parameter block of a MolFile

let t = await grok.data.getDemoTable("chem/smiles.csv");
let bs = await grok.chem.searchSubstructure(
  t.col("canonical_smiles"),
  `
  MJ201900\n
  6  6  1  0  0  0  0  0  0  0999 V2000
   -1.9420    0.9589    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564    0.5464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6564   -0.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9420   -0.6911    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275   -0.2786    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2275    0.5464    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  2  2  0  0  0  0
  6  1  1  0  0  0  0
  1 T    2   6   7
M  ALS   1  2 T C   N   \n 
M  END
`
);

t.filter.copyFrom(bs);
grok.shell.addTableView(t);
