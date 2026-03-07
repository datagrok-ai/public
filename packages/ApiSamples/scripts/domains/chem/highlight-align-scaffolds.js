// Highlights or align scaffolds

const t = grok.data.demo.molecules();

// Highlight multiple scaffolds
const highlightOptions = [
  {"color":"#00FF00","molecule":"\nActelion Java MolfileCreator 1.0\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n   11.1250   -8.3125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   11.1250   -9.8125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   12.4240  -10.5625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   13.7231   -9.8125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   13.7231   -8.3125   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   12.4240   -7.5625   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  1  1  0  0  0  0\nM  END\n"},
  {"color":"#FF0000","molecule": "CO"}
];
t.col('smiles').setTag('.%chem-scaffold-highlight', JSON.stringify(highlightOptions));

// Align to a scaffold
const alignedCol = t.columns.addNew('aligned', DG.TYPE.STRING).init((i) => t.col('smiles').get(i));
alignedCol.setTag('.%chem-scaffold-align', '\nActelion Java MolfileCreator 1.0\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n   13.1250  -12.3750   -0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   14.4240  -11.6250   -0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\nM  END\n');
grok.shell.addTableView(t);