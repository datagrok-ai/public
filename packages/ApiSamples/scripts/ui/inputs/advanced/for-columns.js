// Create inputs dynamically for columns

let t = grok.data.demo.molecules();
ui.divV(t.columns.toList().map(c => DG.InputBase.forColumn(c)));