// HTML cells in a dataframe

let t = DG.DataFrame.fromColumns([
  DG.Column.fromStrings('test', [
    '<button>Click me</button>',
    '<div style="background: red">Div</div>',
    '<input type="range">'
  ])
]);
t.col('test').tags['cell.renderer'] = 'html';

grok.shell.addTableView(t);
