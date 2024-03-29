// Splitters

let t = grok.data.demo.demog();
let block = ui.splitV([
  ui.splitH([
    ui.inputs([
      ui.stringInput('Input', ''),
      ui.input.choice('Select', {items: ['item 1', 'item 2'], value: 'item 1'}),
      ui.textInput('text', ' multi line /n text input'),
    ]),
    DG.Viewer.grid(t).root
  ], {}, true),
  ui.splitH([
    DG.Viewer.scatterPlot(t).root,
    DG.Viewer.histogram(t).root
  ], {}, true)
], {}, true);

let view = grok.shell.newView('Splitters');
view.box = true;
view.append(block);