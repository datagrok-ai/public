// Bioreactor simulation wireframe

let view = grok.shell.newView("Flowrate");
let progressDiv = ui.divText('Ready');
view.setRibbonPanels([[progressDiv]]);

let data = grok.data.demo.randomWalk();

let reactorInput = ui.choiceInput('Reactor', 'Biotek 15kL SS', ['Biotek 15kL SS']);
let productInput = ui.choiceInput('Product', 'BAR FOO56545', ['BAR FOO56545']);
let controlRunInput = ui.choiceInput('Control Run', '15FDX006', ['15FDX006']);
let daysInput = ui.intInput('Process days', 12);

let topGrid = ui.divV([
  ui.divH([
    //ui.divText('FEED BY TIME'),
    ui.divH([ui.button('Use historical feed')])
  ]),
  DG.Viewer.grid(data).root,
]);
topGrid.style.flexGrow = '1';
topGrid.style.marginLeft = '20px';
topGrid.style.minHeight = '300px';

let top = ui.divH([
  ui.inputs([reactorInput, productInput, controlRunInput, daysInput]),
  topGrid
]);

let resultsDiv = ui.tabControl();
resultsDiv.addPane('FLOW BY TIME', () => DG.Viewer.fromType('Line chart', data).root);
resultsDiv.addPane('DATA', () => DG.Viewer.grid(data).root);
resultsDiv.root.style.height = '400px';

view.append(ui.divV([
  top,
  ui.div([ui.bigButton('SIMULATE FLOW', () => progressDiv.innerText = 'Calculating')]),
  resultsDiv.root
]));