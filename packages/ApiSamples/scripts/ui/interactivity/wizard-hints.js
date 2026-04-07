const df = grok.data.demo.demog();
const view = grok.shell.addTableView(df);

// You can describe a series of visual components in the wizard. Each wizard page is associated
// with the [showNextTo] element. Provide either [text] or [root] parameter to populate the page, the
// other parameters are optional. The wizard header is shown only if [title] or [helpUrl] are provided.
// The user can use the arrow buttons to navigate the set of instructions. The wizard can be closed
// from the dialog header (the "x" icon), via the "Cancel" button, or via the [Wizard.close()] method.
const sp = view.scatterPlot();
const hist = view.histogram();

const wizard = ui.hints.addTextHint({title: 'Viewers', pages: [
  {caption: 'First Hint', text: 'Welcome! This is the 1st hint about a scatter plot.', showNextTo: sp.root},
  {caption: 'Second Hint', root: ui.divText('This is the 2nd hint about a histogram.'), showNextTo: hist.root},
]});
