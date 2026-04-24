// FuncCall.getResultViews() — views produced from the call's dataframe/graphics outputs.
// Call on a succeeded FuncCall; returns a TableView per dataframe output (new or existing).

const fc = DG.Func.byName('OpenServerFile')
  .prepare({fullPath: 'System:DemoFiles/demog.csv'});
await fc.call();

const resultView = fc.getResultViews()[0];

// Embed the result view's root as an inner element of a regular custom view
const host = grok.shell.newView('Result view host', [
  ui.h2(`Embedded: ${resultView.name}`),
  ui.box(resultView.root),
]);
