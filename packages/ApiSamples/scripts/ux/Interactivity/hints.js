const df = grok.data.demo.demog();
const view = grok.shell.addTableView(df);

// Place a hint on a visual component in your application
const addNCIcon = $('div.d4-ribbon-item').has('i.svg-add-new-column')[0];
ui.hints.addHint(addNCIcon);

// Remove it once certain conditions are met, e.g., mouse events, platform events
// available via the grok.events entrypoint (use rxjs.operators.first(), if needed).
// Hints can be removed both programmatically and by the user (via the close icon).
$(addNCIcon).on('click', () => {
  ui.hints.remove(addNCIcon);
});

grok.events.onDialogShown.pipe(
  rxjs.operators.filter((dlg) => dlg.title === 'Add New Column'),
  rxjs.operators.first())
  .subscribe((dlg) => {
    const input = dlg.inputs[0].root;
    ui.hints.addHint(input);
    setTimeout(() => ui.hints.remove(input), 4000);
  });

// You can also describe a series of visual components in the wizard. Each wizard page is associated
// with the [showNextTo] element. Provide either [text] or [root] parameter to populate the page, the
// other parameters are optional. The wizard header is shown only if [title] or [helpUrl] are provided.
// The user can use the arrow buttons to navigate the set of instructions. The wizard can be closed
// from the dialog header (the "x" icon), via the "Cancel" button, or via the [Wizard.close()] method.
const sp = view.scatterPlot();
const hist = view.histogram();

const wizard = ui.hints.addTextHint({title: 'Viewers', pages: [
 {caption: 'First Hint', text: 'Welcome! This is the 1st hint about a scatter plot.', showNextTo: sp.root},
 {caption: 'Second Hint', root: ui.divText('This is the 2nd hint about a histogram.'),  showNextTo: hist.root},
]});
