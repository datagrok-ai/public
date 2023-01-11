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
  rxjs.operators.first(),
  rxjs.operators.filter((dlg) => dlg.title === 'Add New Column'))
  .subscribe((dlg) => {
    const input = dlg.inputs[0].root;
    ui.hints.addHint(input);
    setTimeout(() => ui.hints.remove(input), 4000);
  });
