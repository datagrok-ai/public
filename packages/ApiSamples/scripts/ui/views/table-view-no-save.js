//api: DG.TableView.showSaveButton
// A table view whose app saves the data its own way: no project Save button on the ribbon,
// the rest of the ribbon stays, and the app adds its own actions.

const tv = grok.shell.addTableView(grok.data.demo.demog());
tv.showSaveButton = false;
tv.setRibbonPanels([...tv.getRibbonPanels(), [ui.bigButton('Submit', () => grok.shell.info('Submitted'))]]);
