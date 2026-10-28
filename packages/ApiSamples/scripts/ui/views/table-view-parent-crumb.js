//api: DG.TableView.showParentCrumb, DG.ViewBase.parentView
// An app's table view with the app's home as a crumb: the ribbon reads "Demo app › demog".
// A click on "Demo app" makes it current; closing the table goes back to it.

const home = grok.shell.newView('Demo app', [ui.divText('The app home')]);
const tv = grok.shell.addTableView(grok.data.demo.demog());
tv.parentView = home;
tv.showParentCrumb = true;
