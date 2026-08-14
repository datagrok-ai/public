import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DfExportEntry, prepareExportEntries} from './project-export-data';

export {DfExportEntry};

// Project.showSaveDialog is newer than the last released js-api, so it is reached via a
// cast and feature-detected at runtime; older platforms may still have the Dart binding
export const canSaveProject = () =>
  typeof (DG.Project as any).showSaveDialog === 'function' ||
  typeof (window as any).grok_Project_OpenSaveDialog === 'function';

const showProjectSaveDialog = (
  tables: DG.DataFrame[], layouts: string[], name: string, description: string,
): Promise<DG.Project | null> => {
  const jsApiShowSaveDialog = (DG.Project as any).showSaveDialog;
  if (typeof jsApiShowSaveDialog === 'function')
    return jsApiShowSaveDialog.call(DG.Project, {tables, layouts, name, description});
  return (window as any).grok_Project_OpenSaveDialog(
    tables.map((t) => t.dart), [], layouts, name, description, '');
};

/** Publishes the visible dataframes of a run as a project: one layout per dataframe with
 *  its viewers (platform default docking), handed to the standard Save-project dialog,
 *  which does upload, data-sync toggles, layout linking, and sharing. */
export async function saveCallToProject(call: DG.FuncCall, entries: DfExportEntry[]): Promise<void> {
  if (entries.length === 0) {
    grok.shell.warning('Nothing to export: the run has no visible dataframes');
    return;
  }
  const cloned = prepareExportEntries(call, entries);

  // The dialog gets serialized layouts, not live views: detached view roots are not in the
  // document, so the dialog's html2canvas thumbnail capture fails on them ("Unable to find
  // element in cloned iframe"). Views still must be attached to the shell while the layout
  // is taken: a never-attached view has an uninitialized dock layout that serializes as
  // '{}' and crashes DockGraphDeserializer on project open, and view.grid is null.
  const currentView = grok.shell.v;
  const layouts = cloned.map((entry) => {
    const view = DG.TableView.create(entry.df, false);
    grok.shell.addView(view);
    for (const {type, options} of entry.viewers) {
      if (type === DG.VIEWER.GRID)
        view.grid.setOptions(options);
      else
        view.addViewer(type, options);
    }
    // the dialog's layouts parameter takes ViewLayout.viewState strings, not toJson()
    const layout = view.saveLayout().viewState;
    view.close();
    return layout;
  });
  if (currentView)
    grok.shell.v = currentView;
  const name = call.options['title'] ?? call.func.friendlyName ?? call.func.name;
  try {
    const saved = await showProjectSaveDialog(cloned.map((e) => e.df), layouts, name, call.func.description ?? '');
    if (saved?.id)
      grok.shell.info(`Project "${saved.friendlyName ?? name}" saved`);
  } catch (e: any) {
    grok.shell.error(e);
  } finally {
    for (const entry of cloned)
      grok.shell.closeTable(entry.df);
  }
}
