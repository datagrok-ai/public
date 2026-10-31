import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {DfExportEntry, prepareExportEntries} from './project-export-data';
import {getRunName} from '@datagrok-libraries/compute-utils/shared-utils/utils';
import {modelLabel} from './components/RunStorage/run-storage';

export type {DfExportEntry};

export const showProjectSaveDialog = (
  tables: DG.DataFrame[], layouts: string[], name: string, description = '',
): Promise<DG.Project | null> => DG.Project.showSaveDialog({tables, layouts, name, description});

// The save dialog gets serialized layouts, not live views: detached view roots are not in
// the document, so the dialog's thumbnail capture fails on them (zero-height canvas). The
// view still must be attached to the shell while the layout is taken: a never-attached view
// has an uninitialized dock layout that serializes as '{}' and crashes DockGraphDeserializer
// on project open, and view.grid is null.
export function captureTableLayout(df: DG.DataFrame, viewers: {type: string, options: Record<string, any>}[]): string {
  const currentView = grok.shell.v;
  const view = DG.TableView.create(df, false);
  grok.shell.addView(view);
  for (const {type, options} of viewers) {
    if (type === DG.VIEWER.GRID)
      view.grid.setOptions(options);
    else
      view.addViewer(type, options);
  }
  // the dialog's layouts parameter takes ViewLayout.viewState strings, not toJson()
  const layout = view.saveLayout().viewState;
  view.close();
  if (currentView)
    grok.shell.v = currentView;
  return layout;
}

/** Publishes the visible dataframes of a run as a project: one layout per dataframe with
 *  its viewers (platform default docking), handed to the standard Save-project dialog,
 *  which does upload, data-sync toggles, layout linking, and sharing. */
export async function saveCallToProject(call: DG.FuncCall, entries: DfExportEntry[]): Promise<void> {
  if (entries.length === 0) {
    grok.shell.warning('Nothing to export: the run has no visible dataframes');
    return;
  }
  const cloned = prepareExportEntries(call, entries);
  const layouts = cloned.map((entry) => captureTableLayout(entry.df, entry.viewers));
  const name = getRunName(call) || modelLabel(call.func);
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
