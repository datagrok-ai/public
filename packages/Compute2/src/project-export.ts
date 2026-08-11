import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

export interface DfExportEntry {
  name: string;
  isInput: boolean;
  df: DG.DataFrame;
  viewers: {type: string, options: Record<string, any>}[];
}

// Project.showSaveDialog is newer than the last released js-api, so it is reached via a
// cast and feature-detected at runtime; older platforms may still have the Dart binding
export const canSaveProject = () =>
  typeof (DG.Project as any).showSaveDialog === 'function' ||
  typeof (window as any).grok_Project_OpenSaveDialog === 'function';

const NON_SYNCABLE_INPUT_TYPES = new Set<string>([
  DG.TYPE.DATA_FRAME, DG.TYPE.DATA_FRAME_LIST, DG.TYPE.FILE, DG.TYPE.BLOB,
  DG.TYPE.COLUMN, DG.TYPE.COLUMN_LIST, DG.TYPE.GRAPHICS, DG.TYPE.OBJECT, DG.TYPE.DYNAMIC,
]);

/** The producing call serialized by the platform itself (`prepare().toString()`), or null
 *  when the call cannot be replayed from literals — outputs then publish as static snapshots. */
function serializeProducingCall(call: DG.FuncCall): string | null {
  const inputProps = [...call.inputParams.values()].map((p) => p.property);
  if (inputProps.some((p) => NON_SYNCABLE_INPUT_TYPES.has(p.propertyType as string)))
    return null;
  try {
    const inputs: Record<string, any> = {};
    for (const prop of inputProps) {
      const value = call.inputs[prop.name];
      if (value != null)
        inputs[prop.name] = value;
    }
    return call.func.prepare(inputs).toString();
  } catch (e) {
    console.warn('Compute2: could not serialize the producing call', e);
    return null;
  }
}

/** Stamp each output clone with the producing call (`.script` + `.VariableName` tags) so the
 *  core Save-project dialog offers data sync — identical producing calls dedup on project open,
 *  so the function runs once and every table binds its own output. */
function stampCreationScripts(call: DG.FuncCall, entries: DfExportEntry[]) {
  const outputs = entries.filter((e) => !e.isInput);
  if (outputs.length === 0)
    return;
  const callStr = serializeProducingCall(call);
  if (callStr == null)
    return;
  const accessor = [...call.outputParams.values()].length > 1;
  let ts = Date.now();
  for (const entry of outputs) {
    entry.df.setTag(DG.Tags.VariableName, entry.name);
    entry.df.setTag(DG.Tags.CreationScript,
      `${entry.name} = ${callStr}${accessor ? '.' + entry.name : ''} //{"timestamp": ${ts++}}`);
  }
}

const showProjectSaveDialog = (
  tables: DG.DataFrame[], views: DG.TableView[], name: string, description: string,
): Promise<DG.Project | null> => {
  const jsApiShowSaveDialog = (DG.Project as any).showSaveDialog;
  if (typeof jsApiShowSaveDialog === 'function')
    return jsApiShowSaveDialog.call(DG.Project, {tables, views, name, description});
  return (window as any).grok_Project_OpenSaveDialog(
    tables.map((t) => t.dart), views.map((v) => v.dart), [], name, description, '');
};

/** Publishes the visible dataframes of a run as a project: one detached TableView per
 *  dataframe with its annotated viewers (platform default docking), handed to the standard
 *  Save-project dialog, which does upload, data-sync toggles, layout linking, and sharing. */
export async function saveCallToProject(call: DG.FuncCall, entries: DfExportEntry[]): Promise<void> {
  if (entries.length === 0) {
    grok.shell.warning('Nothing to export: the run has no visible dataframes');
    return;
  }
  const usedNames = new Set<string>();
  const cloned = entries.map((entry) => {
    const df = entry.df.clone();
    let name = entry.df.name || entry.name;
    for (let i = 2; usedNames.has(name); i++)
      name = `${entry.df.name || entry.name} (${i})`;
    usedNames.add(name);
    df.name = name;
    return {...entry, df};
  });
  stampCreationScripts(call, cloned);

  // Views must be attached to the shell: a never-attached view has an uninitialized dock
  // layout that serializes as '{}' and crashes DockGraphDeserializer on project open
  const currentView = grok.shell.v;
  const views = cloned.map((entry) => {
    const view = DG.TableView.create(entry.df, false);
    for (const {type, options} of entry.viewers) {
      if (type !== DG.VIEWER.GRID)
        view.addViewer(type, options);
    }
    grok.shell.addView(view);
    return view;
  });
  if (currentView)
    grok.shell.v = currentView;
  const name = call.options['title'] ?? call.func.friendlyName ?? call.func.name;
  try {
    const saved = await showProjectSaveDialog(cloned.map((e) => e.df), views, name, call.func.description ?? '');
    if (saved?.id)
      grok.shell.info(`Project "${saved.friendlyName ?? name}" saved`);
  } catch (e: any) {
    grok.shell.error(e);
  } finally {
    for (const view of views)
      view.close();
    for (const entry of cloned)
      grok.shell.closeTable(entry.df);
  }
}
