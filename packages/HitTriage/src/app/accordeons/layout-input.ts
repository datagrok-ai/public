import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

/** A template's layout: the view state of a .layout file, and the file's path when it lives on a share, where
 * "Apply template layout" reads it again. */
export interface TemplateLayout {
  viewState: string;
  localPath: string | null;
}

/** The platform's file input a template's layout is picked with, from the computer or a share. */
export function layoutFileInput(): DG.InputBase<DG.FileInfo | null> {
  return DG.InputBase.forProperty(DG.Property.fromOptions({name: 'layout', caption: 'Layout', type: 'file'}));
}

/** The layout `file` holds; null when it holds none. */
export async function readLayout(file: DG.FileInfo): Promise<TemplateLayout | null> {
  let layout: DG.ViewLayout | null;
  try {
    layout = DG.ViewLayout.fromJson(await file.readAsString());
  }
  catch {
    return null;
  }
  if (!layout)
    return null;
  // a share's path is namespaced by its connection ("System:AppData/…"); a file from the computer or memory is not
  const onShare = !!file.fullPath?.includes(':') && await grok.dapi.files.exists(file.fullPath).catch(() => false);
  return {viewState: layout.viewState, localPath: onShare ? file.fullPath : null};
}

/** A file in memory holding `layout`, under its share path or `<name>.layout`: what a file input shows for a layout
 * that is already a template's, read back by {@link readLayout} as it was. */
export function layoutFile(layout: TemplateLayout, name: string): DG.FileInfo {
  return DG.FileInfo.fromString(layout.localPath ?? `${name}.layout`,
    DG.ViewLayout.fromViewState(layout.viewState).toJson());
}

/** The layout input of the Hit Triage and PepTriage editors; a clone's layout (`preset`) is there from the start. */
export function getLayoutInput(preset?: {name: string, layoutViewState?: string, localLayoutPath?: string}) {
  const dataFileInput = layoutFileInput();
  let layout: TemplateLayout | null = preset?.layoutViewState ?
    {viewState: preset.layoutViewState, localPath: preset.localLayoutPath ?? null} : null;
  let error: string | null = null;
  let reads = 0;
  dataFileInput.onChanged.subscribe(async (_) => {
    const read = ++reads;
    const f: DG.FileInfo | null = dataFileInput.value;
    const found = f ? await readLayout(f) : null;
    if (read !== reads)
      return;
    layout = found;
    error = f && !found ? 'Invalid layout file' : null;
    dataFileInput.input.classList.toggle('d4-invalid', error !== null);
  });
  ui.tooltip.bind(dataFileInput.root, () => error ? error : 'Select a layout file');
  if (layout !== null)
    dataFileInput.value = layoutFile(layout, preset!.name);
  return {dataFileInput, getLayoutViewState: () => layout?.viewState ?? null};
}
