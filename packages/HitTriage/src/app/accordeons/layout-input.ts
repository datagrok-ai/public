import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

export function getLayoutInput() {
  const dataFileProp = DG.Property.fromOptions({name: 'layout', caption: 'Layout', type: 'file'});
  const dataFileInput = DG.InputBase.forProperty(dataFileProp);
  let layout: DG.ViewLayout | null = null;
  let layoutViewState: string | null = null;
  // if the file is local, we can store it in a variable
  let localFilePath: string | null = null;
  let error: string | null = null;
  dataFileInput.onChanged.subscribe(async (_) => {
    const f: DG.FileInfo = dataFileInput.value;
    if (f) {
      try {
        const exists = f.fullPath ? await grok.dapi.files.exists(f.fullPath) : false;
        localFilePath = exists ? f.fullPath : null;
      } catch (e) {
        localFilePath = null;
      }
      try {
        const contentString = await f.readAsString();
        const l = DG.ViewLayout.fromJson(contentString);
        if (l) {
          layout = l;
          layoutViewState = l.viewState;
          error = null;
        } else {
          layout = null;
          layoutViewState = null;
          localFilePath = null;
          error = 'Invalid layout file';
        }
      } catch (e) {
        layout = null;
        layoutViewState = null;
        localFilePath = null;
        error = 'Invalid layout file';
      }
      dataFileInput.input.classList[(error ? 'add' : 'remove')]('d4-invalid');
    } else {
      layout = null;
      layoutViewState = null;
      localFilePath = null;
    }
  });
  ui.tooltip.bind(dataFileInput.root, () => error ? error : 'Select a layout file');
  return {
    dataFileInput, getLayout: () => layout, getLayoutViewState: () => layoutViewState,
    getLocalFilePath: () => localFilePath,
  };
}
