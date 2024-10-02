import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

export function getLayoutInput() {
  const dataFileProp = DG.Property.fromOptions({name: 'layout', caption: 'Layout', type: 'file'});
  const dataFileInput = DG.InputBase.forProperty(dataFileProp);
  let layout: DG.ViewLayout | null = null;
  let layoutViewState: string | null = null;
  let error: string | null = null;
  dataFileInput.addValidator((_) => error);
  dataFileInput.onChanged.subscribe(async (value) => {
    const f: DG.FileInfo = value;
    if (f) {
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
          error = 'Invalid layout file';
        }
      } catch (e) {
        layout = null;
        layoutViewState = null;
        error = 'Invalid layout file';
      }
    }
    dataFileInput.setTooltip('Select a layout file');
  });
  return {
    dataFileInput, getLayout: () => layout, getLayoutViewState: () => layoutViewState,
  };
}
