import {after, before, category, delay, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isViewerPresent, waitForElement} from '../gui-utils';


category('Viewers: Leaflet', () => {
  let tv: DG.TableView;
  let demog: DG.DataFrame;

  before(async () => {
    demog = grok.data.demo.demog(1000);
    tv = grok.shell.addTableView(demog);
  });

  
  test('leaflet.visual', async () => {
    tv.addViewer('Leaflet');
    await waitForElement('.d4-leaflet', 'cannot find Leaflet');
    isViewerPresent(Array.from(tv.viewers), 'Leaflet');
      
    const leaflet = document.querySelector('[name=viewer-Leaflet]') as HTMLElement;
    const zoomIn = leaflet.getElementsByClassName('leaflet-control-zoom-in')[0] as HTMLElement;
    const zoomOut = leaflet.getElementsByClassName('leaflet-control-zoom-out')[0] as HTMLElement;
    
    zoomIn.click();
    await delay(500);
    zoomOut.click();
    await delay(500);
  });


  test('leaflet.serialization', async () => {
    tv.addViewer('Leaflet');
    await waitForElement('.d4-leaflet', 'cannot find Leaflet');
    const layout = tv.saveLayout();
    tv.resetLayout();
    tv.loadLayout(layout);
    await waitForElement('.d4-leaflet', 'cannot find Leaflet');
    isViewerPresent(Array.from(tv.viewers), 'Leaflet');
  });


  after(async () => {
    tv.close();
    grok.shell.closeTable(demog);
  }); 
});
