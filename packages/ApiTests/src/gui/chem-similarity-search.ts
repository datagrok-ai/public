import {after, before, category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElement} from '../ui/utils';
import {isColumnPresent, isViewerPresent, isDialogPresent, returnDialog, setDialogInputValue} from './gui-utils';
import { Viewer } from 'datagrok-api/dg';

category('Chem: Similarity Search', () => {
  let v: DG.TableView;
  const smiles = grok.data.demo.molecules(20);

  before(async () => {
    v = grok.shell.addTableView(smiles);
  });

  test('chem.similaritySearch', async () => {
    
    smiles.currentRowIdx = 0;

    grok.shell.topMenu.find('Chem').find('Similarity Search...').root.click(); await delay(2000);

    isViewerPresent(Array.from(v.viewers), 'SimilaritySearchViewer');

    let structuresInViewer = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('d4-flex-wrap ui-div')[0] as HTMLElement;
    if (structuresInViewer.childElementCount != 10)
        throw 'molecules number inside the viewer is different than expected'; 

    let similarityViewer:DG.Viewer;
    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'SimilaritySearchViewer') {
            similarityViewer = Array.from(v.viewers)[i];
            break;
        }
    }

    similarityViewer!.props.distanceMetric = 'Dice';
    similarityViewer!.props.limit = 5    

    let similarityLable = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('ui-label')[0] as HTMLElement;
    if (similarityLable.innerText != '0.22')
        throw 'Expected Similarity Lable for 2nd molecule does not match the "Dice" metric';
    
        structuresInViewer = document.getElementsByClassName('d4-layout-root d4-root d4-viewer d4-similaritysearchviewer ui-box')[0].getElementsByClassName('d4-flex-wrap ui-div')[0] as HTMLElement;
    if (structuresInViewer.childElementCount != 5)
        throw 'molecules number inside Similarity viewer is different than expected after change "Limit" property'; 

    let protpertiesBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-settings')[0] as HTMLElement;
    protpertiesBtn.click(); await delay(1000);
    if (document.getElementsByClassName('property-grid-base property-grid-disable-selection').length == 0)
        throw 'Properties table does not open'

    let closeBtn = document.getElementsByClassName('panel-titlebar disable-selection panel-titlebar-tabhost')[0].getElementsByClassName('grok-icon grok-font-icon-close')[0] as HTMLElement;
    closeBtn.click(); await delay(1000);

    for (let i:number = 0; i < Array.from(v.viewers).length; i++) {
        if (Array.from(v.viewers)[i].type == 'SimilaritySearchViewer') {
            throw 'SimilaritySearch viewer was not closed'
        }
    }
  });  
  after(async () => {
    grok.shell.closeAll();
  }); 
});
