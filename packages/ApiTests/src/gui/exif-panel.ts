import {category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {checkHTMLElementbyInnerText, getHTMLElementbyInnerText} from './gui-utils';

category('File Panels: EXIF', () => {

  test('panel.exif', async () => {
    let pictures = await grok.dapi.files.list('Demo:Files/images', true, "jpg");

    let picCar1;
    for(let i = 0; i < pictures.length; i++){
        if (pictures[i].name == "car1.jpg"){
            picCar1 = pictures[i];
            break;
        }
    }

    grok.shell.o = picCar1; await delay(500);

    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'EXIF');

    let exifPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'EXIF')
    exifPanel!.click(); await delay(3000);

    if (document.getElementsByClassName('d4-pane-exif expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in EXIF Panel'  

    if (document.getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-exif expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'table with EXIF content was not rendered in the panel'  

    exifPanel!.click(); await delay(500);
  });
});
