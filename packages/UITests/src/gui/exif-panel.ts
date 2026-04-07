import {category, delay, test} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import {checkHTMLElementbyInnerText, getHTMLElementbyInnerText} from './gui-utils';

category('File Panels: EXIF', () => {
  test('panel.exif', async () => {
    const pictures = await grok.dapi.files.list('Demo:Files/images', true, 'jpg');
    let picCar1;
    for (let i = 0; i < pictures.length; i++) {
      if (pictures[i].name == 'car1.jpg') {
        picCar1 = pictures[i];
        break;
      }
    }
    grok.shell.o = picCar1; await delay(500);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'EXIF');
    const exifPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'EXIF');
    exifPanel!.click(); await delay(3000);

    if (document.getElementsByClassName('d4-pane-exif expanded')[0].getElementsByClassName('d4-error').length != 0)
      throw new Error('Error in EXIF Panel');

    if (document
      .getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-exif expanded')[0]
      .getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
      throw new Error('table with EXIF content was not rendered in the panel');

    exifPanel!.click(); await delay(500);
  });
});
