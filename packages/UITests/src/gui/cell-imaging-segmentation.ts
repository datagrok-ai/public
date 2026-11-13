import {category, test, waitUntil} from '@datagrok-libraries/test/src/test';
import * as grok from 'datagrok-api/grok';
import {checkHTMLElementbyInnerText, getHTMLElementbyInnerText} from './gui-utils';

category('File Panels: Cell Imaging Segmentation', () => {
  test('panel.cellImagingSegmentation', async () => {
    const pictures = await grok.dapi.files.list('Demo:Files/images', true, 'jpeg');

    let picCell1;
    for (let i = 0; i < pictures.length; i++) {
      if (pictures[i].name == 'cells1.jpeg') {
        picCell1 = pictures[i];
        break;
      }
    }

    grok.shell.o = picCell1;
    await waitUntil(() => getHTMLElementbyInnerText('d4-accordion-pane-header',
      'Cell Imaging Segmentation') !== undefined, 'Cell Imaging Segmentation panel not found');

    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Cell Imaging Segmentation');

    const cellsPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Cell Imaging Segmentation');
    cellsPanel!.click();
    await waitUntil(() => {
      const pane = document
        .getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-cell_imaging_segmentation expanded')[0];
      return pane !== undefined && (pane.getElementsByClassName('d4-error').length != 0 ||
        pane.getElementsByClassName('grok-scripting-image-container-info-panel').length != 0);
    }, 'Cell Imaging Segmentation panel content was not loaded');

    if (document
      .getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-cell_imaging_segmentation expanded')[0]
      .getElementsByClassName('d4-error').length != 0)
      throw new Error('Error in Cell Imaging Segmentation Panel');

    if (document
      .getElementsByClassName('d4-accordion-pane-content ui-div d4-pane-cell_imaging_segmentation expanded')[0]
      .getElementsByClassName('grok-scripting-image-container-info-panel').length != 1)
      throw new Error('Cell Imaging Segmentation content was not rendered in the panel');

    cellsPanel!.click();
  });
});
