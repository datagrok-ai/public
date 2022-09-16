import {category, delay, expect, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {checkHTMLElementbyInnerText, getHTMLElementbyInnerText} from './gui-utils';

category('File Panels: NLP', () => {
    
  test('nlp.textStatistics', async () => {
    let txtFiles = await grok.dapi.files.list('Demo:Files/texts', true, "txt");
    let pdfFiles = await grok.dapi.files.list('Demo:Files/texts', true, "pdf");
    let docFiles = await grok.dapi.files.list('Demo:Files/texts', true, "pdf");
    
    let txtTestFile;
    for(let i = 0; i < txtFiles.length; i++){
        if (txtFiles[i].name == "datr.txt"){
            txtTestFile = txtFiles[i];
            break;
        }
    }

    let pdfTestFile;
    for(let i = 0; i < pdfFiles.length; i++){
        if (pdfFiles[i].name == "da-sdg.pdf"){
            pdfTestFile = pdfFiles[i];
            break;
        }
    }

    let docTestFile;
    for(let i = 0; i < docFiles.length; i++){
        if (docFiles[i].name == "datr.txt"){
            docTestFile = docFiles[i];
            break;
        }
    }

    grok.shell.o = txtTestFile; await delay(500);

    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Text Statistics');

    let textStatisticsPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Text Statistics')
    textStatisticsPanel!.click(); await delay(5000);

    if (document.getElementsByClassName('d4-pane-exif expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in Text Statistics Panel for TXT file'  

    if (document.getElementsByClassName('d4-pane-text_statistics expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Text Statistic panel content was not rendered for TXT file'  

    grok.shell.o = pdfTestFile; await delay(5000);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Text Statistics');

    if (document.getElementsByClassName('d4-pane-exif expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in Text Statistics Panel for PDF file'  

    if (document.getElementsByClassName('d4-pane-text_statistics expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Text Statistic panel content was not rendered for PDF file'  

    grok.shell.o = docTestFile; await delay(5000);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Text Statistics');
    
    if (document.getElementsByClassName('d4-pane-exif expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in Text Statistics Panel for DOC file'  

    if (document.getElementsByClassName('d4-pane-text_statistics expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Text Statistic panel content was not rendered for DOC file'  

    textStatisticsPanel!.click(); await delay(500);
  });

  test('nlp.Summary', async () => {
    let txtFiles = await grok.dapi.files.list('Demo:Files/texts', true, "txt");
    let pdfFiles = await grok.dapi.files.list('Demo:Files/texts', true, "pdf");
    let docFiles = await grok.dapi.files.list('Demo:Files/texts', true, "pdf");
    
    let txtTestFile;
    for(let i = 0; i < txtFiles.length; i++){
        if (txtFiles[i].name == "datr.txt"){
            txtTestFile = txtFiles[i];
            break;
        }
    }

    let pdfTestFile;
    for(let i = 0; i < pdfFiles.length; i++){
        if (pdfFiles[i].name == "da-sdg.pdf"){
            pdfTestFile = pdfFiles[i];
            break;
        }
    }

    let docTestFile;
    for(let i = 0; i < docFiles.length; i++){
        if (docFiles[i].name == "datr.txt"){
            docTestFile = docFiles[i];
            break;
        }
    }

    grok.shell.o = txtTestFile; await delay(500);

    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Summary');

    let summaryPanel = getHTMLElementbyInnerText('d4-accordion-pane-header', 'Summary')
    summaryPanel!.click(); await delay(5000);

    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in Summary Panel for TXT file'  

    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Summary panel content was not rendered for TXT file'  

    grok.shell.o = pdfTestFile; await delay(5000);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Summary');

    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Error in Summary Panel for PDF file'  

    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Summary panel content was not rendered for PDF file'  

    grok.shell.o = docTestFile; await delay(5000);
    checkHTMLElementbyInnerText('d4-accordion-pane-header', 'Summary');
    
    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-error').length != 0)
        throw 'Summary Panel for DOC file'  

    if (document.getElementsByClassName('d4-pane-summary expanded')[0].getElementsByClassName('d4-table d4-item-table d4-info-table').length != 1)
        throw 'Summary content was not rendered for DOC file'  

    summaryPanel!.click(); await delay(500);
  });  
});
