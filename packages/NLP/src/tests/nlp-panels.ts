import {before, category, delay, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import {FileInfo} from 'datagrok-api/src/entities';
import {checkHTMLElementByInnerText} from './test-utils';


category('File Panels: NLP', () => {
  const testFiles: FileInfo[] = [];

  before(async () => {
    const txtFiles = await grok.dapi.files.list('Demo:Files/texts', true, 'txt');
    const pdfFiles = await grok.dapi.files.list('Demo:Files/texts', true, 'pdf');
    const docFiles = await grok.dapi.files.list('Demo:Files/texts', true, 'doc');
    const txtTestFile = txtFiles.find((f) => f.name === 'dart.txt');
    const pdfTestFile = pdfFiles.find((f) => f.name === 'da-sdg.pdf');
    const docTestFile = docFiles.find((f) => f.name === 'en-sdg.doc');
    [txtTestFile, pdfTestFile, docTestFile].forEach((f) => {if (f) testFiles.push(f);});
  });

  test('nlp.textStatistics', () => checkPane('Text Statistics'));

  test('nlp.Summary', () => checkPane('Summary'));

	async function checkPane(name: string) {
    if (!testFiles.length)
      throw 'Failed to find files for the test';

    for (const testFile of testFiles) {
      grok.shell.o = testFile;
      await delay(500);
      const pane = checkHTMLElementByInnerText('d4-accordion-pane-header', name)!;
      pane.click();
      await delay(2000);
      if (pane.querySelector('.d4-error'))
        throw `Error while processing ${testFile.name}`;
      else if (!pane.querySelector('.d4-table.d4-item-table.d4-info-table'))
        throw `The panel content was not rendered for file ${testFile.name}`;
			pane.click();
    }
  }
});
