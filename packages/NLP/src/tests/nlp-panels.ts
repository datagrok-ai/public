import {before, category, delay, test} from '@datagrok-libraries/utils/src/test';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {checkHTMLElementByInnerText} from './test-utils';
import {_package} from '../package-test';


category('File Panels: NLP', () => {
  const testFiles: DG.FileInfo[] = [];

  before(async () => {
    testFiles.push(...(await _package.files.list('', true)));
  });

  test('nlp.textStatistics', () => checkPane('Text Statistics'), {skipReason: 'GROK-11391'});

  test('nlp.Summary', () => checkPane('Summary'), {skipReason: 'GROK-11391'});

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
