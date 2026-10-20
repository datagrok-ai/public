import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {after, awaitCheck, before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {DEFAULT_FORMATS} from '../apps/common/model/const';
import {convert} from '../apps/translator/model/conversion-utils';
import {NUCLEOTIDES_FORMAT} from '../apps/translator/view/const';
import {ITranslationHelper} from '../types';

import {_package} from '../package-test';

category('Translator: bulk', () => {
  let th: ITranslationHelper;

  before(async () => {
    th = await _package.getTranslationHelper();
  });

  test('bulk-translation-axolabs.csv', async () => {
    const df = DG.DataFrame.fromCsv(await _package.files.readAsText('samples/bulk-translation-axolabs.csv'));
    // the sixth row is empty
    const sequences: string[] = df.getCol('AxolabsSequences').toList().slice(0, 5);
    expect(sequences.map((s) => convert(s, DEFAULT_FORMATS.AXOLABS, NUCLEOTIDES_FORMAT, th)).join(', '),
      'ACGUACGU, UCGUACGU, CCGUACGU, ACUACGU, ACGUACGU');
    const helm = sequences.map((s) => convert(s, DEFAULT_FORMATS.AXOLABS, DEFAULT_FORMATS.HELM, th));
    expect(helm.every((h) => h?.startsWith('RNA1{')), true, helm.join(' | '));
  });
});

category('Translator: app', () => {
  after(async () => {
    grok.shell.closeAll();
  });

  test('A HELM with a monomer the oligo library lacks refreshes the translations (GROK-20958)', async () => {
    const view: DG.ViewBase = await grok.functions.call(`${_package.name}:oligoTranslatorApp`);
    grok.shell.addView(view);
    const textArea = view.root.querySelector('textarea')!;
    // the first row is the FORMAT / SEQUENCE header
    const translations = (): string => Array.from(view.root.querySelectorAll('.st-translator-output-table tr')).slice(1)
      .map((tr) => Array.from(tr.querySelectorAll('td')).map((td) => td.textContent).join('=')).join(', ');
    const type = async (sequence: string, shown: (translations: string) => boolean): Promise<void> => {
      textArea.value = sequence;
      textArea.dispatchEvent(new Event('input'));
      await awaitCheck(() => shown(translations()), '', 5000).catch(() => {
        throw new Error(`the translations of ${sequence}: ${translations()}`);
      });
    };

    await type('RNA1{r(A)p.r(C)p.r(G)p.r(U)}$$$$', (t) => t.includes('Nucleotides=ACGU'));
    await type('RNA1{r(A)p.r(C)p.[meI]}$$$$', (t) => !t.includes('ACGU'));
    expect(translations(), 'Axolabs=<?><?><?>, BioSpring=AC<?>, Mermade12=<?><?><?>, Nucleotides=AC');
  });
});
