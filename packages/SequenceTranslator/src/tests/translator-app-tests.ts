import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {after, waitUntil, category, expect, test} from '@datagrok-libraries/test/src/test';

import {DEFAULT_FORMATS} from '../apps/common/model/const';
import {DEFAULT_AXOLABS_INPUT} from '../apps/common/view/const';
import {_package} from '../package-test';
import {getHelm} from './utils';

/** The `<select>` of the input captioned [caption] in [root]. */
function selectOf(root: HTMLElement, caption: string): HTMLSelectElement {
  const input = $(root).find('.ui-input-root').filter((_, e) => $(e).find('label').first().text().trim() === caption)[0];
  if (!input)
    throw new Error(`the view has no "${caption}" input`);
  return $(input).find('select')[0] as HTMLSelectElement;
}

/** Picks the option shown as [text] in [select], as a click on it does. */
function choose(select: HTMLSelectElement, text: string): void {
  const option = Array.from(select.options).find((o) => o.text === text);
  if (!option)
    throw new Error(`no "${text}" among ${Array.from(select.options).map((o) => o.text).join(', ')}`);
  select.selectedIndex = option.index;
  select.dispatchEvent(new Event('input', {bubbles: true}));
  select.dispatchEvent(new Event('change', {bubbles: true}));
}

async function openTranslator(): Promise<DG.ViewBase> {
  const view: DG.ViewBase = await grok.functions.call('SequenceTranslator:oligoTranslatorApp');
  grok.shell.addView(view);
  await waitUntil(() => $(view.root).find('textarea').length > 0, 'the Translator did not appear', 5000);
  return view;
}

function convertIn(view: DG.ViewBase): void {
  $(view.root).find('button').filter((_, b) => b.textContent === 'Convert').trigger('click');
}

function axolabsTable(name: string): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.fromStrings('AxolabsSequences', [DEFAULT_AXOLABS_INPUT])]);
  df.name = name;
  return df;
}

category('Oligo Translator', () => {
  after(async () => {
    DG.Balloon.closeAll();
    grok.shell.closeAll();
  });

  test('Get SDF saves the drawn structure of a HELM input', async () => {
    const th = await _package.getTranslationHelper();
    const view: DG.ViewBase = await grok.functions.call('SequenceTranslator:oligoTranslatorApp');
    grok.shell.addView(view);
    await waitUntil(() => $(view.root).find('textarea').length > 0, 'Sequence input did not appear', 5000);

    const textArea = $(view.root).find('textarea')[0] as HTMLTextAreaElement;
    textArea.value = getHelm(DEFAULT_AXOLABS_INPUT, DEFAULT_FORMATS.AXOLABS, th);
    textArea.dispatchEvent(new Event('input'));
    await waitUntil(() => ($(view.root).find('select')[0] as HTMLSelectElement).value === DEFAULT_FORMATS.HELM,
      'HELM input was not detected', 5000);

    DG.Balloon.closeAll();
    $(view.root).find('button').filter((_, b) => b.textContent === 'Get SDF').trigger('click');
    expect($('.d4-balloon.warning').length, 0);
  });
  // The bulk choices lived in one EventBus for the whole page: a Translator converted with the formats, the table or
  // the column another one had picked, while showing its own.
  test('Two Translators keep their own bulk choices: the formats and the table', async () => {
    const a = await openTranslator();
    const b = await openTranslator();
    const first = axolabsTable('bulk-first');
    grok.shell.addTableView(first);
    await waitUntil(() => selectOf(b.root, 'Table').options.length > 0, 'the Translators were not offered the table', 5000);

    choose(selectOf(a.root, 'Output format'), DEFAULT_FORMATS.HELM);
    convertIn(b);
    await waitUntil(() => first.col('AxolabsSequences (Nucleotides)') !== null,
      `the second Translator, showing Nucleotides, did not write them: ${first.columns.names().join(', ')}`, 5000);
    expect(first.col('AxolabsSequences (HELM)'), null);
    convertIn(a);
    await waitUntil(() => first.col('AxolabsSequences (HELM)') !== null,
      `the first Translator did not write HELM: ${first.columns.names().join(', ')}`, 5000);

    const second = axolabsTable('bulk-second');
    grok.shell.addTableView(second);   // both Translators take the table just opened
    await waitUntil(() => Array.from(selectOf(a.root, 'Table').options).some((o) => o.text === 'bulk-second'),
      'the first Translator was not offered the second table', 5000);
    choose(selectOf(a.root, 'Table'), 'bulk-first');
    convertIn(b);
    await waitUntil(() => second.columns.length > 1,
      `the second Translator did not convert its own table: ${second.columns.names().join(', ')}`, 5000);
    expect(first.columns.length, 3);
  });

  // A Translator's table controls followed the workspace for ever: a closed one still took every table opened later.
  test('A closed Translator no longer follows the tables opened after it', async () => {
    const view = await openTranslator();
    const root = view.root;
    view.close();
    grok.shell.addTableView(axolabsTable('bulk-after-close'));
    await new Promise((resolve) => setTimeout(resolve, 1000));
    const offered = Array.from(selectOf(root, 'Table').options).map((o) => o.text);
    expect(offered.includes('bulk-after-close'), false);
  });

  // In the Oligo Toolkit the Translator is a tab of the toolkit's view: the shell removes that one, not the
  // Translator's own.
  test('A closed Oligo Toolkit\'s Translator no longer follows the tables opened after it', async () => {
    const view: DG.ViewBase = await grok.functions.call('SequenceTranslator:oligoToolkitApp');
    grok.shell.addView(view);
    await waitUntil(() => $(view.root).find('textarea').length > 0, 'the toolkit\'s Translator did not appear', 10000);
    const root = view.root;
    view.close();
    grok.shell.addTableView(axolabsTable('bulk-after-toolkit-close'));
    await new Promise((resolve) => setTimeout(resolve, 1000));
    const offered = Array.from(selectOf(root, 'Table').options).map((o) => o.text);
    expect(offered.includes('bulk-after-toolkit-close'), false);
  });
});
