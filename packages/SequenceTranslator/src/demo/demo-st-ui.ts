import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PackageFunctions} from '../package';
import {tryCatch} from '../apps/common/model/helpers';

export async function demoOligoTranslatorUI() {
  await tryCatch(async () => {
    const view = await PackageFunctions.oligoTranslatorApp();
    grok.shell.addView(view);
  });
}

export async function demoOligoPatternUI() {
  await tryCatch(async () => {
    const view = await PackageFunctions.oligoPatternApp();
    grok.shell.addView(view);
  });
}

export async function demoOligoStructureUI() {
  await tryCatch(async () => {
    async function setInputValue(idx: number, sequence: string): Promise<void> {
      // the app builds its inputs after its view is added
      const find = () => document.querySelectorAll<HTMLTextAreaElement>('.st-colored-text-input > textarea')[idx];
      const textarea = find() ?? await new Promise<HTMLTextAreaElement>((resolve) => {
        const observer = new MutationObserver(() => {
          const found = find();
          if (found) {
            observer.disconnect();
            resolve(found);
          }
        });
        observer.observe(document.body, {childList: true, subtree: true});
      });
      textarea.value = sequence;
      const event = new Event('input');
      textarea.dispatchEvent(event);
    }

    const view = await PackageFunctions.oligoStructureApp();
    grok.shell.addView(view);
    const inputSequences = ['Afcgacsu', 'Afcgacsu', 'Afcgacsu'];
    inputSequences.forEach(async (sequence, idx) => {
      await setInputValue(idx, sequence);
    });
  });
}
