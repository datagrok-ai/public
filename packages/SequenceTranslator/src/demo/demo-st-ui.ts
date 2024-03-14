import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';
import {loadJsonData} from '../apps/common/data-loader/json-loader';
import {_package, oligoTranslatorApp, oligoPatternApp, oligoStructureApp} from '../package';
import {tryCatch} from '../apps/common/model/helpers';

export async function demoOligoTranslatorUI() {
  await tryCatch(async () => oligoTranslatorApp());
}

export async function demoOligoPatternUI() {
  await tryCatch(async () => {
    async function emulateUserInput(value: string, idx: number, idxUpdate: (idx: number) => number): Promise<void> {
      await delay(3000);

      // warning: this redefinition is necessary because
      // the ids of the elements can dynamically change
      const choiceInputs: NodeListOf<HTMLSelectElement> = document.querySelectorAll('.st-pattern-choice-input > select');
      len = choiceInputs.length;
      const selectElement = choiceInputs[idxUpdate(idx)];
      selectElement.value = value;
      const event = new Event('input');
      selectElement.dispatchEvent(event);
    }

    await oligoPatternApp();

    let len: number;

    const ssNewValues = ['DNA', 'invAb', 'Z-New'];
    ssNewValues.forEach(async (value, idx) => {
      emulateUserInput(value, idx, (i) => 2 * i);
    });

    const asNewValues = ['2\'-O-Methyl', '2\'-Fluoro', '2\'-O-MOE'];
    asNewValues.forEach(async (value, idx) => {
      emulateUserInput(value, idx, (i) => (len - 2 - 2 * i));
    });
  });
}

export async function demoOligoStructureUI() {
  await tryCatch(async () => {
    async function setInputValue(idx: number, sequence: string): Promise<void> {
      await delay(500);
      const textInputs: NodeListOf<HTMLTextAreaElement> = document.querySelectorAll('.colored-text-input > textarea');
      const textarea = textInputs[idx];
      textarea.value = sequence;
      const event = new Event('input');
      textarea.dispatchEvent(event);
    }
    await oligoStructureApp();
    const inputSequences = ['Afcgacsu', 'Afcgacsu', 'Afcgacsu'];
    inputSequences.forEach(async (sequence, idx) => {
      await setInputValue(idx, sequence);
    });
  });
}
