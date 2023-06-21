import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {handleError} from './handle-error';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';
import {delay} from '@datagrok-libraries/utils/src/test';
import {getJsonData} from '../model/data-loading-utils/json-loader';
import {SequenceTranslatorUI} from '../view/view';
import {_package} from '../package';

export async function demoSequenceTranslatorUI() {
  try {
    const demoScript = new DemoScript(
      'Sequence Design',
      'Sequence Translator is an application for design and visualization of oligonucleotide sequences', undefined, {autoStartFirstStep: true});

    let tabControl: DG.TabControl;
    let panes: DG.TabPane[];

    await demoScript
      .step(`Translate`, async () => {
        await getJsonData();
        await _package.initMonomerLib();
        const v = new SequenceTranslatorUI();
        await v.createLayout();
        tabControl = (await v.tabs.getControl());
        panes = tabControl.panes;
        tabControl.currentPane = panes[0];
      }, {
        description: `Translate sequences across various formats aceepted by synthesizers`,
        delay: 2000,
      })
      .step(`Create pattern`, async () => {
        tabControl.currentPane = panes[1];
        const ssNewValues = ['DNA', 'invAb', 'Z-New'];

        let len: number;

        async function emulateUserInput(value: string, idx: number, idxUpdate: (idx: number) => number) {
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

        ssNewValues.forEach(async (value, idx) => {
          emulateUserInput(value, idx, (i) => 2 * i);
        });

        const asNewValues = ['2\'-O-Methyl', '2\'-Fluoro', '2\'-O-MOE'];
        asNewValues.forEach(async (value, idx) => {
          emulateUserInput(value, idx, (i) => (len - 2 - 2 * i));
        })
      }, {
        description: `Create a modification pattern for a dimer`,
        delay: 0,
      })
      .step(`Visualize Duplex`, async () => {
        tabControl.currentPane = panes[2];
      }, {
        description: `Get atomic-level structure for the duplex`,
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
