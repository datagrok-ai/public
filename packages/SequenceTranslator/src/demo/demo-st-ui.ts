import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {delay} from '@datagrok-libraries/utils/src/test';
import {getJsonData} from '../model/data-loading-utils/json-loader';
import {_package, oligoTranslatorApp} from '../package';
import {tryCatch} from '../model/helpers';

export async function demoOligoTranslatorUI() {
  await tryCatch(async () => oligoTranslatorApp());
}

// export async function demoDesignPatternUI() {
//   try {
//     async function emulateUserInput(value: string, idx: number, idxUpdate: (idx: number) => number) {
//       await delay(3000);

//       // warning: this redefinition is necessary because
//       // the ids of the elements can dynamically change
//       const choiceInputs: NodeListOf<HTMLSelectElement> = document.querySelectorAll('.st-pattern-choice-input > select');
//       len = choiceInputs.length;
//       const selectElement = choiceInputs[idxUpdate(idx)];
//       selectElement.value = value;
//       const event = new Event('input');
//       selectElement.dispatchEvent(event);
//     }

//     openSequenceTranslatorOnPane(1);

//     let len: number;

//     const ssNewValues = ['DNA', 'invAb', 'Z-New'];
//     ssNewValues.forEach(async (value, idx) => {
//       emulateUserInput(value, idx, (i) => 2 * i);
//     });

//     const asNewValues = ['2\'-O-Methyl', '2\'-Fluoro', '2\'-O-MOE'];
//     asNewValues.forEach(async (value, idx) => {
//       emulateUserInput(value, idx, (i) => (len - 2 - 2 * i));
//     })
//   } catch (err: any) {
//     handleError(err);
//   }
// }

// export async function demoVisualizeDuplexUI() {
//   try {
//     await openSequenceTranslatorOnPane(2);
//   } catch (err: any) {
//     handleError(err);
//   }
// }

// async function openSequenceTranslatorOnPane(paneNumber: number): Promise<void> {
//   let tabControl: DG.TabControl;
//   let panes: DG.TabPane[];
//   await getJsonData();
//   await _package.initMonomerLib();
//   const v = new UnifiedUI();
//   tabControl = (await v.tabs.getControl());
//   panes = tabControl.panes;
//   tabControl.currentPane = panes[paneNumber];
// }
