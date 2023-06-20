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
  // let view: DG.TableView;
  // let df: DG.DataFrame;

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
        const senseStrandInputs: NodeListOf<HTMLSelectElement> = document.querySelectorAll('.st-pattern-choice-input > select');
        await delay(1000);
        senseStrandInputs[1].value = 'DNA';
        console.log('value changed:', senseStrandInputs[1].value);
        await delay(1000);
        senseStrandInputs[2].value = 'invAb';
        console.log('value changed:', senseStrandInputs[2].value);
      }, {
        description: `Create a modification pattern for a dimer`,
        delay: 2000,
      })
      .start();
  } catch (err: any) {
    handleError(err);
  }
}
