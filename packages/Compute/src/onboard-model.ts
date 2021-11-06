import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';

export function onboardModel() {
  let name = ui.stringInput('Name', '');
  let type = ui.choiceInput('Type', 'Bio', ['Bio', 'Chem']);

  new DG.Wizard({ title: 'Onboard a model' })
    .page({
      caption: 'Details',
      root: ui.inputs([name, type]),
      validate: () => name.value !== '' ? null : 'Please provide a name'
    })
    .page({
      caption: 'Model Source',
      root: ui.divText('Second page') })
    .page({
      caption: 'Finalize',
      root: ui.divText('Finish') })
    .onOK(() => {})
    .show({width: 500, height: 300});
}