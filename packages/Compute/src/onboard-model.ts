import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {Property, PropertyOptions} from "datagrok-api/dg";

/* eslint-disable */

interface IModelTemplate {
  name: string;
  props: PropertyOptions[];
}

export async function onboardModel(): Promise<void> {

  let modelPropertiesHost = ui.divV([]);
  let name = ui.stringInput('Name', '');
  let type = ui.choiceInput('Type', templates[0].name, templates.map((t) => t.name));
  let modelProperties: any = {};

  function updateProps() {
    let template: IModelTemplate = templates.filter((t) => t.name == type.value)[0];
    ui.empty(modelPropertiesHost);
    modelPropertiesHost.appendChild(ui.input.form(modelProperties, template.props.map((p) => Property.fromOptions(p))));
  }

  updateProps();
  type.onChanged(updateProps);

  new DG.Wizard({ title: 'Onboard a model' })
    .page({
      caption: 'Details',
      root: ui.divV([
        name,
        type,
        modelPropertiesHost
      ]),
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

let templates: IModelTemplate[] = [
  {
    "name": "Bioreactor simulation",
    "props": [
      {
        "name": "reactorType",
        "type": "string",
        "choices": ["Experimental", "Production"]
      },
      {
        "name": "department",
        "type": "string",
        "choices": ["", "Production"]
      }
    ]
  },
  {
    "name": "Biosensors",
    "props": [
      {
        "name": "sensor",
        "type": "string",
        "choices": ["Accelerometer", "GSR", "EEG", "ECG", "Eye Tracker"]
      }
    ]
  }
];