import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {createPlateTemplate, initPlates, PlateProperty, PlateTemplate} from '../plates-crud';
import {SchemaEditor} from '@datagrok-libraries/utils/src/schema-editor';
import { merge } from 'rxjs';


function platePropertyToOptions(p: Partial<PlateProperty>): DG.IProperty {
  return { name: p.name, type: p.value_type}
}

export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = 'Templates / ' + template.name;

  const platePropEditor = new SchemaEditor({properties: template.plateProperties.map(platePropertyToOptions)});
  const wellPropEditor = new SchemaEditor({properties: template.wellProperties.map(platePropertyToOptions)});
  const nameEditor = ui.input.string('Name', {value: template.name});
  const descriptionEditor = ui.input.string('Description', {value: template.description});

  const saveButton = ui.bigButton('SAVE', async() => {
    template.name = nameEditor.value;
    template.description = descriptionEditor.value;
    template.plateProperties = platePropEditor.properties.map(p => ({...p, value_type: p.type}));
    template.wellProperties = wellPropEditor.properties.map(p => ({...p, value_type: p.type}));
    await createPlateTemplate(template);
    await initPlates(true);
    grok.shell.info('Template saved');
  });
 
  const updateSaveVisibility = () => saveButton.style.display = template.id !== -1 ? 'none' : 'block';
  updateSaveVisibility();

  merge(platePropEditor.table.onChanged, wellPropEditor.table.onChanged).subscribe((_) => updateSaveVisibility());

  view.root.appendChild(ui.divV([
    ui.h2('Template'),
    ui.form([
      nameEditor,
      descriptionEditor,
    ]),

    ui.h2('Plate properties'),
    platePropEditor,

    ui.h2('Well properties'),
    wellPropEditor,
  ]));

  view.setRibbonPanels([[saveButton]]);

  return view;
}