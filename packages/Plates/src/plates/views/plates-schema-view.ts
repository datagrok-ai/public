/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {createPlateTemplate, initPlates, PlateTemplate} from '../plates-crud';
import {SchemaEditor} from '@datagrok-libraries/utils/src/schema-editor';
import {merge} from 'rxjs';


export function propertySchemaView(template: PlateTemplate): DG.View {
  const view = DG.View.create();
  view.name = 'Templates / ' + template.name;

  const extraPropertiesDiv = ui.div([]);
  const platePropEditor =new SchemaEditor({properties: template.plateProperties, extraPropertiesDiv: extraPropertiesDiv});
  const wellPropEditor =new SchemaEditor({properties: template.wellProperties, extraPropertiesDiv: extraPropertiesDiv});
  const nameEditor = ui.input.string('Name', {value: template.name});
  const descriptionEditor = ui.input.string('Description', {value: template.description});

  const saveButton = ui.bigButton('SAVE', async () => {
    template.name = nameEditor.value;
    template.description = descriptionEditor.value;
    template.plateProperties = platePropEditor.properties;
    template.wellProperties = wellPropEditor.properties;
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

    ui.divH([
      ui.divV([
        ui.h2('Plate properties'),
        platePropEditor,
        ui.h2('Well properties'),
        wellPropEditor,
      ]),
      ui.div([], {style: {width: '20px'}}),
      extraPropertiesDiv
    ]),
  ]));

  view.setRibbonPanels([[saveButton]]);

  return view;
}
