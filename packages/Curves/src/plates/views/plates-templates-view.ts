import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {plateTemplates} from '../plates-crud';
import {propertySchemaView} from './plates-schema-view';
import {PlateTemplateHandler} from '../objects/plate-template-handler';

export function createTemplatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Templates';

  const templatesTable = ui.table(plateTemplates, (t) => [
    t.name,
    t.description,
    ui.button('Edit', () => PlateTemplateHandler.editTemplate(t)),
    ui.button('Clone', () => PlateTemplateHandler.openCloneTemplateView(t)),
  ]);

  view.root.appendChild(templatesTable);

  view.setRibbonPanels([[
    ui.bigButton('CREATE NEW', async () => {
      const newTemplate = {
        name: 'New Template',
        description: '',
        plateProperties: [],
        wellProperties: [
          {name: 'Concentration', type: DG.COLUMN_TYPE.FLOAT},
          {name: 'Volume', type: DG.COLUMN_TYPE.FLOAT},
          {name: 'Sample', type: DG.COLUMN_TYPE.STRING},
        ],
        id: -1,
        plate_layout_id: -1
      };

      grok.shell.addPreview(propertySchemaView(newTemplate));
    }),
  ]]);

  return view;
}
