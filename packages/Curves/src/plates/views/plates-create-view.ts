import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {
  getPlateById,
  initPlates,
  PlateTemplate,
  plateTemplates,
  PlateType,
  plateTypes,
  savePlate,
  savePlateAsTemplate
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';

/** Creates a new plate for a user to edit. Does not add it to the database. */
async function createNewPlateForTemplate(plateType: PlateType, plateTemplate: PlateTemplate){
  if (!plateTemplate.plate_layout_id) {
    const plate = new Plate(plateType.rows, plateType.cols);
    for (const property of plateTemplate.wellProperties) {
      plate.data.columns.addNew(property.name!, property.value_type! as DG.ColumnType);
    }
    return plate;
  }

  const p = await getPlateById(plateTemplate.plate_layout_id);
  p.id = undefined;
  return p;
}

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';
  const platePropertiesHost = ui.divV([]);
  let plateTemplateValues = {};
  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];

  const updatePlate = async () => {
    plateWidget.plate = await createNewPlateForTemplate(plateType, plateTemplate);
    plateWidget.editable = true;
  }

  const setTemplate = (template: PlateTemplate) => {
    plateTemplate = template;
    plateTemplateValues = {};
    const form = ui.input.form(plateTemplateValues, template.plateProperties.map(p => DG.Property.js(p.name!, p.value_type! as DG.TYPE)));
    ui.empty(platePropertiesHost);
    platePropertiesHost.appendChild(form);
    updatePlate().then(_ => {});
  }

  const plateTypeSelector = ui.input.choice('Plate type', {
      items: plateTypes.map(pt => pt.name),
      value: plateTypes[0].name,
      onValueChanged: (v) => {
        plateType = plateTypes.find(pt => pt.name === v)!;
        updatePlate().then(_ => {});
      }
  });

  const plateTemplateSelector = ui.input.choice('Template', {
    items: plateTemplates.map(pt => pt.name),
    value: plateTemplates[0].name,
    onValueChanged: (v) => setTemplate(plateTemplates.find(pt => pt.name === v)!)
  });

  // Create plate viewer
  let plate = new Plate(plateType.rows, plateType.cols);
  const plateWidget = PlateWidget.fromPlate(plate);
  plateWidget.editable = true;
  plateWidget.root.style.height = '400px';
  setTemplate(plateTemplates[0]);

  // Add components to view
  view.root.appendChild(ui.divV([
      ui.form([
        plateTypeSelector,
        plateTemplateSelector
      ]),
      platePropertiesHost,
      plateWidget.root
  ]));

  const getPlate = () => {
    const plate = plateWidget.plate;
    plate.details = plateTemplateValues;
    return plate;
  }

  view.setRibbonPanels([[
    ui.bigButton('CREATE', async() => {
      const plate = getPlate();
      await savePlate(plate);
      grok.shell.info(`Plate created: ${plate.id}`);
    }),

    ui.button('SAVE TEMPLATE', async() => {
      const plate = getPlate();
      await savePlateAsTemplate(plate, plateTemplate);
      await initPlates();
      grok.shell.info(`Plate template saved: ${plate.id}`);
    })
  ]]);

  return view;
}
