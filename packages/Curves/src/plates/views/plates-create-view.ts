import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {PlateTemplate, plateTemplates, plateTypes, savePlate} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';
  const platePropertiesHost = ui.divV([]);
  let plateTemplateValues = {};
  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];

  const updatePlate = () => {
    plate = new Plate(plateType.rows, plateType.cols);
    for (const property of plateTemplate.wellProperties) {
      plate.data.columns.addNew(property.name!, property.value_type! as DG.ColumnType);
    }
    plateWidget.plateData = plate.data;
    plateWidget.editable = true;
  }

  const setTemplate = (template: PlateTemplate) => {
    plateTemplate = template;
    plateTemplateValues = {};
    const form = ui.input.form(plateTemplateValues, template.plateProperties.map(p => DG.Property.js(p.name!, p.value_type! as DG.TYPE)));
    ui.empty(platePropertiesHost);
    platePropertiesHost.appendChild(form);
    updatePlate();
  }

  const plateTypeSelector = ui.input.choice('Plate type', {
      items: plateTypes.map(pt => pt.name),
      value: plateTypes[0].name,
      onValueChanged: (v) => {
          plateType = plateTypes.find(pt => pt.name === v)!;
          updatePlate();
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

  view.setRibbonPanels([
    [ui.bigButton('CREATE', async() => {
      const plate = plateWidget.plate;
      plate.details = plateTemplateValues;
      console.log(plateWidget.plate.data.toCsv());
      await savePlate(plate);
      grok.shell.info(`Plate created: ${plate.id}`);
    })]
  ]);

  return view;
}
