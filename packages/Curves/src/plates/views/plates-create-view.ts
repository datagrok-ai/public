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

  const setTemplate = (template: PlateTemplate) => {
    plateTemplateValues = {};
    const form = ui.input.form(plateTemplateValues, template.plateProperties.map(p => DG.Property.js(p.name!, p.value_type! as DG.TYPE)));
    ui.empty(platePropertiesHost);
    platePropertiesHost.appendChild(form);
  }

  const plateTypeSelector = ui.input.choice('Plate type', {
      items: plateTypes.map(pt => pt.name),
      value: plateTypes[0].name,
      onValueChanged: (v) => {
          const selectedType = plateTypes.find(pt => pt.name === v)!;
          plate = new Plate(selectedType.rows, selectedType.cols);
          plate.data.columns.addNew('activity', DG.TYPE.FLOAT).init(i => Math.random() * 100);
          plateWidget.plateData = plate.data;
      }
  });

  const plateTemplateSelector = ui.input.choice('Template', {
    items: plateTemplates.map(pt => pt.name),
    value: plateTemplates[0].name,
    onValueChanged: (v) => setTemplate(plateTemplates.find(pt => pt.name === v)!)
  });
    
  // Create plate viewer
  let plate = new Plate(plateTypes[0].rows, plateTypes[0].cols);
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