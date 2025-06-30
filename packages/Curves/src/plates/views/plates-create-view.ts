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
  savePlateAsTemplate,
  createNewPlateForTemplate
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';


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
    const form = ui.input.form(plateTemplateValues, template.plateProperties.map(p => DG.Property.js(p.name!, p.type! as DG.TYPE)));
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

  const fileInput = ui.input.file('Import from Excel', {
    nullable: true,
    onValueChanged: async (file) => {
      const plate = await Plate.fromExcel(await file.readAsBytes());
      if (plate.rows !== plateType.rows || plate.cols !== plateType.cols) {
        grok.shell.error(`Plate dimensions do not match the template: ${plate.rows}x${plate.cols} vs ${plateType.rows}x${plateType.cols}`);
        return;
      }

      for (const layer of plate.getLayerNames().filter(l => plateTemplate.wellProperties.find(p => p.name?.toLowerCase() === l.toLowerCase()))) {
        const excelCol = plate.data.col(layer)!;
        plateWidget.plate.data.col(layer)!.init(i => excelCol.get(i));
      }

      const missingLayers = plateTemplate.wellProperties.filter(p => !plate.getLayerNames().includes(p.name!));
      if (missingLayers.length > 0) 
        grok.shell.warning(`Plate is missing the following layers: ${missingLayers.map(p => p.name!).join(', ')}`);

      plateWidget.refresh();
    }
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
        plateTemplateSelector,
        fileInput
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
