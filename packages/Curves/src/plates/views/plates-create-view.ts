/* eslint-disable max-len */
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {Plate} from '../../plate/plate';
import {
  PlateTemplate,
  plateTemplates,
  PlateType,
  plateTypes,
  savePlate,
  savePlateAsTemplate,
  createNewPlateForTemplate,
  initPlates
} from '../plates-crud';
import {PlateWidget} from '../../plate/plate-widget';
import {parsePlateFromCsv} from '../../plate/csv-plates';
import {renderValidationResults} from './plates-validation-panel';

export function createPlatesView(): DG.View {
  const view = DG.View.create();
  view.name = 'Create Plate';

  const templateState = new Map<number, {plate: Plate, file: DG.FileInfo}>();

  const platePropertiesHost = ui.divV([]);
  const validationHost = ui.divV([]);
  const importHost = ui.div();

  let plateTemplateValues = {};
  let plateType = plateTypes[0];
  let plateTemplate = plateTemplates[0];
  if (!plateTemplate) {
    grok.shell.warning('No plate templates found. Please create a template first.');
    throw new Error('No plate templates found. Please create a template first.');
  }

  const plate = new Plate(plateType.rows, plateType.cols);
  const plateWidget = PlateWidget.fromPlate(plate);
  plateWidget.editable = true;
  plateWidget.root.style.height = '400px';

  let csvFileInput: DG.InputBase<any>;
  let excelFileInput: DG.InputBase<any>;

  const setTemplate = async (template: PlateTemplate) => {
    plateTemplate = template;
    ui.empty(validationHost);
    ui.empty(platePropertiesHost);

    const state = templateState.get(template.id);

    if (state) {
      plateTemplateValues = state.plate.details;
      plateWidget.plate = state.plate;
      validationHost.appendChild(renderValidationResults(state.plate, template));
      updateImportHost(state.file);
    } else {
      plateTemplateValues = {};
      plateWidget.plate = await createNewPlateForTemplate(plateType, template);
      updateImportHost(null);
    }

    const form = ui.input.form(plateTemplateValues, template.plateProperties.map((p) => DG.Property.js(p.name!, p.type! as DG.TYPE)));
    platePropertiesHost.appendChild(form);
    plateWidget.refresh();
  };

  const updateImportHost = (file: DG.FileInfo | null) => {
    ui.empty(importHost);
    if (file) {
      const clearIcon = ui.iconFA('times', () => {
        templateState.delete(plateTemplate.id);
        setTemplate(plateTemplate);
      }, 'Clear imported data');
      clearIcon.style.marginLeft = '10px';
      clearIcon.style.cursor = 'pointer';

      const fileLabel = ui.divText(file.name);
      fileLabel.style.overflow = 'hidden';
      fileLabel.style.textOverflow = 'ellipsis';
      fileLabel.style.whiteSpace = 'nowrap';
      fileLabel.style.maxWidth = '200px';

      importHost.appendChild(ui.divH([
        ui.label('Import from CSV'),
        fileLabel,
        clearIcon
      ], {style: {alignItems: 'center'}}));
    } else {
      excelFileInput = ui.input.file('Import from Excel', {});
      csvFileInput = ui.input.file('Import from CSV', {
        onValueChanged: async (file) => {
          const plate = await parsePlateFromCsv(await file.readAsString());

          console.log(`Parsed Plate Data from ${file.name}:`);
          console.table(plate.data.toJson());

          templateState.set(plateTemplate.id, {plate, file});
          await setTemplate(plateTemplate);
        }
      });
      importHost.appendChild(ui.divH([excelFileInput.root, csvFileInput.root], {style: {'gap': '10px'}}));
    }
  };

  const plateTypeSelector = ui.input.choice('Plate type', {
    items: plateTypes.map((pt) => pt.name),
    value: plateTypes[0].name,
    onValueChanged: (v) => {
      plateType = plateTypes.find((pt) => pt.name === v)!;
      setTemplate(plateTemplate);
    }
  });

  const plateTemplateSelector = ui.input.choice('Template', {
    items: plateTemplates.map((pt) => pt.name),
    value: plateTemplates[0]?.name,
    onValueChanged: (v) => setTemplate(plateTemplates.find((pt) => pt.name === v)!)
  });

  setTemplate(plateTemplates[0]);

  view.root.appendChild(ui.divV([
    ui.form([
      plateTypeSelector,
      plateTemplateSelector,
    ]),
    importHost,
    validationHost,
    platePropertiesHost,
    plateWidget.root
  ]));

  const getPlate = () => {
    const plate = plateWidget.plate;
    plate.details = plateTemplateValues;
    return plate;
  };

  view.setRibbonPanels([[
    ui.bigButton('CREATE', async () => {
      const plate = getPlate();
      await savePlate(plate);
      grok.shell.info(`Plate created: ${plate.id}`);
    }),
    ui.button('SAVE TEMPLATE', async () => {
      const plate = getPlate();
      await savePlateAsTemplate(plate, plateTemplate);
      await initPlates(true);
      grok.shell.info(`Plate template saved: ${plate.id}`);
    })
  ]]);

  return view;
}
