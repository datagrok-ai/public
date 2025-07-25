import { PlateProperty, plateTypes, savePlate } from "./plates-crud";
import { initPlates } from "./plates-crud";
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import { createPlateTemplate, createNewPlateForTemplate } from "./plates-crud";
import { Plate } from '../plate/plate';


export async function __createDummyPlateData() {
  await initPlates(true);

  await createDummyPlates();
  await createDummyPlatesFromExcel();

  await initPlates(true);
}

/** Demonstrates construction of templates and plates on the fly */
async function createDummyPlates() {
  const cellCountingTemplate = await createPlateTemplate({
    name: 'Cell counting',
    description: 'Microscopy-based cell counting',
    plateProperties: [
      {name: 'Imaging device', choices: ['Kodak', 'Nikon'], type: DG.COLUMN_TYPE.STRING},
      {name: 'Status', choices: ['Pending', 'Filling', 'Measuring', 'Done'], type: DG.COLUMN_TYPE.STRING},
      {name: 'Plate cell count', min: 0, max: 10000, type: DG.COLUMN_TYPE.INT}
    ],
    wellProperties: [
      {name: 'Well cell count', min: 0, max: 100, type: DG.COLUMN_TYPE.INT},
      {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], type: DG.COLUMN_TYPE.STRING}
    ]
  });

  const doseResponseTemplate = await createPlateTemplate({
    name: 'Dose-response',
    description: 'Dose-response campaign',
    plateProperties: [
      {name: 'Project', type: DG.COLUMN_TYPE.STRING, choices: ['Modulators of mGluR5', 'Agonists for GPCR GPR139', 'Glutaminase Inhibitors for TNBC']},
      {name: 'Stage', type: DG.COLUMN_TYPE.STRING, choices: ['Lead generation', 'Lead optimization']},
      {name: 'Chemist', type: DG.COLUMN_TYPE.STRING, choices: ['John Marlowski', 'Mary Hopton']},
      {name: 'Biologist', type: DG.COLUMN_TYPE.STRING, choices: ['Anna Fei', 'Joan Dvorak']},
      {name: 'QC Passed', type: DG.COLUMN_TYPE.BOOL},
      {name: 'Z-Score', min: 0, max: 3, type: DG.COLUMN_TYPE.FLOAT},
    ],
    wellProperties: [
      {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], type: DG.COLUMN_TYPE.STRING},
      {name: 'Role', choices: ['Control', 'Treatment'], type: DG.COLUMN_TYPE.STRING},
      {name: 'Concentration', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Volume', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Activity', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT},
    ]
  });

  await initPlates(true);

  const getDemoValue = (property: PlateProperty) =>
    property.choices ? DG.Utils.random(property.choices) :
    property.type === DG.COLUMN_TYPE.BOOL ? Math.random() > 0.5 :
    property.min !== undefined && property.max !== undefined ?
      property.type === DG.COLUMN_TYPE.INT ?
        Math.floor(property.min + Math.random() * (property.max - property.min)) :
        property.min + Math.random() * (property.max - property.min) :
    null;

  for (const template of [cellCountingTemplate, doseResponseTemplate]) {
    for (let i = 0; i < 20; i++) {
      const plate = await createNewPlateForTemplate(plateTypes[0], template);
      plate.details = Object.fromEntries(
        template.plateProperties.map(p => [p.name!, getDemoValue(p as PlateProperty)])
      );

      // Initialize well properties
      for (const property of template.wellProperties)
        plate.data.col(property!.name!)?.init((_) => getDemoValue(property as PlateProperty));

      await savePlate(plate);
    }
  }
}

/**
 * Demonstrates importing plates from Excel files
 * and forcing them to comply to the template
*/
async function createDummyPlatesFromExcel() {
  // create a template for excel plates (you only need to do this once)
  const excelPlateTemplate = await createPlateTemplate({
    name: 'Excel',
    description: 'Excel plates',
    plateProperties: [ {name: 'Barcode', type: DG.COLUMN_TYPE.STRING } ],
    wellProperties: [
      {name: 'raw data', type: DG.COLUMN_TYPE.FLOAT },
      {name: 'plate layout', type: DG.COLUMN_TYPE.STRING},
      {name: 'concentrations', min: 0, max: 100, type: DG.COLUMN_TYPE.FLOAT}
    ]
  });
  await initPlates(true);

  // import excel plates from the directory
  const xlsxFiles = await grok.dapi.files.list('System.DemoFiles/hts/xlsx_plates', false, '*.xlsx');
  for (const file of xlsxFiles) {
    const plate = await Plate.fromExcelFile(file);
    plate.plateTemplateId = excelPlateTemplate.id;
    plate.details = { Barcode: (Math.round(Math.random() * 100000)).toString().padStart(8, '0') };
    await savePlate(plate, {autoCreateProperties: false});
  }
}
