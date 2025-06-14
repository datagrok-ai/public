import { PlateProperty, PlateTemplate, plateTypes, savePlate } from "./plates-crud";

import { initPlates } from "./plates-crud";

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from "../plate/plate";
import { createPlateTemplate, createNewPlateForTemplate } from "./plates-crud";


export async function __createDummyPlateData() {
  await initPlates(true);

  const cellCountingTemplate = await createPlateTemplate({
    name: 'Cell counting',
    description: 'Microscopy-based cell counting',
    plateProperties: [
      {name: 'Imaging device', choices: ['Kodak', 'Nikon'], value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Status', choices: ['Pending', 'Filling', 'Measuring', 'Done'], value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Plate cell count', min: 0, max: 10000, value_type: DG.COLUMN_TYPE.INT}
    ],
    wellProperties: [
      {name: 'Well cell count', min: 0, max: 100, value_type: DG.COLUMN_TYPE.INT},
      {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], value_type: DG.COLUMN_TYPE.STRING}
    ]
  });

  const doseResponseTemplate = await createPlateTemplate({
    name: 'Dose-response',
    description: 'Dose-response campaign',
    plateProperties: [
      {name: 'Project', value_type: DG.COLUMN_TYPE.STRING, choices: ['Modulators of mGluR5', 'Agonists for GPCR GPR139', 'Glutaminase Inhibitors for TNBC']},
      {name: 'Stage', value_type: DG.COLUMN_TYPE.STRING, choices: ['Lead generation', 'Lead optimization']},
      {name: 'Chemist', value_type: DG.COLUMN_TYPE.STRING, choices: ['John Marlowski', 'Mary Hopton']},
      {name: 'Biologist', value_type: DG.COLUMN_TYPE.STRING, choices: ['Anna Fei', 'Joan Dvorak']},
      {name: 'QC Passed', value_type: DG.COLUMN_TYPE.BOOL},
      {name: 'Z-Score', min: 0, max: 3, value_type: DG.COLUMN_TYPE.FLOAT},
    ],
    wellProperties: [
      {name: 'Sample', choices: ['GRK-1', 'GRK-2', 'GRK-3', 'GRK-4', 'GRK-5', 'GRK-6'], value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Role', choices: ['Control', 'Treatment'], value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Concentration', min: 0, max: 100, value_type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Volume', min: 0, max: 100, value_type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Activity', min: 0, max: 100, value_type: DG.COLUMN_TYPE.FLOAT},
    ]
  });

  await initPlates(true);

  const getDemoValue = (property: PlateProperty) => 
    property.choices ? DG.Utils.random(property.choices) :
    property.value_type === DG.COLUMN_TYPE.BOOL ? Math.random() > 0.5 :
    property.min !== undefined && property.max !== undefined ? 
      property.value_type === DG.COLUMN_TYPE.INT ?
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

  await initPlates(true);
}