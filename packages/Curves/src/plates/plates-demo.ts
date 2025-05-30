import { savePlate } from "./plates-crud";

import { initPlates } from "./plates-crud";

import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {Plate} from "../plate/plate";
import { createPlateTemplate } from "./plates-crud";


export async function __createDummyPlateData() {
  await initPlates();

  await createPlateTemplate({
    name: 'Cell counting',
    description: 'Microscopy-based cell counting',
    plateProperties: [
      {name: 'Imaging device', value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Status', value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Plate cell count', value_type: DG.COLUMN_TYPE.INT}
    ],
    wellProperties: [
      {name: 'Well cell count', value_type: DG.COLUMN_TYPE.INT},
      {name: 'Sample', value_type: DG.COLUMN_TYPE.STRING}
    ]
  });

  await createPlateTemplate({
    name: 'Dose-response',
    description: 'Dose-response campaign',
    plateProperties: [
      {name: 'QC Passed', value_type: DG.COLUMN_TYPE.BOOL},
    ],
    wellProperties: [
      {name: 'Sample', value_type: DG.COLUMN_TYPE.STRING},
      {name: 'Concentration', value_type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Volume', value_type: DG.COLUMN_TYPE.FLOAT},
      {name: 'Activity', value_type: DG.COLUMN_TYPE.FLOAT},
    ]
  });


  for (let i = 0; i < 50; i++) {
    const plate = Plate.demo();
    plate.details = {
      'Project': DG.Utils.random([
        'Allosteric Modulators of mGluR5',
        'Agonists for Orphan GPCR GPR139',
        'Glutaminase Inhibitors for TNBC']),
      'Stage': DG.Utils.random(['Lead generation', 'Lead optimization']),
      'Chemist': DG.Utils.random(['John Marlowski', 'Mary Hopton']),
      'Passed QC': DG.Utils.random([true, false]),
      'Z-score': Math.random() * 3,
    }
    await savePlate(plate);
  }

  for (let i = 0; i < 50; i++) {
    const plate = Plate.demo();
    plate.details = {
      'Project': DG.Utils.random([
        'NaV1.7 Blockers for Pain Relief',
        'TLR4 Antagonists for Sepsis',]),
      'Stage': DG.Utils.random(['Stage 1', 'Stage 2', 'Stage 3']),
      'Chemist': DG.Utils.random(['John Marlowski', 'Andrew Smith']),
      'Biologist': DG.Utils.random(['Anna Fei', 'Joan Dvorak']),
      'Cells': Math.ceil(Math.random() * 100)
    }
    await savePlate(plate);
  }

  grok.shell.info('100 plates saved');
}
