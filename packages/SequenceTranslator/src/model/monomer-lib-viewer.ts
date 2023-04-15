/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {drawZoomedInMolecule} from '../view/utils/draw-molecule';

import {LIB_PATH, DEFAULT_LIB_FILENAME} from './data-loader/const';
import {SYNTHESIZERS, TECHNOLOGIES} from './const';

import {Monomer} from '@datagrok-libraries/bio/src/types';

type CodesField = {
  [synthesizer: string]: {
    [technology: string]: string[]
  }
}

type MonomerExtension = {
  [key: string]: string | CodesField
}

type ExtendedMonomer = Monomer & MonomerExtension;

type FormattedMonomer = {
  [key: string]: string
}

const enum RELEVANT_FIELD {
  NAME = 'name',
  MOLFILE = 'molfile',
  CODES = 'codes',
}

export async function viewMonomerLib(): Promise<void> {
  const table = await parseMonomerLib(LIB_PATH, DEFAULT_LIB_FILENAME);
  table.name = 'Monomer Library';
  const view = grok.shell.addTableView(table);
  view.grid.props.allowEdit = false;
  const onDoubleClick = view.grid.onCellDoubleClick;
  onDoubleClick.subscribe(async (gridCell: DG.GridCell) => {
    const molfile = gridCell.cell.value;
    if (gridCell.tableColumn?.semType === 'Molecule')
      await drawZoomedInMolecule(molfile);
  });
}

async function parseMonomerLib(path: string, fileName: string): Promise<DG.DataFrame> {
  const fileSource = new DG.FileSource(path);
  const file = await fileSource.readAsText(fileName);
  const objList = JSON.parse(file);
  const formattedObjectsList = new Array(objList.length);
  for (let i = 0; i < objList.length; i++)
    formattedObjectsList[i] = formatMonomerObject(objList[i]);
  const df = DG.DataFrame.fromObjects(formattedObjectsList)!;
  return df;
}

function formatMonomerObject(sourceObj: ExtendedMonomer): FormattedMonomer {
  const formattedObject: FormattedMonomer = {};
  formattedObject[RELEVANT_FIELD.NAME] = sourceObj[RELEVANT_FIELD.NAME];
  formattedObject[RELEVANT_FIELD.MOLFILE] = sourceObj[RELEVANT_FIELD.MOLFILE];
  const codes = sourceObj[RELEVANT_FIELD.CODES] as CodesField;
  for (const synthesizer of Object.values(SYNTHESIZERS)) {
    const fieldName = synthesizer;
    const valuesList = [];
    // const technologySet = new Set();
    for (const technology of Object.values(TECHNOLOGIES)) {
      if (codes[synthesizer] !== undefined) {
        if (codes[synthesizer][technology] !== undefined) {
          valuesList.push(codes[synthesizer][technology].toString());
          // technologySet.add(technology);
        }
      }
    }
    // formattedObject['technologies'] = [...technologySet].toString();
    formattedObject[fieldName] = valuesList.toString();
  }
  return formattedObject;
}
