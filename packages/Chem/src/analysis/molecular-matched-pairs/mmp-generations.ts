import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import { MMP_COLNAME_FROM, MMP_COLNAME_TO, columnsDescriptions } from './mmp-constants';
import { FILTER_TYPES, chemSubstructureSearchLibrary } from '../../chem-searches';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { SubstructureSearchType } from '../../constants';

export function getGenerations(molecules: DG.Column, frags: [string, string][][],
  allPairsGrid: DG.Grid, activityMeanNames: Array<string>, activities: DG.ColumnList, module:RDModule):
  DG.Grid {
  const rulesColumns = allPairsGrid.dataFrame.columns;
  const rulesFrom = rulesColumns.byName(MMP_COLNAME_FROM);
  const rulesTo = rulesColumns.byName(MMP_COLNAME_TO);
  const activityN = activityMeanNames.length;

  const structures = molecules.toList();
  const cores = new Array<string>(activityN * structures.length);
  const from = new Array<string>(activityN * structures.length);
  const to = new Array<string>(activityN * structures.length);
  const prediction = new Array<number>(activityN * structures.length).fill(0);
  const allStructures = Array(structures.length * activityN);
  const allInitActivities: Array<number> = Array(structures.length * activityN);
  const activityName: Array<string> = Array(structures.length * activityN);

  for (let i = 0; i < structures.length; i ++) {
    const singleActivities = new Array<number>(activityN);
    for (let j = 0; j < activityN; j++) {
      const name = activityMeanNames[j].replace('Mean Difference ', '');
      const initActivity = activities.byName(name).get(i);
      singleActivities[j] = initActivity;
      allStructures[j * structures.length + i] = structures[i];
      allInitActivities[j * structures.length + i] = initActivity;
      activityName[j * structures.length + i] = name;
    }

    for (let j = 0; j < frags[i].length; j++) {
      const core = frags[i][j][0];
      const subst = frags[i][j][1];
      if (core != '') {
        for (let k = 0; k < rulesFrom.length; k++) {
          if (subst === rulesFrom.get(k)) {
            for (let kk = 0; kk < activityN; kk++) {
              const activity = singleActivities[kk] + rulesColumns.byName(activityMeanNames[kk]).get(k);
              if (activity > prediction[kk * structures.length + i]) {
                prediction[kk * structures.length + i] = activity;
                cores[kk * structures.length + i] = core;
                from[kk * structures.length + i] = subst;
                to[kk * structures.length + i] = rulesTo.get(k);
              }
            }
          }
        }
      }
    }
  }


  const generation = new Array<string>(cores.length);
  let mol = null;
  for (let i = 0; i < cores.length; i++) {
    let smilesGen = '';
    try {
      const smi = `${cores[i]}.${to[i]}`.replaceAll('([*:1])', '9').replaceAll('[*:1]', '9');
      mol = module.get_mol(smi);
      smilesGen = mol.get_smiles();
      generation[i] = smilesGen;
    } catch (e: any) {

    } finally {
      mol?.delete();
      generation[i] = smilesGen;
    }
  }

  const cols = [];
  cols.push(createColumnWithDescription('string', 'Structure', allStructures, DG.SEMTYPE.MOLECULE));
  cols.push(createColumnWithDescription('double', `Initial value`, allInitActivities));
  cols.push(createColumnWithDescription('string', `Activity`, activityName));
  cols.push(createColumnWithDescription('string', `Core`, cores, DG.SEMTYPE.MOLECULE));
  cols.push(createColumnWithDescription('string', `From`, from, DG.SEMTYPE.MOLECULE));
  cols.push(createColumnWithDescription('string', `To`, to, DG.SEMTYPE.MOLECULE));
  cols.push(createColumnWithDescription('double', `Prediction`, prediction));
  cols.push(createColumnWithDescription('string', `Generation`, generation, DG.SEMTYPE.MOLECULE));
  const grid = DG.DataFrame.fromColumns(cols).plot.grid();
  createIsMoleculeExistingCol(molecules, generation, grid);
  return grid;
}

export function createColumnWithDescription(colType: any, colName: string, list: any[], semType?: DG.SemType): DG.Column {
  const col = DG.Column.fromList(colType, colName, list);
  if (columnsDescriptions[colName])
    col.setTag('description', columnsDescriptions[colName]);
  if (semType)
    col.semType = semType;
  return col;
}

export async function createIsMoleculeExistingCol(molecules: DG.Column, generation: string[], grid: DG.Grid): Promise<void> {
  const progressBar = DG.TaskBarProgressIndicator.create(`Calculating generations...`);
  const promises = [];
  for (let i = 0; i < generation.length; i++)
    promises.push(await chemSubstructureSearchLibrary(molecules, generation[i], '', FILTER_TYPES.substructure,
      false, true, SubstructureSearchType.EXACT_MATCH));
  Promise.all(promises).then((res: BitArray[]) => {
    const molExistsColNew = createColumnWithDescription('bool', `Existing`, res.map((it, i) => generation[i] && !it.allFalse));
    grid.dataFrame.columns.add(molExistsColNew);
    grid.col(molExistsColNew.name)!.editable = false;
    grid.invalidate();
    progressBar.close();
  });
}
