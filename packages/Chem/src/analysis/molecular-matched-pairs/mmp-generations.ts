import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';

export function getGenerations(molecules: DG.Column, frags: [string, string][][],
  allPairsGrid: DG.Grid, activityMeanNames: Array<string>, activities: DG.ColumnList, module:RDModule):
  DG.Grid {
  const rulesColumns = allPairsGrid.dataFrame.columns;
  const rulesFrom = rulesColumns.byName('From');
  const rulesTo = rulesColumns.byName('To');
  const activityN = activityMeanNames.length;

  const structures = molecules.toList();
  const cores = new Array<Array<string>>(activityN);
  const from = new Array<Array<string>>(activityN);
  const to = new Array<Array<string>>(activityN);
  const prediction = new Array<Array<number>>(activityN);

  for (let i = 0; i < activityN; i ++) {
    cores[i] = new Array<string>(structures.length);
    from[i] = new Array<string>(structures.length);
    to[i] = new Array<string>(structures.length);
    prediction[i] = new Array<number>(structures.length);
    for (let j = 0; j < structures.length; j++)
      prediction[i][j] = 0;
  }

  for (let i = 0; i < structures.length; i ++) {
    const singleActivities = new Array<number>(activityN);
    for (let j = 0; j < activityN; j++)
      singleActivities[j] = activities.byName(activityMeanNames[j].replace('Mean Difference ', '')).get(i);

    for (let j = 0; j < frags[i].length; j++) {
      const core = frags[i][j][0];
      const subst = frags[i][j][1];
      if (core != '') {
        for (let k = 0; k < rulesFrom.length; k++) {
          if (subst === rulesFrom.get(k)) {
            for (let kk = 0; kk <activityN; kk++) {
              const activity = singleActivities[kk] + rulesColumns.byName(activityMeanNames[kk]).get(k);
              if (activity > prediction[kk][i]) {
                prediction[kk][i] = activity;
                cores[kk][i] = core;
                from[kk][i] = subst;
                to[kk][i] = rulesTo.get(k);
              }
            }
          }
        }
      }
    }
  }

  const allStructures = Array(structures.length * activityN).fill(0).map((_, i) => structures[i % structures.length]);

  const colStructure = DG.Column.fromList('string', 'structure', allStructures);
  colStructure.semType = DG.SEMTYPE.MOLECULE;
  const cols = [colStructure];
  let allCores: Array<string> = [];
  let allFrom: Array<string> = [];
  let allTo: Array<string> = [];
  let allInits: Array<number> = [];
  let activityName: Array<string> = [];
  let allPreds: Array<number> = [];

  for (let i = 0; i < activityN; i ++) {
    const name = activityMeanNames[i].replace('Mean Difference ', '');
    allInits = allInits.concat(activities.byName(name).toList());
    activityName = activityName.concat(Array(structures.length).fill(name));

    allCores = allCores.concat(cores[i]);
    allFrom = allFrom.concat(from[i]);
    allTo = allTo.concat(to[i]);
    allPreds = allPreds.concat(prediction[i]);
  }

  const generation = new Array<string>(allCores.length);
  let mol = null;
  for (let i = 0; i < allCores.length; i++) {
    let smilesGen = '';
    try {
      const smi = `${allCores[i]}.${allTo[i]}`.replaceAll('([*:1])', '9').replaceAll('[*:1]', '9');
      mol = module.get_mol(smi);
      smilesGen = mol.get_smiles();
      generation[i] = smilesGen;
    } catch (e: any) {

    } finally {
      mol?.delete();
      generation[i] = smilesGen;
    }
  }

  cols.push(DG.Column.fromList('double', `Initial value`, allInits));
  cols.push(DG.Column.fromList('string', `Activity`, activityName));
  const coreCol = DG.Column.fromList('string', `Core`, allCores);
  coreCol.semType = DG.SEMTYPE.MOLECULE;
  cols.push(coreCol);
  const fromCol = DG.Column.fromList('string', `From`, allFrom);
  fromCol.semType = DG.SEMTYPE.MOLECULE;
  cols.push(fromCol);
  const toCol = DG.Column.fromList('string', `To`, allTo);
  toCol.semType = DG.SEMTYPE.MOLECULE;
  cols.push(toCol);
  cols.push(DG.Column.fromList('double', `Prediction`, allPreds));
  const generationCol = DG.Column.fromList('string', `Generation`, generation);
  generationCol.semType = DG.SEMTYPE.MOLECULE;
  cols.push(generationCol);

  const generations = DG.DataFrame.fromColumns(cols);
  return generations.plot.grid();
}
