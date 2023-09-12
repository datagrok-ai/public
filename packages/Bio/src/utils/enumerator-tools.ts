/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {NotationConverter} from '@datagrok-libraries/bio/src/utils/notation-converter';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibHelper} from '../utils/monomer-lib';
import {_package} from '../package';

const LEFT_HELM_WRAPPER = 'PEPTIDE1{';
const RIGHT_HELM_WRAPPER = '}$$$$';
const ALL_MONOMERS = '<All>';

function addCommonTags(col: DG.Column):void {
  col.setTag('quality', DG.SEMTYPE.MACROMOLECULE);
  col.setTag('aligned', ALIGNMENT.SEQ);
  col.setTag('alphabet', ALPHABET.PT);
}

export function _setPeptideColumn(col: DG.Column): void {
  addCommonTags(col);
  col.setTag('units', NOTATION.SEPARATOR);
  col.setTag('separator', '-');
  // col.setTag('cell.renderer', 'sequence');
}

async function enumerator(molColumn: DG.Column, leftTerminal: string, rightTerminal: string): Promise<void> {
  function hasSpecifiedTerminals(helm: string, leftTerminal: string, rightTerminal: string): boolean {
    if (leftTerminal === ALL_MONOMERS || rightTerminal === ALL_MONOMERS)
      return true;
    return helm.includes(LEFT_HELM_WRAPPER + leftTerminal) && helm.includes(rightTerminal + RIGHT_HELM_WRAPPER);
  }
  function getCyclicHelm(helm: string): string {
    if (!hasSpecifiedTerminals(helm, leftTerminal, rightTerminal))
      return helm;
    const seq = helm.replace(LEFT_HELM_WRAPPER, '').replace(RIGHT_HELM_WRAPPER, '');
    const lastMonomerNumber = seq.split('.').length;
    const result = helm.replace(RIGHT_HELM_WRAPPER,
      `}$PEPTIDE1,PEPTIDE1,${lastMonomerNumber}:R2-1:R1${'$'.repeat(6)}`);
    console.log('result:', result);
    return result;
  }

  const df = molColumn.dataFrame;
  const nc = new NotationConverter(molColumn);
  const sourceHelmCol = nc.convert(NOTATION.HELM);
  // df.columns.add(sourceHelmCol);
  const targetList = sourceHelmCol.toList().map((helm) => getCyclicHelm(helm));
  const colName = df.columns.getUnusedName('Enumerator(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', colName, targetList);
  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);
  targetHelmCol.setTag('cell.renderer', 'helm');
  df.columns.add(targetHelmCol);
  await grok.data.detectSemanticTypes(df);
}

export function _getEnumeratorWidget(molColumn: DG.Column): DG.Widget {
  const monomerLib = MonomerLibHelper.instance.getBioLib();
  const peptideList: string[] = [ALL_MONOMERS].concat(monomerLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE));
  const leftTerminalChoice = ui.choiceInput('', peptideList[0], peptideList);
  const rightTerminalChoice = ui.choiceInput('', peptideList[0], peptideList);

  const rGroups = ['term N', 'term O'];
  const r1groupChoice = ui.choiceInput('R1', rGroups[0], rGroups);
  const r2groupChoice = ui.choiceInput('R2', rGroups[1], rGroups);

  const modifications = ['Cyclization'];

  const modificationChoice = ui.choiceInput('Modification', modifications[0], modifications);

  const btn = ui.bigButton('Run', async () =>
    enumerator(molColumn, leftTerminalChoice.value!, rightTerminalChoice.value!)
  );

  const div = ui.div([
    modificationChoice,
    ui.divH([r1groupChoice.root, leftTerminalChoice.root]),
    ui.divH([r2groupChoice.root, rightTerminalChoice.root]),
    btn
  ]);

  return new DG.Widget(div);
}
