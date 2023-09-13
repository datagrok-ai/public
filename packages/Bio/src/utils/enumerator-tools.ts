
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibHelper} from '../utils/monomer-lib';
import {toAtomicLevel} from '../package';
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

async function enumerator(
  molColumn: DG.Column, leftTerminal: string, rightTerminal: string): Promise<void> {
  function isFullCycle(helm: string, leftTerminal: string, rightTerminal: string): boolean {
    if (leftTerminal === ALL_MONOMERS || rightTerminal === ALL_MONOMERS)
      return true;
    return helm.includes(LEFT_HELM_WRAPPER + leftTerminal) && helm.includes(rightTerminal + RIGHT_HELM_WRAPPER);
  }

  function applyModification(helm: string): string {
    if (isFullCycle(helm, leftTerminal, rightTerminal))
      return getCycle(helm, getLinkedPositions(helm));
    return helm;
  }

  function getLinkedPositions(helm: string): [number, number] {
    const seq = helm.replace(LEFT_HELM_WRAPPER, '').replace(RIGHT_HELM_WRAPPER, '');
    const lastMonomerNumber = seq.split('.').length;
    return [1, lastMonomerNumber];
  }

  function getCycle(helm: string, position: [number, number]): string {
    const result = helm.replace(RIGHT_HELM_WRAPPER,
      `}$PEPTIDE1,PEPTIDE1,${position[1]}:R2-${position[0]}:R1${'$'.repeat(6)}`);
    // console.log('result:', result);
    return result;
  }

  const df = molColumn.dataFrame;
  const uh = UnitsHandler.getOrCreate(molColumn);
  const sourceHelmCol = uh.convert(NOTATION.HELM);
  const targetList = sourceHelmCol.toList().map((helm) => applyModification(helm));
  const colName = df.columns.getUnusedName('Enumerator(' + molColumn.name + ')');
  const targetHelmCol = DG.Column.fromList('string', colName, targetList);

  addCommonTags(targetHelmCol);
  targetHelmCol.setTag('units', NOTATION.HELM);
  targetHelmCol.setTag('cell.renderer', 'helm');

  df.columns.add(targetHelmCol);
  // if (getAtomic)
  //   toAtomicLevel(df, targetHelmCol);
  await grok.data.detectSemanticTypes(df);
}

export function _getEnumeratorWidget(molColumn: DG.Column): DG.Widget {
  const monomerLib = MonomerLibHelper.instance.getBioLib();
  // const monomerList: string[] = [ALL_MONOMERS].concat(monomerLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE));
  const monomerList: string[] = [ALL_MONOMERS].concat(
    monomerLib.getMonomerSymbolsByRGroup(3, HELM_POLYMER_TYPE.PEPTIDE)
  );
  const leftTerminalChoice = ui.choiceInput('', monomerList[0], monomerList);
  const rightTerminalChoice = ui.choiceInput('', monomerList[0], monomerList);

  // const rGroups = ['N', 'O'];
  const rGroups = ['S', 'S'];
  const r1groupChoice = ui.choiceInput('R1', rGroups[0], rGroups);
  const r2groupChoice = ui.choiceInput('R2', rGroups[1], rGroups);

  const modifications = ['Cyclization'];

  const modificationChoice = ui.choiceInput('Modification', modifications[0], modifications);

  const selectAtomicStructure = ui.boolInput('Get mols', false);
  const btn = ui.bigButton('Run', async () =>
    enumerator(molColumn, leftTerminalChoice.value!, rightTerminalChoice.value!)
  );

  const div = ui.div([
    modificationChoice,
    ui.divH([r1groupChoice.root, leftTerminalChoice.root]),
    ui.divH([r2groupChoice.root, rightTerminalChoice.root]),
    // ui.divH([btn, atomicStructureInput.root]),
    btn
  ]);

  return new DG.Widget(div);
}
