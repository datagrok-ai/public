
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {NOTATION, ALIGNMENT, ALPHABET} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {UnitsHandler} from '@datagrok-libraries/bio/src/utils/units-handler';
import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibHelper} from '../utils/monomer-lib';
import {_package} from '../package';
import * as rxjs from 'rxjs';

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
  function hasTerminals(helm: string, leftTerminal: string, rightTerminal: string): boolean {
    if (leftTerminal === ALL_MONOMERS || rightTerminal === ALL_MONOMERS)
      return true;
    const positions = getLinkedPositions(helm);
    return positions.every((el) => el > 0);
  }

  function applyModification(helm: string): string {
    if (hasTerminals(helm, leftTerminal, rightTerminal))
      return getCycle(helm, getLinkedPositions(helm));
    return helm;
  }

  function getLinkedPositions(helm: string): [number, number] {
    const seq = helm.replace(LEFT_HELM_WRAPPER, '').replace(RIGHT_HELM_WRAPPER, '');
    const monomers = seq.split('.');
    const start = monomers.findIndex((el) => el === leftTerminal);
    const end = monomers.findIndex((el, idx) => el === rightTerminal && idx > start);
    return [start + 1, end + 1];
  }

  function getCycle(helm: string, position: [number, number]): string {
    const result = helm.replace(RIGHT_HELM_WRAPPER,
      `}$PEPTIDE1,PEPTIDE1,${position[0]}:R3-${position[1]}:R3${'$'.repeat(6)}`);
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
  function updateMonomerList(): void {
    console.log('hi from update:');
    if (cyclizationTypeChoice.value === cyclizationTypes[0]) {
      console.log('hi from first branch:');
      monomerList = [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE)
      );
    } else if (cyclizationTypeChoice.value === cyclizationTypes[1]) {
      monomerList = [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByRGroup(3, HELM_POLYMER_TYPE.PEPTIDE)
      );
      console.log('hi from second branch:');
    }
    leftTerminalChoice = ui.choiceInput('R1:', monomerList[0], monomerList);
    rightTerminalChoice = ui.choiceInput('R2:', monomerList[0], monomerList);
    ui.empty(terminalControls);
    [leftTerminalChoice, rightTerminalChoice].forEach((el) => { terminalControls.appendChild(el.root); });
  }

  const onCyclizationChoice = new rxjs.Subject<string>();
  onCyclizationChoice.subscribe(() => updateMonomerList());

  const modifications = ['Cyclization'];
  const modificationChoice = ui.choiceInput('Modification', modifications[0], modifications);

  const cyclizationTypes = ['N-O', 'R3-R3'];
  const cyclizationTypeChoice = ui.choiceInput(
    'Type', cyclizationTypes[0], cyclizationTypes, () => { onCyclizationChoice.next(); }
  );

  const monomerLib = MonomerLibHelper.instance.getBioLib();
  let monomerList: string[] = [];
  let leftTerminalChoice = ui.choiceInput('R1:', monomerList[0], monomerList);
  let rightTerminalChoice = ui.choiceInput('R2:', monomerList[0], monomerList);
  const terminalControls = ui.divV([leftTerminalChoice.root, rightTerminalChoice.root]);

  updateMonomerList();

  const btn = ui.bigButton('Run', async () =>
    enumerator(molColumn, leftTerminalChoice.value!, rightTerminalChoice.value!)
  );

  const div = ui.div([
    modificationChoice,
    cyclizationTypeChoice,
    terminalControls,
    btn
  ]);

  return new DG.Widget(div);
}
