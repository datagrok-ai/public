/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {HELM_POLYMER_TYPE} from '@datagrok-libraries/bio/src/utils/const';
import {MonomerLibHelper} from '../../utils/monomer-lib';
import {ALL_MONOMERS, CYCLIZATION_TYPE, TRANSFORMATION_TYPE} from './const';
import {addTransformedColumn} from './transformation';
import * as rxjs from 'rxjs';
import {MetaData} from './types';

export function getPolyToolDialog(): DG.Dialog {
  function getMonomerList(cyclizationType: CYCLIZATION_TYPE): string[] {
    if (cyclizationType === cyclizationTypes[0]) {
      return [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByType(HELM_POLYMER_TYPE.PEPTIDE)
      );
    }
    if (cyclizationType === cyclizationTypes[1]) {
      return [ALL_MONOMERS].concat(
        monomerLib.getMonomerSymbolsByRGroup(3, HELM_POLYMER_TYPE.PEPTIDE)
      );
    }
    return ['C'];
  }

  function updateMonomerList(): void {
    if (cyclizationTypeChoice.value === CYCLIZATION_TYPE.NCys) {
      monomerList1 = getMonomerList(CYCLIZATION_TYPE.NO);
      monomerList2 = getMonomerList(CYCLIZATION_TYPE.NCys);
    } else {
      monomerList1 = getMonomerList(cyclizationTypeChoice.value as CYCLIZATION_TYPE);
      monomerList2 = [...monomerList1];
    }

    leftTerminalChoice = ui.choiceInput(
      'R1:', monomerList1[0], monomerList1, () => { onRGroupValueChange.next(); }
    );
    rightTerminalChoice = ui.choiceInput('R2:', monomerList2[0], monomerList2, () => { onRGroupValueChange.next(); });
    onRGroupValueChange.next();
    ui.empty(terminalControls);
    [leftTerminalChoice, rightTerminalChoice].forEach((el) => { terminalControls.appendChild(el.root); });
  }

  function updateMeta() {
    meta.cyclizationType = cyclizationTypeChoice.value!;
    meta.leftTerminal = leftTerminalChoice.value!;
    meta.rightTerminal = rightTerminalChoice.value!;
    meta.transformationType = transformationChoice.value!;
  }


  const onCyclizationChoice = new rxjs.Subject<string>();
  const onRGroupValueChange = new rxjs.Subject<string>();
  onCyclizationChoice.subscribe(() => {
    meta.cyclizationType = cyclizationTypeChoice.value!;
    updateMonomerList();
  });
  onRGroupValueChange.subscribe(() => {
    meta.rightTerminal = rightTerminalChoice.value!;
    meta.leftTerminal = leftTerminalChoice.value!;
  });

  const meta = {} as MetaData;
  const transformations = [TRANSFORMATION_TYPE.CYCLIZATION];
  const transformationChoice = ui.choiceInput(
    'Modification', transformations[0], transformations, () => meta.transformationType = transformationChoice.value!
  );

  const cyclizationTypes = [CYCLIZATION_TYPE.NO, CYCLIZATION_TYPE.R3, CYCLIZATION_TYPE.NCys];
  const cyclizationTypeChoice = ui.choiceInput(
    'Type', cyclizationTypes[2], cyclizationTypes, () => { onCyclizationChoice.next(); }
  );

  const monomerLib = MonomerLibHelper.instance.getBioLib();
  let monomerList1: string[] = [];
  let monomerList2: string[] = [];
  let leftTerminalChoice = ui.choiceInput(
    'R1:', monomerList1[0], monomerList1, () => {
      meta.leftTerminal = leftTerminalChoice.value!;
    }
  );
  let rightTerminalChoice = ui.choiceInput('R2:', monomerList2[0], monomerList2, () => {
    meta.rightTerminal = rightTerminalChoice.value!;
  });
  const terminalControls = ui.divV([leftTerminalChoice.root, rightTerminalChoice.root]);
  updateMonomerList();

  updateMeta();

  const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
  if (!targetColumns)
    throw new Error('No dataframe with maceomolecule columns open');

  const targetColumnInput = ui.columnInput(
    'Column', grok.shell.t, targetColumns[0], null,
    {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE}
  );

  const generateHelmChoiceInput = ui.boolInput('Get HELM', true);
  ui.tooltip.bind(generateHelmChoiceInput.root, 'Add HELM column');

  const div = ui.div([
    targetColumnInput,
    transformationChoice,
    cyclizationTypeChoice,
    terminalControls,
    generateHelmChoiceInput,
  ]);

  const dialog = ui.dialog('Poly Tool')
    .add(div)
    .onOK(async () => {
      const molCol = targetColumnInput.value;
      if (!molCol) {
        grok.shell.warning('No marcomolecule column chosen!');
        return;
      }
      addTransformedColumn(molCol!, meta, generateHelmChoiceInput.value!);
    }
    );

  return dialog;
}
