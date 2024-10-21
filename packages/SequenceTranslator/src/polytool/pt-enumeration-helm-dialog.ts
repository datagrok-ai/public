import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import wu from 'wu';
import {fromEvent, Unsubscribable} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {HelmAtom, HelmMol} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';
import {getHelmHelper, HelmInputBase} from '@datagrok-libraries/bio/src/helm/helm-helper';
import {getMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {HelmType, PolymerType} from '@datagrok-libraries/bio/src/helm/types';
import {helmTypeToPolymerType} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {getSeqHelper, ISeqHelper} from '@datagrok-libraries/bio/src/utils/seq-helper';
import '@datagrok-libraries/bio/src/types/input';
import {errInfo} from '@datagrok-libraries/bio/src/utils/err-info';
import {InputColumnBase} from '@datagrok-libraries/bio/src/types/input';

import {
  PolyToolEnumeratorParams, PolyToolEnumeratorType, PolyToolEnumeratorTypes
} from './types';
import {getLibrariesList} from './utils';
import {doPolyToolEnumerateHelm, PT_HELM_EXAMPLE} from './pt-enumeration-helm';
import {PolyToolPlaceholdersInput} from './pt-placeholders-input';
import {defaultErrorHandler} from '../utils/err-info';
import {PolyToolPlaceholdersBreadthInput} from './pt-placeholders-breadth-input';
import {PT_UI_DIALOG_ENUMERATION} from './const';

import {_package} from '../package';

type PolyToolEnumerateInputs = {
  macromolecule: HelmInputBase;
  placeholders: PolyToolPlaceholdersInput;
  placeholdersBreadth: PolyToolPlaceholdersBreadthInput;
  enumeratorType: DG.ChoiceInput<PolyToolEnumeratorType>
  trivialNameCol: InputColumnBase,
  keepOriginal: DG.InputBase<boolean>;
  toAtomicLevel: DG.InputBase<boolean>;
  toHarmonizedSequence: DG.InputBase<boolean>;
};

type PolyToolEnumerateHelmSerialized = {
  macromolecule: string;
  placeholders: string;
  placeholdersBreadth: string;
  enumeratorType: PolyToolEnumeratorType;
  trivialNameCol: string;
  keepOriginal: boolean;
  toAtomicLevel: boolean;
  toHarmonizedSequence: boolean;
};

export async function polyToolEnumerateHelmUI(cell?: DG.Cell): Promise<void> {
  const maxWidth = window.innerWidth;
  const maxHeight = window.innerHeight;

  try {
    let dialog: DG.Dialog;
    const resizeInputs = () => {
      if (dialog == null) return;

      const dialogContentEl = $(dialog.root).find('div.d4-dialog-contents').get(0)! as HTMLElement;
      const contentHeight = dialogContentEl.clientHeight;

      const fitInputs: { [idx: number]: number } = {0: 1 /*, 3: 0.5*/};
      const fitInputsSumHeight = Object.values(fitInputs).reduce((sum, h) => sum + h, 0);

      const otherInputsHeight: number = wu.count(0).take(dialogContentEl.children.length)
        .filter((i) => !(i in fitInputs))
        .map((idx) => dialogContentEl.children[idx])
        .filter((el) => el instanceof HTMLElement)
        .map((el) => (el as HTMLElement).offsetHeight).reduce((sum, h) => sum + h, 0);
      const remainFitHeight = contentHeight - otherInputsHeight - 38;
      dialog.inputs.forEach((input, idx) => {
        if (idx in fitInputs) {
          const inputFitHeight = remainFitHeight * fitInputs[idx] / fitInputsSumHeight;
          input.root.style.height = `${inputFitHeight}px`;
        }
      });
    };
    dialog = await getPolyToolEnumerateDialog(cell, resizeInputs);

    let isFirstShow = true;
    ui.onSizeChanged(dialog.root).subscribe(() => {
      if (isFirstShow) {
        const dialogInputList = dialog.inputs;
        const dialogRootCash = $(dialog.root);
        const contentMaxHeight = maxHeight -
          dialogRootCash.find('div.d4-dialog-header').get(0)!.offsetHeight -
          dialogRootCash.find('div.d4-dialog-footer').get(0)!.offsetHeight;

        // dialog.inputs2.macromolecule.root.style.backgroundColor = '#CCFFCC';

        const dialogWidth = maxWidth * 0.7;
        const dialogHeight = maxHeight * 0.7;

        // Centered, but resizable dialog
        dialog.root.style.width = `${Math.min(maxWidth, dialogWidth)}px`;
        dialog.root.style.height = `${Math.min(maxHeight, dialogHeight)}px`;
        dialog.root.style.left = `${Math.floor((maxWidth - dialog.root.offsetWidth) / 2)}px`;
        dialog.root.style.top = `${Math.floor((maxHeight - dialog.root.offsetHeight) / 2)}px`;

        isFirstShow = false;
      }

      resizeInputs();
    });
    resizeInputs();

    _package.logger.debug('PolyToolEnumerateHelmUI: dialog before show');
    const res = dialog.show({width: Math.max(350, maxWidth * 0.7), /* center: true,*/ resizable: true});
    _package.logger.debug('PolyToolEnumerateHelmUI: dialog after show');
  } catch (err: any) {
    const [errMsg, errStack] = errInfo(err);
    grok.shell.warning('To run PolyTool Enumeration, sketch the macromolecule and select monomers to vary');
    _package.logger.error(errMsg, undefined, errStack);
  }
}


async function getPolyToolEnumerateDialog(
  cell?: DG.Cell, resizeInputs?: () => void
): Promise<DG.Dialog> {
  const logPrefix = `ST: PT: HelmDialog()`;
  let inputs: PolyToolEnumerateInputs;
  const subs: Unsubscribable[] = [];
  const destroy = () => {
    for (const sub of subs) sub.unsubscribe();
    inputs.placeholders.detach();
  };
  try {
    const monomerLib = (await getMonomerLibHelper()).getMonomerLib();
    const seqHelper = await getSeqHelper();
    const emptyDf: DG.DataFrame = DG.DataFrame.fromColumns([]);

    const [libList, helmHelper] = await Promise.all([getLibrariesList(), getHelmHelper()]);

    let srcId: { value: string, colName: string } | null = null;
    const trivialNameSampleDiv = ui.divText('', {style: {marginLeft: '8px', marginTop: '2px'}});
    const warningsTextDiv = ui.divText('', {style: {color: 'red'}});
    inputs = {
      macromolecule: helmHelper.createHelmInput(
        'Macromolecule', {editable: false}),
      placeholders: await PolyToolPlaceholdersInput.create(
        'Placeholders', {
          showAddNewRowIcon: true,
          showRemoveRowIcon: true,
          showRowHeader: false,
          showCellTooltip: false,
        }),
      enumeratorType: ui.input.choice<PolyToolEnumeratorType>(
        'Enumerator type', {
          value: PolyToolEnumeratorTypes.Single,
          items: Object.values(PolyToolEnumeratorTypes)
        }) as DG.ChoiceInput<PolyToolEnumeratorType>,
      placeholdersBreadth: await PolyToolPlaceholdersBreadthInput.create(
        'Breadth', {
          showAddNewRowIcon: true,
          showRemoveRowIcon: true,
          showRowHeader: false,
          showCellTooltip: false,
        }),
      keepOriginal: ui.input.bool(
        'Keep original', {value: false}),
      toAtomicLevel: ui.input.bool(
        'To atomic level', {value: false}),
      toHarmonizedSequence: ui.input.bool(
        'To harmonized sequence', {value: false}),
      trivialNameCol: ui.input.column2(
        'Trivial name', {
          table: cell?.dataFrame,
          filter: (col: DG.Column): boolean => {
            return col.type === DG.COLUMN_TYPE.STRING && col != cell?.column; /* id */
          },
          onValueChanged: (): void => {
            const valueCol = inputs.trivialNameCol.value;
            let newSrcId: typeof srcId = null;
            if (cell && valueCol) {
              const originalId = valueCol.get(cell.rowIndex)!;
              newSrcId = {value: originalId, colName: valueCol.name};
            }
            srcId = newSrcId;
            trivialNameSampleDiv.textContent = srcId ? `Original ID: ${srcId.value}` : '';
          },
          nullable: true,
        }),
    };

    inputs.toHarmonizedSequence.root.style.display = 'none';
    inputs.trivialNameCol.addOptions(trivialNameSampleDiv);

    let placeholdersValidity: string | null = null;
    inputs.placeholders.addValidator((value: string): string | null => {
      const errors: string[] = [];
      try {
        const missedMonomerList: { polymerType: PolymerType, symbol: string }[] = [];
        for (const ph of inputs.placeholders.placeholdersValue) {
          const pos = ph.position;
          if (pos >= inputs.macromolecule.molValue.atoms.length) {
            errors.push(`There is no monomer at position ${pos + 1}.`);
            continue;
          }
          const a = inputs.macromolecule.molValue.atoms[pos];
          const helmType: HelmType = a.biotype()!;
          const polymerType = helmTypeToPolymerType(helmType);
          for (const symbol of ph.monomers) {
            const substituteMonomer = monomerLib.getMonomer(polymerType, symbol)!;
            // TODO: Check substitution monomer is presented in the library
            if (!substituteMonomer || !substituteMonomer.lib)
              missedMonomerList.push({polymerType, symbol});
          }
        }

        const byType: { [polymerType: string]: string[] } = {};
        for (const sm of missedMonomerList) {
          let byTypeList = byType[sm.polymerType];
          if (!byTypeList) byTypeList = byType[sm.polymerType] = [];
          byTypeList.push(sm.symbol);
        }
        const byTypeStr: string = Object.entries(byType)
          .map(([polymerType, symbolList]) => `${polymerType}: ${symbolList.join(', ')}`)
          .join('\n');
        if (Object.keys(byTypeStr).length > 0)
          errors.push(`Placeholders contain missed monomers: ${byTypeStr}`);
        placeholdersValidity = errors.length > 0 ? errors.join('\n') : null;
      } catch (err: any) {
        const [errMsg, errStack] = defaultErrorHandler(err, false);
        placeholdersValidity = errMsg;
      }
      setTimeout(() => { updateWarnings(); }, 0);
      return placeholdersValidity;
    });

    subs.push(inputs.macromolecule.onMouseMove.subscribe((e: MouseEvent) => {
      try {
        _package.logger.debug(`${logPrefix}, placeholdersInput.onMouseMove()`);

        const argsX = e.offsetX;
        const argsY = e.offsetY;
        const mol = inputs.macromolecule.molValue;
        const hoveredAtom = helmHelper.getHoveredAtom(argsX, argsY, mol, inputs.macromolecule.root.clientHeight);
        if (hoveredAtom) {
          const hoveredAtomContIdx = hoveredAtom._parent.atoms.indexOf(hoveredAtom);
          const substitutingMonomers = inputs.placeholders.placeholdersValue
            .find((ph) => ph.position === hoveredAtomContIdx)?.monomers;

          if (substitutingMonomers) {
            const cnt = ui.divText(substitutingMonomers.join(', '));
            inputs.macromolecule.showTooltip(cnt, hoveredAtom);
            e.preventDefault();
            e.stopPropagation();
          }
        }
      } catch (err: any) {
        defaultErrorHandler(err, false);
      }
    }));
    subs.push(inputs.macromolecule.onClick.subscribe((e: MouseEvent) => {
      try {
        _package.logger.debug(`${logPrefix}, placeholdersInput.onClick()`);

        const argsX = e.offsetX;
        const argsY = e.offsetY;
        const mol = inputs.macromolecule.molValue;
        const clickedAtom = helmHelper.getHoveredAtom(argsX, argsY, mol, inputs.macromolecule.root.clientHeight);
        if (clickedAtom) {
          const clickedAtomContIdx = clickedAtom._parent.atoms.indexOf(clickedAtom);
          const clickedAtomContIdxStr = (clickedAtomContIdx + 1).toString();

          const phDf = inputs.placeholders.grid.dataFrame;
          const posList = phDf.columns.byName('Position').toList();
          let rowIdx = posList.indexOf(clickedAtomContIdxStr);
          if (rowIdx === -1) {
            rowIdx = posList.findIndex((v) => isNaN(v));
            if (rowIdx === -1)
              rowIdx = phDf.rows.addNew([clickedAtomContIdxStr, '']).idx;

            phDf.set('Position', rowIdx, clickedAtomContIdxStr);
            // const tgtCell = inputs.placeholders.grid.cell('Monomers', rowIdx);
          }
          phDf.currentCell = phDf.cell(rowIdx, 'Monomers');
          //const gridRowIdx = inputs.placeholders.grid.tableRowToGrid(rowIdx);
          //const monomersGCell = inputs.placeholders.grid.cell('Monomers', gridRowIdx);
          const k = 42;
        }
      } catch (err: any) {
        defaultErrorHandler(err);
      }
    }));
    subs.push(inputs.placeholders.grid.dataFrame.onDataChanged.subscribe(() => {
      updateMolView();
    }));
    subs.push(fromEvent<KeyboardEvent>(inputs.placeholders.grid.root, 'keydown')
      .subscribe((e: KeyboardEvent) => {
        if (e.key === 'Enter') e.stopPropagation();
      }));

    // TODO: suspect
    subs.push(ui.onSizeChanged(inputs.placeholders.root).subscribe(() => {
      if (resizeInputs) resizeInputs();
    }));

    // Displays the molecule from a current cell (monitors changes)
    subs.push(grok.events.onCurrentCellChanged.subscribe(() => {
      const cell = grok.shell.tv.dataFrame.currentCell;
      if (cell.column.semType !== DG.SEMTYPE.MACROMOLECULE) return;

      fillForCurrentCell(cell);
    }));

    inputs.macromolecule.root.style.setProperty('min-width', '250px', 'important');
    // inputs.macromolecule.root.style.setProperty('max-height', '300px', 'important');

    const updateMolView = () => {
      const mol = inputs.macromolecule.molValue;
      for (let aI = 0; aI < mol.atoms.length; aI++) {
        const a = mol.atoms[aI];
        a.highlighted = aI in inputs.placeholders.placeholdersValue;
      }
      inputs.macromolecule.redraw();
    };

    const updateWarnings = () => {
      const warnings = placeholdersValidity;
      // const iw = inputs.warnings;
      const w = warningsTextDiv;
      if (!!warnings) {
        // iw.value = warnings; // <- breaks dialog resize
        // iw.enabled = true;
        // iw.root.style.removeProperty('display');

        w.innerText = warnings;
        w.style.removeProperty('display');
      } else {
        // iw.value = ''; // <- breaks dialog resize
        // iw.enabled = false;
        // iw.root.style.setProperty('display', 'none');

        w.innerText = '';
        w.style.setProperty('display', 'none');
      }
      //resizeInputs();
    };

    const fillTrivialNameList = (table?: DG.DataFrame) => {
      // const colList: DG.Column[] = [];
      // const colCount: number = cell.dataFrame.columns.length;
      //
      // // TODO: by semType?
      // for (let colI: number = 0; colI < colCount; ++colI) {
      //   const col = cell.dataFrame.columns.byIndex(colI);
      //   if (col.type === DG.COLUMN_TYPE.STRING) colList.push(col);
      // }

      if (table) {
        inputs.trivialNameCol.setColumnInputTable(table);
        inputs.trivialNameCol.root.style.removeProperty('display');
      } else {
        inputs.trivialNameCol.setColumnInputTable(emptyDf);
        inputs.trivialNameCol.root.style.setProperty('display', 'none');
      }
      if (resizeInputs) resizeInputs();
    };

    const fillForCurrentCell = async (cell?: DG.Cell): Promise<void> => {
      let helmValue: string;
      let table: DG.DataFrame | undefined = undefined;
      if (cell && cell.rowIndex >= 0 && cell?.column.semType == DG.SEMTYPE.MACROMOLECULE) {
        const sh = _package.seqHelper.getSeqHandler(cell.column);
        const helmSemValue = await sh.getHelm(cell.rowIndex);
        helmValue = helmSemValue.value;
        table = cell.dataFrame;
      } else {
        helmValue = PT_HELM_EXAMPLE;
      }

      inputs.macromolecule.stringValue = helmValue;
      fillTrivialNameList(table);
    };

    await fillForCurrentCell(cell);

    const exec = async (): Promise<void> => {
      try {
        const srcHelm = inputs.macromolecule.stringValue;
        const helmSelections: number[] = wu.enumerate<HelmAtom>(inputs.macromolecule.molValue.atoms)
          .filter(([a, aI]) => a.highlighted)
          .map(([a, aI]) => aI).toArray();
        if (srcHelm === undefined || srcHelm === '') {
          grok.shell.warning('PolyTool: no molecule was provided');
        } else /* if (helmSelections === undefined || helmSelections.length < 1) {
          grok.shell.warning('PolyTool: no selection was provided');
        } else /**/ {
          if (Object.keys(inputs.placeholders.placeholdersValue).length === 0 &&
            Object.keys(inputs.placeholdersBreadth.placeholdersBreadthValue).length === 0
          ) {
            grok.shell.warning(`${PT_UI_DIALOG_ENUMERATION}: placeholders are empty`);
            return;
          }
          await getHelmHelper(); // initializes JSDraw and org
          const params: PolyToolEnumeratorParams = {
            placeholders: inputs.placeholders.placeholdersValue,
            type: inputs.enumeratorType.value!,
            placeholdersBreadth: inputs.placeholdersBreadth.placeholdersBreadthValue,
            keepOriginal: inputs.keepOriginal.value,
          };
          const enumeratorResDf = await polyToolEnumerateHelm(srcHelm, srcId, params, inputs.toAtomicLevel.value, seqHelper);
          grok.shell.addTableView(enumeratorResDf);
        }
      } catch (err: any) {
        defaultErrorHandler(err);
      }
    };

    const dialog = ui.dialog({title: PT_UI_DIALOG_ENUMERATION, showFooter: true})
      .add(inputs.macromolecule)
      .add(inputs.placeholders)
      .add(inputs.enumeratorType)
      .add(inputs.placeholdersBreadth)
      .add(inputs.trivialNameCol)
      .add(ui.divH([inputs.keepOriginal.root, inputs.toAtomicLevel.root, inputs.toHarmonizedSequence.root],
        {style: {width: '100%'}}))
      .add(warningsTextDiv)
      // .addButton('Enumerate', () => {
      //   execDialog()
      //     .then(() => {});
      // }, 0, 'Keeps the dialog open')
      .onOK(() => { exec(); });
    subs.push(dialog.onClose.subscribe(() => {
      destroy();
    }));
    dialog.history(
      /* getInput */ (): PolyToolEnumerateHelmSerialized => {
        return {
          macromolecule: inputs.macromolecule.stringValue,
          placeholders: inputs.placeholders.stringValue,
          enumeratorType: inputs.enumeratorType.value,
          placeholdersBreadth: inputs.placeholdersBreadth.stringValue,
          trivialNameCol: inputs.trivialNameCol.stringValue,
          keepOriginal: inputs.keepOriginal.value,
          toAtomicLevel: inputs.toAtomicLevel.value,
          toHarmonizedSequence: inputs.toHarmonizedSequence.value,
        };
      },
      /* applyInput */ (x: PolyToolEnumerateHelmSerialized): void => {
        inputs.macromolecule.stringValue = x.macromolecule;
        inputs.placeholders.stringValue = x.placeholders;
        inputs.enumeratorType.value = x.enumeratorType;
        inputs.placeholdersBreadth.stringValue = x.placeholdersBreadth;
        inputs.trivialNameCol.stringValue = x.trivialNameCol;
        inputs.keepOriginal.value = x.keepOriginal;
        inputs.toAtomicLevel.value = x.toAtomicLevel;
        inputs.toHarmonizedSequence.value = x.toHarmonizedSequence;
      });
    return dialog;
  } catch (err: any) {
    destroy(); // on failing to build a dialog
    throw err;
  }
}

async function polyToolEnumerateHelm(
  srcHelm: string, srcId: { value: string, colName: string } | null, params: PolyToolEnumeratorParams,
  toAtomicLevel: boolean, seqHelper: ISeqHelper,
): Promise<DG.DataFrame> {
  const pi = DG.TaskBarProgressIndicator.create('PolyTool enumerating...');
  try {
    await getHelmHelper(); // initializes JSDraw and org

    const resList = doPolyToolEnumerateHelm(srcHelm, srcId?.value ?? '', params);
    const enumHelmCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, 'Enumerated', resList.length)
      .init((rowIdx: number) => resList[rowIdx][0]);
    const enumeratorResDf = DG.DataFrame.fromColumns([enumHelmCol]);

    if (toAtomicLevel) {
      const seqHelper: ISeqHelper = await getSeqHelper();
      const toAtomicLevelRes = await seqHelper.helmToAtomicLevel(enumHelmCol, true, true);
      toAtomicLevelRes.molCol!.semType = DG.SEMTYPE.MOLECULE;
      enumeratorResDf.columns.add(toAtomicLevelRes.molCol!, false);
    }

    if (srcId) {
      const enumIdCol = DG.Column.fromType(DG.COLUMN_TYPE.STRING, srcId.colName, resList.length)
        .init((rowIdx: number) => resList[rowIdx][1]);
      enumeratorResDf.columns.add(enumIdCol);
    }

    return enumeratorResDf;
  } finally {
    pi.close();
  }
}
