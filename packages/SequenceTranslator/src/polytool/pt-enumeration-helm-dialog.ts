import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {Unsubscribable} from 'rxjs';

import {NOTATION} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {HelmAtom, HelmMol} from '@datagrok-libraries/helm-web-editor/src/types/org-helm';
import {getHelmHelper, HelmInputBase} from '@datagrok-libraries/bio/src/helm/helm-helper';

import {
  PolyToolEnumeratorParams, PolyToolEnumeratorType, PolyToolEnumeratorTypes, PolyToolPlaceholders
} from './types';

import {PT_UI_DIALOG_ENUMERATION} from './const';
import {getLibrariesList} from './utils';
import {getPtEnumeratorHelm, PT_HELM_EXAMPLE} from './pt-enumeration-helm';
import {PolyToolPlaceholdersInput} from './pt-placeholders-input';
import {Dialog} from 'datagrok-api/dg';
import {defaultErrorHandler} from '../utils/err-info';

import {_package} from '../package';

type PolyToolEnumerateInputs = {
  enumeratorType: DG.ChoiceInput<PolyToolEnumeratorType>
  macromolecule: HelmInputBase;
  placeholders: PolyToolPlaceholdersInput;
  toAtomicLevel: DG.InputBase<boolean>;
};


export class PolyToolEnumerateDialog extends DG.Dialog {
  protected constructor(
    // public readonly inputs2: PolyToolEnumerateInputs
  ) {
    const dlg = ui.dialog({title: PT_UI_DIALOG_ENUMERATION});
    super(dlg.dart);

    // for (const [key, value] of Object.entries(this.inputs2)) { this.add(value); }
  }

  public override show(options?: { modal?: boolean; resizable?: boolean; fullScreen?: boolean; center?: boolean; centerAt?: Element; x?: number; y?: number; width?: number; height?: number; backgroundColor?: string; showNextTo?: HTMLElement }): Dialog {
    return super.show(options);
  }

  public static async create2(cell?: DG.Cell, resizeInputs?: () => void): Promise<[PolyToolEnumerateDialog, PolyToolEnumerateInputs]> {
    const logPrefix = `ST: PT: HelmDialog()`;

    const [libList, helmHelper] = await Promise.all([getLibrariesList(), getHelmHelper()]);

    const helmValue = cell ? cell.value : PT_HELM_EXAMPLE;
    const posDf = DG.DataFrame.fromObjects([{Position: '', Monomers: ''}, {Position: '', Monomers: ''}]);

    const macromoleculeInput = helmHelper.createHelmInput(
      'Macromolecule', {value: helmValue, editable: false});

    const inputs: PolyToolEnumerateInputs = {
      enumeratorType: ui.input.choice<PolyToolEnumeratorType>(
        'Enumerator type', {
          value: PolyToolEnumeratorTypes.Single,
          items: Object.values(PolyToolEnumeratorTypes)
        }) as DG.ChoiceInput<PolyToolEnumeratorType>,
      macromolecule: macromoleculeInput,

      placeholders: await PolyToolPlaceholdersInput.create(
        'Placeholders', posDf, {
          showAddNewRowIcon: true,
          showRemoveRowIcon: true,
          showRowHeader: false,
          showCellTooltip: false,
        }/*, 2/**/),

      toAtomicLevel: ui.input.bool(
        'To atomic level', {value: true}),
    };

    const subs: Unsubscribable[] = [];
    const destroy = () => { for (const sub of subs) sub.unsubscribe(); };

    subs.push(inputs.macromolecule.onMouseMove.subscribe((e: MouseEvent) => {
      try {
        _package.logger.debug(`${logPrefix}, placeholdersInput.onMouseMove()`);

        const argsX = e.offsetX;
        const argsY = e.offsetY;
        const mol = inputs.macromolecule.molValue;
        const hoveredAtom = helmHelper.getHoveredAtom(argsX, argsY, mol, inputs.macromolecule.root.clientHeight);
        if (hoveredAtom) {
          const hoveredAtomContIdx = hoveredAtom._parent.atoms.indexOf(hoveredAtom);
          const hoveredAtomContIdxStr = (hoveredAtomContIdx + 1).toString();
          const substitutingMonomers = inputs.placeholders.placeholdersValue[hoveredAtomContIdx];

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


    // TODO: suspect
    subs.push(ui.onSizeChanged(inputs.placeholders.root).subscribe(() => {
      if (resizeInputs) resizeInputs();
    }));

    // Displays the molecule from a current cell (monitors changes)
    subs.push(grok.events.onCurrentCellChanged.subscribe(() => {
      const cell = grok.shell.tv.dataFrame.currentCell;

      if (cell.column.semType === DG.SEMTYPE.MACROMOLECULE && cell.column.meta.units === NOTATION.HELM)
        inputs.macromolecule.stringValue = cell.value;
    }));

    inputs.macromolecule.root.style.setProperty('min-width', '250px', 'important');
    // inputs.macromolecule.root.style.setProperty('max-height', '300px', 'important');

    const updateMolView = () => {
      const mol = inputs.macromolecule.molValue;
      for (let aI = 0; aI < mol.atoms.length; aI++) {
        const a = mol.atoms[aI];
        a.highlighted = aI in inputs.placeholders.placeholdersValue;
        inputs.macromolecule.redraw();
      }
    };

    const dialog = new PolyToolEnumerateDialog()
      .add(inputs.enumeratorType)
      .add(inputs.macromolecule)
      .add(inputs.placeholders)
      .add(inputs.toAtomicLevel)
      .onOK(async () => {
        try {
          const helmString = inputs.macromolecule.stringValue;
          const helmSelections: number[] = wu.enumerate<HelmAtom>(inputs.macromolecule.molValue.atoms)
            .filter(([a, aI]) => a.highlighted)
            .map(([a, aI]) => aI).toArray();
          if (helmString === undefined || helmString === '') {
            grok.shell.warning('PolyTool: no molecule was provided');
          } else /* if (helmSelections === undefined || helmSelections.length < 1) {
          grok.shell.warning('PolyTool: no selection was provided');
        } else /**/ {
            await getHelmHelper(); // initializes JSDraw and org

            const params: PolyToolEnumeratorParams = {
              type: inputs.enumeratorType.value!,
              placeholders: inputs.placeholders.placeholdersValue,
            };
            if (Object.keys(params.placeholders).length === 0) {
              grok.shell.warning(`${PT_UI_DIALOG_ENUMERATION}: placeholders are empty`);
              return;
            }

            const enumHelmList = getPtEnumeratorHelm(helmString, params);
            const enumHelmCol = DG.Column.fromStrings('Enumerated', enumHelmList);
            const enumeratorResDf = DG.DataFrame.fromColumns([enumHelmCol]);

            if (inputs.toAtomicLevel.value) {
              const enumMolCol = await grok.functions.call('Bio:getMolFromHelm', {
                'df': enumeratorResDf,
                'helmCol': enumHelmCol,
                'chiralityEngine': true,
              });
              enumMolCol.name = enumeratorResDf.columns.getUnusedName(`molfile(${enumHelmCol.name})`);
              enumMolCol.semType = DG.SEMTYPE.MOLECULE;
              enumeratorResDf.columns.add(enumMolCol, true);
            }
            grok.shell.addTableView(enumeratorResDf);
          }
        } catch (err: any) {
          defaultErrorHandler(err);
        } finally {
          destroy();
        }
      })
      .onCancel(() => { destroy(); }) as PolyToolEnumerateDialog;
    return [dialog, inputs];
  }
}
