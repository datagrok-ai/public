import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isCurrentUserAppAdmin, stringify} from './helpers';
import {STORAGE_NAME, CURRENT_USER, ADDITIONAL_MODS_COL_NAMES, BASE_MODIFICATIONS, TOOLTIPS,
  EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION, USER_IS_NOT_ADMIN_MESSAGE} from './constants';


export async function addModification(modificationsDf: DG.DataFrame): Promise<void> {
  if (!await isCurrentUserAppAdmin()) return grok.shell.warning(USER_IS_NOT_ADMIN_MESSAGE);

  const longName = ui.stringInput(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES, '');
  ui.tooltip.bind(longName.root, TOOLTIPS.LONG_NAME);
  const abbreviation = ui.stringInput(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION, '');
  ui.tooltip.bind(abbreviation.root, TOOLTIPS.ABBREVIATIONS);
  const molecularWeight = ui.floatInput(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT, 0);
  const baseModification = ui.choiceInput(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION,
    BASE_MODIFICATIONS.NO, Object.values(BASE_MODIFICATIONS), (v: string) => {
      if (v != BASE_MODIFICATIONS.NO)
        extCoefficient.value = EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION;
      extCoefficient.enabled = (v == BASE_MODIFICATIONS.NO);
    });
  const extCoefficient = ui.stringInput(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT, '');

  ui.dialog('Add Modification')
    .add(ui.block([
      longName.root,
      abbreviation.root,
      molecularWeight.root,
      baseModification.root,
      extCoefficient.root,
    ]))
    .onOK(async () => {
      if (longName.value.length > 300)
        return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');
      if (abbreviation.value.length > 100)
        return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
      const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
      if (abbreviation.value in entries)
        return grok.shell.warning(`Abbreviation ${abbreviation.value} already exists`);
      const userObj = await grok.dapi.users.current();
      const modifiedLogs = Date() + ' by ' + userObj.firstName + ' ' + userObj.lastName + '; ';
      await grok.dapi.userDataStorage.postValue(
        STORAGE_NAME,
        abbreviation.value,
        JSON.stringify({
          longName: longName.value,
          abbreviation: abbreviation.value,
          molecularWeight: molecularWeight.value,
          extinctionCoefficient: extCoefficient.value,
          baseModification: baseModification.value,
          changeLogs: modifiedLogs,
        }),
        CURRENT_USER,
      ).then(() => grok.shell.info('New modification added'));

      modificationsDf.rows.addNew([
        longName.value, abbreviation.value, molecularWeight.value,
        baseModification.value, extCoefficient.value, modifiedLogs,
      ]);
    })
    .show();
}


export async function editModification(additionalModsDf: DG.DataFrame, additionaModifsGrid: DG.Grid): Promise<void> {
  additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES)!.width = 110;
  additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION)!.width = 80;
  additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT)!.width = 105;
  additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION)!.width = 110;
  additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT)!.width = 100;

  // Hide 'CHANGE_LOGS' column, display its content in tooltip
  additionaModifsGrid.columns.setVisible(additionalModsDf.columns.names().slice(0, -1));
  additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS).name = '~' + ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS;
  additionaModifsGrid.onCellTooltip(function(cell, x, y) {
    if (cell.isTableCell) {
      const v = additionalModsDf.getCol('~' + ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS)
        .get(cell.gridRow).split('; ').slice(0, -1);
      ui.tooltip.show(ui.divText(v), x, y);
      return true;
    }
  });

  additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION)
    .setTag(DG.TAGS.CHOICES, stringify(Object.values(BASE_MODIFICATIONS)));

  const col = additionaModifsGrid.col(ADDITIONAL_MODS_COL_NAMES.ACTION)!;
  col.cellType = 'html';
  additionaModifsGrid.onCellPrepare(function(gc) {
    if (gc.isTableCell && gc.gridColumn.name == ADDITIONAL_MODS_COL_NAMES.ACTION)
      gc.style.element = ui.button(ui.iconFA('trash-alt'), () => deleteModification(additionalModsDf, gc.gridRow));
  });

  let tempValue = '';
  additionalModsDf.onCurrentCellChanged.subscribe(() => {
    tempValue = additionalModsDf.currentCell.value;
  });

  DG.debounce(additionalModsDf.onValuesChanged, 10).subscribe(async (_) => {
    if (!await isCurrentUserAppAdmin())
      return grok.shell.warning(USER_IS_NOT_ADMIN_MESSAGE);
    if (additionalModsDf.currentCol.name == ADDITIONAL_MODS_COL_NAMES.ABBREVIATION) {
      const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
      if (additionalModsDf.currentCell.value.length > 100)
        return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
      if (additionalModsDf.currentCell.value in entries) {
        additionalModsDf.set(additionalModsDf.currentCol.name, additionalModsDf.currentRowIdx, tempValue);
        return grok.shell.warning(`Abbreviation ${additionalModsDf.currentCell.value} already exists`);
      }
    }
    if (additionalModsDf.currentCol.name == ADDITIONAL_MODS_COL_NAMES.LONG_NAMES &&
      additionalModsDf.currentCell.value.length > 300)
      return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');

    const rowIndex = additionalModsDf.currentCell.rowIndex;
    if (additionalModsDf.currentCol.name == ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION) {
      if (additionalModsDf.currentCell.value == BASE_MODIFICATIONS.NO) {
        const extCoefChoiceInput = ui.floatInput('', 0);
        ui.dialog('Enter Extinction Coefficient Value')
          .add(extCoefChoiceInput)
          .onOK(() => {
            const col = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT);
            col.set(rowIndex, String(extCoefChoiceInput.value), false);
            additionaModifsGrid.invalidate();
          })
          .show();
      } else {
        const col = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT);
        col.set(rowIndex, EXT_COEFF_VALUE_FOR_NO_BASE_MODIFICATION, false);
        additionaModifsGrid.invalidate();
      }
    }

    await grok.dapi.userDataStorage.postValue(
      STORAGE_NAME,
      additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).getString(rowIndex),
      JSON.stringify({
        longName: additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES).get(rowIndex),
        abbreviation: additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).get(rowIndex),
        molecularWeight: additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT).get(rowIndex),
        extinctionCoefficient: additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.EXTINCTION_COEFFICIENT).get(rowIndex),
        baseModification: additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION).get(rowIndex),
        changeLogs: additionalModsDf.getCol('~' + ADDITIONAL_MODS_COL_NAMES.CHANGE_LOGS).get(rowIndex),
      }),
      CURRENT_USER,
    );
  });
}


export async function deleteModification(additionalModsDf: DG.DataFrame, rowIndex: number): Promise<void> {
  if (!await isCurrentUserAppAdmin()) return grok.shell.warning(USER_IS_NOT_ADMIN_MESSAGE);
  const keyToDelete = additionalModsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).getString(rowIndex);
  ui.dialog(`Do you want to delete abbreviation ${keyToDelete}?`)
    .onOK(async () => {
      additionalModsDf.rows.removeAt(rowIndex, 1, true);
      await grok.dapi.userDataStorage.remove(STORAGE_NAME, keyToDelete, CURRENT_USER)
        .then(() => grok.shell.info(`${keyToDelete} deleted`));
    })
    .show();
}
