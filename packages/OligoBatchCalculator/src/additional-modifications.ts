import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {isCurrentUserAppAdmin} from './helpers';
import {STORAGE_NAME, CURRENT_USER, ADDITIONAL_MODS_COL_NAMES, BASE_MODIFICATIONS} from './constants';

export async function addModificationButton(modificationsDf: DG.DataFrame): Promise<void> {
  if (!await isCurrentUserAppAdmin())
    return grok.shell.info('You don\'t have permission for this action');

  const longName = ui.stringInput(ADDITIONAL_MODS_COL_NAMES.LONG_NAMES, '');
  ui.tooltip.bind(longName.root, 'Examples: \'Inverted Abasic\', \'Cyanine 3 CPG\', \'5-Methyl dC\'');
  const abbreviation = ui.stringInput(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION, '');
  ui.tooltip.bind(abbreviation.root, 'Examples: \'invabasic\', \'Cy3\', \'5MedC\'');
  const molecularWeight = ui.floatInput(ADDITIONAL_MODS_COL_NAMES.MOLECULAR_WEIGHT, 0);
  const baseModification = ui.choiceInput(ADDITIONAL_MODS_COL_NAMES.BASE_MODIFICATION,
    BASE_MODIFICATIONS.NO, Object.values(BASE_MODIFICATIONS), (v: string) => {
      if (v != BASE_MODIFICATIONS.NO)
        extCoefficient.value = 'Base';
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
        return grok.shell.warning('Abbreviation ' + abbreviation.value + ' already exists');
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
      ).then(() => grok.shell.info('Posted'));
      modificationsDf.rows.addNew([
        longName.value, abbreviation.value, molecularWeight.value,
        baseModification.value, extCoefficient.value, modifiedLogs,
      ]);
    })
    .show();
}

export function deleteAdditionalModification(additionalModificationsDf: DG.DataFrame, rowIndex: number): void {
  const v = additionalModificationsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).getString(rowIndex);
  ui.dialog('Delete Additional Modification')
    .add(ui.divText('Do you want to delete row ' + v + ' ?'))
    .onOK(async () => {
      if (!await isCurrentUserAppAdmin())
        return grok.shell.info('You don\'t have permission for this action');
      additionalModificationsDf.rows.removeAt(rowIndex, 1, true);
      const keyToDelete = additionalModificationsDf.getCol(ADDITIONAL_MODS_COL_NAMES.ABBREVIATION).get(rowIndex);
      await grok.dapi.userDataStorage.remove(STORAGE_NAME, keyToDelete, CURRENT_USER)
        .then(() => grok.shell.info('Deleted'));
    })
    .show();
}
