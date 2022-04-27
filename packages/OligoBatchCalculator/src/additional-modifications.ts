import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

const CURRENT_USER = false;
const STORAGE_NAME = 'oligo-batch-calculator-storage';
export const COL_NAMES = {
  LONG_NAMES: 'Long name',
  ABBREVIATION: 'Abbreviation',
  MOLECULAR_WEIGHT: 'Molecular weight',
  BASE_MODIFICATION: 'Base modification',
  EXTINCTION_COEFFICIENT: 'Ext. coefficient',
  ACTION: 'Action',
  CHANGE_LOGS: 'Change logs',
};
const ADMIN_USERS = ['Baozhong Zhao', 'Sijin Guo', 'Saika Siddiqui', 'Vadym Kovadlo'];

export async function getAdditionalModifications(): Promise<DG.DataFrame> {
  const modifications: any[] = [];
  const entries = await grok.dapi.userDataStorage.get(STORAGE_NAME, CURRENT_USER);
  if (entries != null && Object.keys(entries).length == 0)
    grok.shell.info('Storage is empty. Try to post something to the storage');
  else {
    for (const key of Object.keys(entries))
      modifications.push(JSON.parse(entries[key]));
  }
  const molWeightList = modifications.map((e) => (e.molecularWeight == undefined) ? 0 : e.molecularWeight);
  const extinctionCoefList = modifications.map((e) => String(e.extinctionCoefficient));
  return DG.DataFrame.fromColumns([
    DG.Column.fromStrings(COL_NAMES.LONG_NAMES, modifications.map((e) => e.longName)),
    DG.Column.fromStrings(COL_NAMES.ABBREVIATION, modifications.map((e) => e.abbreviation)), // @ts-ignore
    DG.Column.fromFloat32Array(COL_NAMES.MOLECULAR_WEIGHT, molWeightList),
    DG.Column.fromStrings(COL_NAMES.BASE_MODIFICATION, modifications.map((e) => e.baseModification)),
    DG.Column.fromStrings(COL_NAMES.EXTINCTION_COEFFICIENT, extinctionCoefList),
    DG.Column.fromStrings(COL_NAMES.ACTION, Array(modifications.length)),
    DG.Column.fromStrings(COL_NAMES.CHANGE_LOGS, modifications.map((e) => e.changeLogs)),
  ])!;
}

export async function addModificationButton(modificationsDf: DG.DataFrame): Promise<void> {
  grok.dapi.users.current().then((user) => {
    if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
      const longName = ui.stringInput(COL_NAMES.LONG_NAMES, '');
      ui.tooltip.bind(longName.root, 'Examples: \'Inverted Abasic\', \'Cyanine 3 CPG\', \'5-Methyl dC\'');
      const abbreviation = ui.stringInput(COL_NAMES.ABBREVIATION, '');
      ui.tooltip.bind(abbreviation.root, 'Examples: \'invabasic\', \'Cy3\', \'5MedC\'');
      const molecularWeight = ui.floatInput(COL_NAMES.MOLECULAR_WEIGHT, 0);
      const baseModification = ui.choiceInput(COL_NAMES.BASE_MODIFICATION,
        'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
          if (v != 'NO')
            extCoefficient.value = 'Base';
          extCoefficient.enabled = (v == 'NO');
        });
      const extCoefficient = ui.stringInput(COL_NAMES.EXTINCTION_COEFFICIENT, '');
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
          const modifiedLogs = Date() + ' by ' + user.firstName + ' ' + user.lastName + '; ';
          grok.dapi.userDataStorage.postValue(
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
            extCoefficient.value, baseModification.value, modifiedLogs,
          ]);
        })
        .show();
    } else
      grok.shell.info('You don\'t have permission for this action');
  });
}

export function deleteAdditionalModification(additionalModificationsDf: DG.DataFrame, rowIndex: number): void {
  ui.dialog(
    'Do you want to delete ' + additionalModificationsDf.col(COL_NAMES.ABBREVIATION)!.getString(rowIndex) + ' ?',
  )
    .onOK(() => {
      grok.dapi.users.current().then(async (user) => {
        if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
          additionalModificationsDf.rows.removeAt(rowIndex, 1, true);
          const keyToDelete = additionalModificationsDf.col(COL_NAMES.ABBREVIATION)!.get(rowIndex);
          await grok.dapi.userDataStorage.remove(STORAGE_NAME, keyToDelete, CURRENT_USER);
        } else
          grok.shell.info('You don\'t have permission for this action');
      });
    })
    .show();
}

export function editAdditionalModification(additionalModificationsDf: DG.DataFrame, rowIndex: number): void {
  grok.dapi.users.current().then(async (user) => {
    if (ADMIN_USERS.includes(user.firstName + ' ' + user.lastName)) {
      const longName = ui.stringInput(COL_NAMES.LONG_NAMES,
        additionalModificationsDf.col(COL_NAMES.LONG_NAMES)!.get(rowIndex));
      ui.tooltip.bind(longName.root, 'Examples: \'Inverted Abasic\', \'Cyanine 3 CPG\', \'5-Methyl dC\'');
      const oldAbbreviation = additionalModificationsDf.col(COL_NAMES.ABBREVIATION)!.get(rowIndex);
      const abbreviation = ui.stringInput(COL_NAMES.ABBREVIATION, oldAbbreviation);
      ui.tooltip.bind(abbreviation.root, 'Examples: \'invabasic\', \'Cy3\', \'5MedC\'');
      const molecularWeight = ui.floatInput(COL_NAMES.MOLECULAR_WEIGHT,
        additionalModificationsDf.col(COL_NAMES.MOLECULAR_WEIGHT)!.get(rowIndex));
      const baseModification = ui.choiceInput(COL_NAMES.BASE_MODIFICATION,
        'NO', ['NO', 'rU', 'rA', 'rC', 'rG', 'dA', 'dC', 'dG', 'dT'], (v: string) => {
          if (v != 'NO')
            extinctionCoefficient.value = 'Base';
          extinctionCoefficient.enabled = (v == 'NO');
        });
      const extinctionCoefficient = ui.stringInput(COL_NAMES.EXTINCTION_COEFFICIENT,
        additionalModificationsDf.col(COL_NAMES.EXTINCTION_COEFFICIENT)!.get(rowIndex));
      const changeLogsCol = additionalModificationsDf.col(COL_NAMES.CHANGE_LOGS)!;
      ui.dialog('Edit Modification')
        .add(ui.block([
          longName.root,
          abbreviation.root,
          molecularWeight.root,
          baseModification.root,
          extinctionCoefficient.root,
        ]))
        .onOK(async () => {
          if (longName.value.length > 300)
            return grok.shell.warning('Long Name shouldn\'t contain more than 300 characters');
          if (abbreviation.value.length > 100)
            return grok.shell.warning('Abbreviation shouldn\'t contain more than 100 characters');
          const newLog = changeLogsCol.get(rowIndex) + Date() + ' by ' + user.firstName + ' ' + user.lastName + '; ';
          await grok.dapi.userDataStorage.postValue(
            STORAGE_NAME,
            abbreviation.value,
            JSON.stringify({
              longName: longName.value,
              abbreviation: abbreviation.value,
              molecularWeight: molecularWeight.value,
              extinctionCoefficient: extinctionCoefficient.value,
              baseModification: baseModification.value,
              changeLogs: newLog,
            }),
            CURRENT_USER,
          );
          if (oldAbbreviation != abbreviation.value)
            await grok.dapi.userDataStorage.remove(STORAGE_NAME, oldAbbreviation, CURRENT_USER);
          additionalModificationsDf.set(COL_NAMES.LONG_NAMES, rowIndex, longName.value);
          additionalModificationsDf.set(COL_NAMES.ABBREVIATION, rowIndex, abbreviation.value);
          additionalModificationsDf.set(COL_NAMES.MOLECULAR_WEIGHT, rowIndex, molecularWeight.value);
          additionalModificationsDf.set(COL_NAMES.EXTINCTION_COEFFICIENT, rowIndex, extinctionCoefficient.value);
          additionalModificationsDf.set(COL_NAMES.BASE_MODIFICATION, rowIndex, baseModification.value);
          additionalModificationsDf.set(COL_NAMES.CHANGE_LOGS, rowIndex, newLog);
        })
        .show();
    } else
      grok.shell.info('You don\'t have permission for this action');
  });
}
