/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {addTransformedColumn} from './transformation';
import {RULES_PATH, RULES_STORAGE_NAME} from './transformation';

export type UserRuleSettings = {
  included: string[],
  notIncluded: string[],
}

async function getAllAvailableRuleFiles(): Promise<string[]> {
  const list = await grok.dapi.files.list(RULES_PATH);
  const paths = list.map((fileInfo) => {
    return fileInfo.fullPath;
  });

  return paths;
}

async function getUserRulesSettings(): Promise<UserRuleSettings> {
  const resStr: string = await grok.dapi.userDataStorage.getValue(RULES_STORAGE_NAME, 'Settings', true);
  const res = resStr ? JSON.parse(resStr) : {included: [], enotIncludedxplicit: []};

  res.included = res.included instanceof Array ? res.included : [];
  res.notIncluded = res.notIncluded instanceof Array ? res.notIncluded : [];

  return res!;
}

async function setUserLibSettings(value: UserRuleSettings): Promise<void> {
  await grok.dapi.userDataStorage.postValue(RULES_STORAGE_NAME, 'Settings', JSON.stringify(value), true);
}

export class PolyTool {
  ruleFiles: string[];
  userRuleSettings: UserRuleSettings;
  ruleFilesInputs: HTMLDivElement;// DG.InputBase<boolean | null>[];
  dialog: DG.Dialog;

  constructor() {
    this.ruleFiles = [];
    this.userRuleSettings = {included: [], notIncluded: []};
  }

  private updateRulesSelectionStatus(ruleFileName: string, isSelected: boolean): void {
    const isRuleFileSelected = this.userRuleSettings.included.includes(ruleFileName);

    if (!isRuleFileSelected && isSelected) {
      this.userRuleSettings.included.push(ruleFileName);
      this.userRuleSettings.included = this.userRuleSettings.included.sort();

      const index = this.userRuleSettings.notIncluded.indexOf(ruleFileName);
      if (index > -1)
        this.userRuleSettings.notIncluded.splice(index, 1);
    } else {
      const index = this.userRuleSettings.included.indexOf(ruleFileName);
      if (index > -1)
        this.userRuleSettings.included.splice(index, 1);

      this.userRuleSettings.notIncluded.push(ruleFileName);
      this.userRuleSettings.notIncluded = this.userRuleSettings.notIncluded.sort();
    }

    setUserLibSettings(this.userRuleSettings);
  }

  private getAddButton(): HTMLButtonElement {
    return ui.button('ADD RULES', () => {
      DG.Utils.openFile({
        accept: '.csv',
        open: async (selectedFile) => {
          const content = await selectedFile.text();
          await grok.dapi.files.writeAsText(RULES_PATH + `${selectedFile.name}`, content);
          this.updateRulesSelectionStatus(selectedFile.name, false);
          const cb = ui.boolInput(
            selectedFile.name,
            false,
            (isSelected: boolean) => this.updateRulesSelectionStatus(RULES_PATH + `${selectedFile.name}`, isSelected)
          );
          this.ruleFilesInputs.append(cb.root);
        },
      });
    });
  }

  private async getRuleFilesBlock(): Promise<DG.InputBase<boolean | null>[]> {
    this.ruleFiles = await getAllAvailableRuleFiles();
    this.userRuleSettings = await getUserRulesSettings();
    const cBoxes: DG.InputBase<boolean | null>[] = [];

    for (let i = 0; i < this.ruleFiles.length; i++) {
      const ruleFileName = this.ruleFiles[i];
      const isRuleFileSelected = this.userRuleSettings.included.includes(ruleFileName);
      const cb = ui.boolInput(
        ruleFileName.replace(RULES_PATH, ''),
        isRuleFileSelected,
        (isSelected: boolean) => this.updateRulesSelectionStatus(ruleFileName, isSelected)
      );

      cBoxes.push(cb);
    }
    return cBoxes;
  }

  async getPolyToolDialog(): Promise<DG.Dialog> {
    const targetColumns = grok.shell.t.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE);
    if (!targetColumns)
      throw new Error('No dataframe with macromolecule columns open');

    const targetColumnInput = ui.columnInput(
      'Column', grok.shell.t, targetColumns[0], null,
      {filter: (col: DG.Column) => col.semType === DG.SEMTYPE.MACROMOLECULE}
    );

    const generateHelmChoiceInput = ui.boolInput('Get HELM', true);
    ui.tooltip.bind(generateHelmChoiceInput.root, 'Add HELM column');

    const addButton = this.getAddButton();
    this.ruleFilesInputs = ui.div(await this.getRuleFilesBlock());
    //const rulesFiles = ui.div(this.ruleFilesInputs);

    const div = ui.div([
      targetColumnInput,
      generateHelmChoiceInput,
      'Rules used',
      this.ruleFilesInputs,
      addButton
    ]);

    this.dialog = ui.dialog('Poly Tool')
      .add(div)
      .onOK(async () => {
        const molCol = targetColumnInput.value;
        if (!molCol) {
          grok.shell.warning('No marcomolecule column chosen!');
          return;
        }
        addTransformedColumn(molCol!, generateHelmChoiceInput.value!, this.userRuleSettings.included);
      });

    return this.dialog;
  }
}
