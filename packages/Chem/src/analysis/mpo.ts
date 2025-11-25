/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DesirabilityProfile} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';

export class MpoProfileDialog {
  dataFrame: DG.DataFrame;
  mpoProfileEditor: MpoProfileEditor;
  templateInput: DG.ChoiceInput<string | null>;
  addParetoFront: DG.InputBase<boolean>;
  mpoFiles: DG.FileInfo[] = [];
  currentTemplate: DesirabilityProfile | null = null;
  currentTemplateFileName: string | null = null;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame ?? grok.shell.t;
    this.mpoProfileEditor = new MpoProfileEditor(this.dataFrame);

    this.templateInput = ui.input.choice('Template', {
      onValueChanged: async (value) => await this.loadProfile(value),
      nullable: false,
    });

    this.addParetoFront = ui.input.bool('Pareto front');
  }

  async init(): Promise<void> {
    this.mpoFiles = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
    const defaultTemplate = this.mpoFiles.length > 0 ? this.mpoFiles[0].fileName : null;

    this.templateInput.items = this.mpoFiles.map((f) => f.fileName);
    if (defaultTemplate) {
      this.templateInput.value = defaultTemplate;
      await this.loadProfile(defaultTemplate);
    }
  }

  private async loadProfile(fileName: string | null): Promise<void> {
    if (!fileName) return;
    this.currentTemplateFileName = fileName;
    const templateFile = this.mpoFiles.find((f) => f.fileName === fileName);
    if (!templateFile) return;

    try {
      const templateContent = await templateFile.readAsString();
      this.currentTemplate = JSON.parse(templateContent) as DesirabilityProfile;
      this.mpoProfileEditor.setProfile(this.currentTemplate);
    } catch (e) {
      grok.shell.error(`Failed to load template '${fileName}': ${e instanceof Error ? e.message : e}`);
    }
  }

  private async saveTemplate(): Promise<void> {
    if (!this.currentTemplateFileName || !this.currentTemplate) return;
    try {
      const updatedTemplateString = JSON.stringify(this.currentTemplate, null, 2);
      await grok.dapi.files.writeAsText(
        `${MPO_TEMPLATE_PATH}/${this.currentTemplateFileName}`,
        updatedTemplateString,
      );
      grok.shell.info(`Template '${this.currentTemplateFileName}' updated.`);
    } catch (e) {
      grok.shell.error(`Failed to save template '${this.currentTemplateFileName}': ${e instanceof Error ? e.message : e}`);
    }
  }

  private addParetoFrontViewer(columnNames: string[]) {
    const view = grok.shell.getTableView(this.dataFrame.name);
    const paretoFrontViewer = DG.Viewer.fromType('Pareto front', this.dataFrame);
    view.addViewer(paretoFrontViewer);

    /**
     * Temporary workaround to ensure the Pareto front viewer applies the provided
     * column names. The names are not applied immediately after the viewer
     * is created (issue under investigation).
    */
    setTimeout(() => paretoFrontViewer.setOptions({
      minimizeColumnNames: [],
      maximizeColumnNames: columnNames,
    }), 1000);
  }

  private async runMpoCalculation(): Promise<void> {
    try {
      const [func] = await DG.Func.find({name: 'mpoTransformFunction'});
      const funcCall = await func.prepare({
        df: this.dataFrame,
        currentProperties: this.currentTemplate?.properties,
      }).call(undefined, undefined, {processed: false});

      const columnNames: string[] = funcCall.getOutputParamValue();
      if (columnNames.length && this.addParetoFront.value)
        this.addParetoFrontViewer(columnNames);
    } catch (e) {
      grok.shell.error(`Failed to run MPO calculation: ${e instanceof Error ? e.message : e}`);
    }
  }

  getEditor(): HTMLElement {
    return ui.divV([
      this.templateInput.root,
      this.mpoProfileEditor.root,
      this.addParetoFront.root,
    ]);
  }

  async showDialog(): Promise<void> {
    await this.init();

    ui.dialog('MPO Score')
      .add(this.getEditor())
      .onOK(async () => {
        await this.saveTemplate();
        await this.runMpoCalculation();
      })
      .show();
  }
}
