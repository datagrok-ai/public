/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DesirabilityProfile, PropertyDesirability, WeightedAggregation, WEIGHTED_AGGREGATIONS} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

export const MPO_TEMPLATE_PATH = 'System:AppData/Chem/mpo';

export class MpoProfileDialog {
  dataFrame: DG.DataFrame;
  mpoProfileEditor: MpoProfileEditor;
  aggregationInput: DG.ChoiceInput<WeightedAggregation | null>;
  profileInput: DG.ChoiceInput<string | null>;
  designModeInput: DG.InputBase<boolean>;
  addParetoFront: DG.InputBase<boolean>;
  mpoFiles: DG.FileInfo[] = [];
  currentProfile: DesirabilityProfile | null = null;
  currentProfileFileName: string | null = null;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame ?? grok.shell.t;
    this.mpoProfileEditor = new MpoProfileEditor(this.dataFrame);

    this.aggregationInput = ui.input.choice('Aggregation', {items: [...WEIGHTED_AGGREGATIONS], nullable: false});

    this.profileInput = ui.input.choice('Profile', {
      onValueChanged: async (value) => await this.loadProfile(value),
      nullable: false,
    });

    this.designModeInput = ui.input.bool('Design mode', {value: false, onValueChanged: (v) => this.mpoProfileEditor.setDesignMode(!!v)});
    this.addParetoFront = ui.input.bool('Pareto front');
  }

  async init(): Promise<void> {
    this.mpoFiles = await grok.dapi.files.list(MPO_TEMPLATE_PATH);
    const defaultProfile = this.mpoFiles.length > 0 ? this.mpoFiles[0].fileName : null;

    this.profileInput.items = this.mpoFiles.map((f) => f.fileName);
    if (defaultProfile) {
      this.profileInput.value = defaultProfile;
      await this.loadProfile(defaultProfile);
    }
  }

  private async loadProfile(fileName: string | null): Promise<void> {
    if (!fileName) return;
    this.currentProfileFileName = fileName;
    const profileFile = this.mpoFiles.find((f) => f.fileName === fileName);
    if (!profileFile) return;

    try {
      const profileContent = await profileFile.readAsString();
      this.currentProfile = JSON.parse(profileContent) as DesirabilityProfile;
      this.mpoProfileEditor.setProfile(this.currentProfile);
    } catch (e) {
      grok.shell.error(`Failed to load profile '${fileName}': ${e instanceof Error ? e.message : e}`);
    }
  }

  private async saveProfile(): Promise<void> {
    if (!this.currentProfileFileName || !this.currentProfile) return;
    try {
      const updatedProfileString = JSON.stringify(this.currentProfile, null, 2);
      await grok.dapi.files.writeAsText(
        `${MPO_TEMPLATE_PATH}/${this.currentProfileFileName}`,
        updatedProfileString,
      );
      grok.shell.info(`Profile '${this.currentProfileFileName}' updated.`);
    } catch (e) {
      grok.shell.error(`Failed to save profile '${this.currentProfileFileName}': ${e instanceof Error ? e.message : e}`);
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
      const mapping = this.mpoProfileEditor.columnMapping;

      const mappedProperties: { [key: string]: PropertyDesirability } = {};
      for (const [propName, prop] of Object.entries(this.currentProfile!.properties)) {
        const columnName = mapping[propName] ?? propName;
        mappedProperties[columnName] = prop;
      }

      const [func] = await DG.Func.find({name: 'mpoTransformFunction'});
      const funcCall = await func.prepare({
        df: this.dataFrame,
        profileName: this.currentProfile?.name ?? 'MPO',
        currentProperties: mappedProperties,
        aggregation: this.aggregationInput.value,
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
      this.profileInput.root,
      this.aggregationInput.root,
      this.mpoProfileEditor.root,
      this.designModeInput.root,
      this.addParetoFront.root,
    ]);
  }

  async showDialog(): Promise<void> {
    await this.init();

    ui.dialog('MPO Score')
      .add(this.getEditor())
      .onOK(async () => {
        await this.saveProfile();
        await this.runMpoCalculation();
      })
      .show();
  }
}
