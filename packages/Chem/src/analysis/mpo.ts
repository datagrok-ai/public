import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DesirabilityProfile,
  WeightedAggregation,
  WEIGHTED_AGGREGATIONS_LIST,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {MpoContextPanel} from '../mpo/mpo-context-panel';
import {PackageFunctions} from '../package';
import {computeMpo, MPO_TEMPLATE_PATH, loadMpoProfiles, MpoProfileInfo,
  deepEqual, findSuitableProfiles} from '../mpo/utils';

export class MpoProfileDialog {
  dataFrame: DG.DataFrame;
  mpoProfileEditor: MpoProfileEditor;
  aggregationInput: DG.ChoiceInput<WeightedAggregation | null>;
  profileInput: DG.ChoiceInput<string | null>;
  designModeInput: DG.InputBase<boolean>;
  addParetoFront: DG.InputBase<boolean>;
  mpoProfiles: MpoProfileInfo[] = [];
  currentProfile: DesirabilityProfile | null = null;
  currentProfileFileName: string | null = null;
  manageButton: HTMLElement;
  saveButton: HTMLButtonElement;

  private dialog?: DG.Dialog<{}>;
  private mpoContextPanel?: MpoContextPanel;
  private originalProfile: DesirabilityProfile | null = null;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame ?? grok.shell.t;
    this.mpoProfileEditor = new MpoProfileEditor(this.dataFrame);

    this.aggregationInput = ui.input.choice('Aggregation', {
      items: WEIGHTED_AGGREGATIONS_LIST,
      nullable: false,
      onValueChanged: () => grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {}),
    });

    this.profileInput = ui.input.choice('Profile', {
      nullable: false,
      onValueChanged: async (value) => await this.loadProfile(value),
    });

    this.designModeInput = ui.input.toggle('Design mode', {
      value: false,
      onValueChanged: (v) => this.mpoProfileEditor.setDesignMode(!!v),
    });

    this.addParetoFront = ui.input.bool('Pareto front');

    this.manageButton = ui.button('Manage...', async () => {
      grok.shell.addView(await PackageFunctions.mpoProfilesApp());
      this.dialog?.close();
    });
    this.manageButton.classList.add('chem-mpo-dialog-manage-button');

    this.saveButton = ui.button('Save', async () => {
      await this.saveProfile();
      this.originalProfile = structuredClone(this.currentProfile!);
      this.updateSaveButtonVisibility();
    });
    this.saveButton.classList.add('chem-mpo-dialog-manage-button', 'chem-mpo-d-none');
  }

  async init(): Promise<void> {
    this.mpoProfiles = await loadMpoProfiles();
    const suitableProfiles = findSuitableProfiles(this.dataFrame, this.mpoProfiles).map((p) => p.fileName);

    this.profileInput.items = this.mpoProfiles.map((p) => p.fileName);
    requestAnimationFrame(() => this.highlightSuitableProfiles(suitableProfiles));

    const defaultProfile = suitableProfiles.length > 0 ?
      suitableProfiles[0] :
      this.mpoProfiles[0]?.fileName ?? null;

    if (defaultProfile) {
      this.profileInput.value = defaultProfile;
      await this.loadProfile(defaultProfile);
    }

    this.listenForProfileChanges();
    this.mpoContextPanel = new MpoContextPanel(this.dataFrame);
  }

  private highlightSuitableProfiles(suitableFileNames: string[]): void {
    const select = this.profileInput.input as HTMLSelectElement;
    for (let i = 0; i < select.options.length; i++) {
      const option = select.options[i];
      const text = option.textContent ?? '';

      if (suitableFileNames.includes(option.value) && !text.startsWith('⭐'))
        option.textContent = `⭐ ${text}`;
    }
  }

  private listenForProfileChanges(): void {
    grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      this.updateSaveButtonVisibility();
      if (this.currentProfile && this.mpoContextPanel) {
        await this.mpoContextPanel.render(
          this.currentProfile,
          this.mpoProfileEditor.columnMapping,
          this.aggregationInput.value ?? undefined,
        );
      }
    });
  }

  private async loadProfile(fileName: string | null): Promise<void> {
    if (!fileName)
      return;

    const profileInfo = this.mpoProfiles.find((p) => p.fileName === fileName);
    if (!profileInfo)
      return;

    this.currentProfileFileName = fileName;
    this.currentProfile = structuredClone({
      name: profileInfo.name,
      description: profileInfo.description,
      properties: profileInfo.properties,
    });
    this.mpoProfileEditor.setProfile(this.currentProfile);
    this.originalProfile = structuredClone(this.currentProfile);
    this.updateSaveButtonVisibility();
  }

  private isProfileModified(): boolean {
    if (!this.currentProfile || !this.originalProfile)
      return false;
    return !deepEqual(this.currentProfile, this.originalProfile);
  }

  private updateSaveButtonVisibility(): void {
    this.saveButton.classList.toggle('chem-mpo-d-none', !this.isProfileModified());
  }

  private async saveProfile(): Promise<void> {
    if (!this.currentProfileFileName || !this.currentProfile)
      return;

    try {
      const updatedProfileString = JSON.stringify(this.currentProfile, null, 2);
      await grok.dapi.files.writeAsText(
        `${MPO_TEMPLATE_PATH}/${this.currentProfileFileName}`,
        updatedProfileString,
      );
      grok.shell.info(`Profile '${this.currentProfileFileName}' updated.`);
    } catch (e) {
      grok.shell.error(
        `Failed to save profile '${this.currentProfileFileName}': ${e instanceof Error ? e.message : e}`);
    }
  }

  private addParetoFrontViewer(columnNames: string[]): void {
    const view = grok.shell.getTableView(this.dataFrame.name);
    const paretoFrontViewer = DG.Viewer.fromType('Pareto front', this.dataFrame);
    view.addViewer(paretoFrontViewer);

    // Temporary workaround: set column names after viewer creation
    setTimeout(() => paretoFrontViewer.setOptions({
      minimizeColumnNames: [],
      maximizeColumnNames: columnNames,
    }), 1000);
  }

  private async runMpoCalculation(): Promise<void> {
    try {
      const columnNames = await computeMpo(
        this.dataFrame,
        this.currentProfile!,
        this.mpoProfileEditor.columnMapping,
        this.aggregationInput.value!,
      );

      if (columnNames.length && this.addParetoFront.value)
        this.addParetoFrontViewer(columnNames);
    } catch (e) {
      grok.shell.error(`Failed to run MPO calculation: ${e instanceof Error ? e.message : e}`);
    }
  }

  private getProfileControls(): HTMLElement {
    return ui.divH([this.profileInput.root, this.manageButton, this.saveButton], {style: {gap: '10px'}});
  }

  private getMpoControls(): HTMLElement {
    return ui.divV([
      this.aggregationInput.root,
      this.designModeInput.root,
      this.mpoProfileEditor.root,
      this.addParetoFront.root,
    ]);
  }

  getEditor(): HTMLElement {
    return ui.divV([
      this.getProfileControls(),
      this.getMpoControls(),
    ]);
  }

  async showDialog(): Promise<void> {
    await this.init();

    this.dialog = ui.dialog('MPO Score')
      .add(this.getEditor())
      .onOK(async () => {
        await this.runMpoCalculation();
        this.mpoContextPanel?.detach();
      })
      .onCancel(() => this.mpoContextPanel?.detach())
      .show();
  }
}
