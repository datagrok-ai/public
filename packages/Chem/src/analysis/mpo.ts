import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';
import {
  DesirabilityProfile,
  WeightedAggregation,
  WEIGHTED_AGGREGATIONS_LIST,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {MpoContextPanel} from '../mpo/mpo-context-panel';
import {MpoProfileManager} from '../mpo/mpo-profile-manager';
import {PackageFunctions} from '../package';
import {computeMpo, MpoProfileInfo, deepEqual, findSuitableProfiles} from '../mpo/utils';

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
  private suitableProfileNames: string[] = [];
  private subs: Subscription[] = [];
  private footerWarning: HTMLElement | null = null;

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
    this.mpoProfiles = await MpoProfileManager.ensureLoaded();
    this.suitableProfileNames = findSuitableProfiles(this.dataFrame, this.mpoProfiles).map((p) => p.fileName);

    this.profileInput.items = this.mpoProfiles.map((p) => p.fileName);
    requestAnimationFrame(() => this.highlightSuitableProfiles(this.suitableProfileNames));

    const defaultProfile = this.suitableProfileNames.length > 0 ?
      this.suitableProfileNames[0] :
      this.mpoProfiles[0]?.fileName ?? null;

    if (defaultProfile) {
      this.profileInput.value = defaultProfile;
      await this.loadProfile(defaultProfile);
    }

    this.listenForProfileChanges();
    this.mpoContextPanel = new MpoContextPanel(this.dataFrame);
  }

  private highlightSuitableProfiles(suitableFileNames: string[]): void {
    const STAR_PLACEHOLDER = '\u2003\u2002';
    const select = this.profileInput.input as HTMLSelectElement;

    for (let i = 0; i < select.options.length; i++) {
      const option = select.options[i];
      const text = option.textContent ?? '';
      if (text.startsWith('⭐') || text.startsWith('\u2003'))
        continue;

      const isApplicable = suitableFileNames.includes(option.value);
      option.textContent = `${isApplicable ? '⭐' : STAR_PLACEHOLDER} ${text}`;
    }
  }

  private listenForProfileChanges(): void {
    this.subs.push(grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      this.updateSaveButtonVisibility();
      this.updateOkButtonState();
      if (this.currentProfile && this.mpoContextPanel) {
        await this.mpoContextPanel.render(
          this.currentProfile,
          this.mpoProfileEditor.columnMapping,
          this.aggregationInput.value ?? undefined,
        );
      }
    }));
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
    this.updateOkButtonState();
  }

  private isProfileModified(): boolean {
    if (!this.currentProfile || !this.originalProfile)
      return false;
    return !deepEqual(this.currentProfile, this.originalProfile);
  }

  private updateSaveButtonVisibility(): void {
    this.saveButton.classList.toggle('chem-mpo-d-none', !this.isProfileModified());
  }

  private isProfileApplicable(): boolean {
    if (!this.currentProfile)
      return false;

    const mapping = this.mpoProfileEditor.columnMapping;
    const dfColumnNames = this.dataFrame.columns.names();
    const propertyNames = Object.keys(this.currentProfile.properties);

    if (propertyNames.length === 0)
      return false;

    return propertyNames.every((prop) => {
      if (mapping[prop] != null)
        return true;
      if (prop in mapping)
        return false;
      return dfColumnNames.includes(prop);
    });
  }

  private updateOkButtonState(): void {
    const isApplicable = this.isProfileApplicable();
    const okButton = this.dialog?.getButton('OK');
    if (okButton)
      okButton.disabled = !isApplicable;
    this.updateFooterWarning(!isApplicable);
  }

  private updateFooterWarning(show: boolean): void {
    if (!this.footerWarning) {
      const commandBar = this.dialog?.root.querySelector('.d4-command-bar');
      if (!commandBar)
        return;
      this.footerWarning = ui.divText('⚠️ Some profile properties are not mapped', 'chem-mpo-footer-warning');
      commandBar.append(this.footerWarning);
    }
    this.footerWarning.classList.toggle('chem-mpo-d-none', !show);
  }

  private async saveProfile(): Promise<void> {
    if (!this.currentProfileFileName || !this.currentProfile)
      return;

    await MpoProfileManager.save(this.currentProfile, this.currentProfileFileName);
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
        this.detach();
      })
      .onCancel(() => this.detach())
      .show();

    this.updateOkButtonState();
  }

  private detach(): void {
    this.mpoContextPanel?.close();
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }
}
