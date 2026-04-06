import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {Subscription} from 'rxjs';
import {
  createDefaultNumerical,
  DESIRABILITY_PROFILE_TYPE,
  DesirabilityProfile,
  PropertyDesirability,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MPO_SCORE_CHANGED_EVENT} from '@datagrok-libraries/statistics/src/mpo/utils';

import {MpoContextPanel} from '../mpo/mpo-context-panel';
import {MpoProfileManager} from '../mpo/mpo-profile-manager';
import {PackageFunctions} from '../package';
import {computeMpo, MpoProfileInfo, deepEqual, findSuitableProfiles} from '../mpo/utils';
import {checkPackage} from '../utils/elemental-analysis-utils';

const CREATE_NEW_PROFILE_ITEM = '+ Create New...';
const UNTITLED_PROFILE = 'Untitled Profile';

export class MpoProfileDialog {
  dataFrame: DG.DataFrame;
  mpoProfileEditor: MpoProfileEditor;
  profileInput: DG.ChoiceInput<string | null>;
  designModeInput: DG.InputBase<boolean>;
  addParetoFront: DG.InputBase<boolean>;
  addRadarInCell: DG.InputBase<boolean>;
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
  private isNewProfile = false;

  private pmpoSettingsIcon: HTMLElement;
  private methodInput: DG.ChoiceInput<string | null>;
  private pmpoSettingsContainer: HTMLElement;
  private pmpoSettingsOpened = false;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame ?? grok.shell.t;
    this.mpoProfileEditor = new MpoProfileEditor(this.dataFrame);

    this.profileInput = ui.input.choice('Profile', {
      nullable: false,
      onValueChanged: async (value) => {
        if (value === CREATE_NEW_PROFILE_ITEM) {
          await this.startNewProfile();
          return;
        }
        await this.loadProfile(value);
      },
    });

    this.designModeInput = ui.input.toggle('Design mode', {
      value: false,
      onValueChanged: (v) => this.mpoProfileEditor.setDesignMode(!!v),
    });

    this.addParetoFront = ui.input.bool('Pareto front');
    this.addRadarInCell = ui.input.bool('Add radar');

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

    this.methodInput = ui.input.choice('Method', {
      items: ['Manual', 'Probabilistic'],
      nullable: false,
      value: 'Manual',
      onValueChanged: async (value) => {
        if (value === 'Probabilistic')
          await this.createProbabilisticProfile();
        else
          this.createManualProfile();
      },
    });

    this.pmpoSettingsContainer = ui.divV([this.methodInput.root], 'chem-mpo-d-none');

    this.pmpoSettingsIcon = ui.icons.settings(() => {
      this.pmpoSettingsOpened = !this.pmpoSettingsOpened;
      this.pmpoSettingsContainer.classList.toggle('chem-mpo-d-none', !this.pmpoSettingsOpened);
    }, 'Generate from data (pMPO)');
    this.pmpoSettingsIcon.classList.add('chem-mpo-d-none');
  }

  async init(): Promise<void> {
    this.mpoProfiles = await MpoProfileManager.load();
    this.suitableProfileNames = findSuitableProfiles(this.dataFrame, this.mpoProfiles).map((p) => p.name);

    this.profileInput.items = [CREATE_NEW_PROFILE_ITEM, ...this.mpoProfiles.map((p) => p.name)];
    requestAnimationFrame(() => {
      this.addCreateNewSeparator();
      this.highlightSuitableProfiles(this.suitableProfileNames);
    });

    const defaultProfile = this.suitableProfileNames.length > 0 ?
      this.suitableProfileNames[0] :
      this.mpoProfiles[0]?.name ?? null;

    if (defaultProfile) {
      this.profileInput.value = defaultProfile;
      await this.loadProfile(defaultProfile);
    }

    this.listenForProfileChanges();
    this.mpoContextPanel = new MpoContextPanel(this.dataFrame);
  }

  private addCreateNewSeparator(): void {
    const select = this.profileInput.input as HTMLSelectElement;
    if (select.options.length < 2)
      return;
    const separator = document.createElement('hr');
    select.insertBefore(separator, select.options[1]);
  }

  private highlightSuitableProfiles(suitableFileNames: string[]): void {
    const STAR_PLACEHOLDER = '\u2003\u2002';
    const select = this.profileInput.input as HTMLSelectElement;

    for (let i = 0; i < select.options.length; i++) {
      const option = select.options[i];
      if (option.value === CREATE_NEW_PROFILE_ITEM)
        continue;

      const text = option.textContent ?? '';
      if (text.startsWith('✓') || text.startsWith('\u2003'))
        continue;

      const isApplicable = suitableFileNames.includes(option.value);
      option.textContent = `${isApplicable ? '✓' : STAR_PLACEHOLDER} ${text}`;
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
          this.mpoProfileEditor.aggregationInput.value ?? undefined,
        );
      }
    }));
  }

  private async loadProfile(profileName: string | null): Promise<void> {
    if (!profileName)
      return;

    const profileInfo = this.mpoProfiles.find((p) => p.name === profileName);
    if (!profileInfo)
      return;

    this.isNewProfile = false;
    this.pmpoSettingsOpened = false;
    this.pmpoSettingsIcon.classList.add('chem-mpo-d-none');
    this.pmpoSettingsContainer.classList.add('chem-mpo-d-none');

    const select = this.profileInput.input as HTMLSelectElement;
    if (select.options.length > 0)
      select.options[0].textContent = CREATE_NEW_PROFILE_ITEM;

    this.currentProfileFileName = profileInfo.fileName;
    this.currentProfile = structuredClone(profileInfo);
    this.designModeInput.value = false;
    this.mpoProfileEditor.setProfile(this.currentProfile);
    this.originalProfile = structuredClone(this.currentProfile);
    this.updateSaveButtonVisibility();
    this.updateOkButtonState();
  }

  private async startNewProfile(): Promise<void> {
    this.isNewProfile = true;
    this.pmpoSettingsOpened = false;
    this.pmpoSettingsContainer.classList.add('chem-mpo-d-none');
    this.methodInput.value = 'Manual';

    const select = this.profileInput.input as HTMLSelectElement;
    if (select.options.length > 0)
      select.options[0].textContent = UNTITLED_PROFILE;

    const hasBoolColumns = [...this.dataFrame.columns].some((c) => c.type === DG.COLUMN_TYPE.BOOL);
    this.pmpoSettingsIcon.classList.toggle('chem-mpo-d-none', !hasBoolColumns);

    this.createManualProfile();
  }

  private createManualProfile(): void {
    const props: {[key: string]: PropertyDesirability} = {};
    for (const col of this.dataFrame.columns.numerical)
      props[col.name] = createDefaultNumerical(1, col.min, col.max);

    this.currentProfile = {type: DESIRABILITY_PROFILE_TYPE, name: UNTITLED_PROFILE, description: '', properties: props};
    this.currentProfileFileName = null;
    this.originalProfile = null;
    this.isNewProfile = true;

    this.mpoProfileEditor.setProfile(this.currentProfile);
    this.designModeInput.value = true;
    this.saveButton.classList.remove('chem-mpo-d-none');
    this.updateOkButtonState();
  }

  private async createProbabilisticProfile(): Promise<void> {
    const tableView = grok.shell.getTableView(this.dataFrame.name);
    if (!tableView) {
      grok.shell.error('No table view found for the current dataframe');
      return;
    }

    try {
      const pMpoItems = await grok.functions.call('EDA:getPmpoAppItems', {view: tableView});
      if (!pMpoItems?.profile)
        throw new Error('pMPO is not applicable for this dataset');

      this.currentProfile = pMpoItems.profile;
      this.currentProfileFileName = null;
      this.originalProfile = null;
      this.isNewProfile = true;
      this.mpoProfileEditor.setProfile(this.currentProfile!);
      this.designModeInput.value = true;
      this.saveButton.classList.remove('chem-mpo-d-none');
      this.updateOkButtonState();
    } catch (e) {
      grok.shell.warning(`pMPO training failed: ${e instanceof Error ? e.message : e}`);
      this.methodInput.value = 'Manual';
      this.createManualProfile();
    }
  }

  private isProfileModified(): boolean {
    if (!this.currentProfile || !this.originalProfile)
      return false;
    return !deepEqual(this.currentProfile, this.originalProfile);
  }

  private updateSaveButtonVisibility(): void {
    this.saveButton.classList.toggle('chem-mpo-d-none', !this.isNewProfile && !this.isProfileModified());
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
    if (!this.currentProfile)
      return;

    const result = await MpoProfileManager.saveProfile(this.currentProfile, this.currentProfileFileName);
    if (result.saved) {
      this.currentProfileFileName = result.fileName;
      this.isNewProfile = false;
      this.pmpoSettingsOpened = false;
      this.pmpoSettingsIcon.classList.add('chem-mpo-d-none');
      this.pmpoSettingsContainer.classList.add('chem-mpo-d-none');
      await this.refreshProfilesDropdown(this.currentProfile!.name);
      this.originalProfile = structuredClone(this.currentProfile!);
      this.updateSaveButtonVisibility();
    }
  }

  private async refreshProfilesDropdown(selectName: string): Promise<void> {
    this.mpoProfiles = await MpoProfileManager.load();
    this.suitableProfileNames = findSuitableProfiles(this.dataFrame, this.mpoProfiles).map((p) => p.name);
    this.profileInput.items = [CREATE_NEW_PROFILE_ITEM, ...this.mpoProfiles.map((p) => p.name)];
    requestAnimationFrame(() => {
      this.addCreateNewSeparator();
      this.highlightSuitableProfiles(this.suitableProfileNames);
    });
    this.profileInput.value = selectName;
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
      const radarRequested = this.addRadarInCell.value;
      const profile = this.currentProfile!;
      const profileName = profile.name || 'MPO';
      const scoreColumnNames = await computeMpo(
        this.dataFrame,
        profile,
        this.mpoProfileEditor.columnMapping,
        this.mpoProfileEditor.aggregationInput.value!,
        false,
        false,
        radarRequested,
      );

      if (scoreColumnNames.length && this.addParetoFront.value)
        this.addParetoFrontViewer(scoreColumnNames);

      if (radarRequested) {
        const mapping = this.mpoProfileEditor.columnMapping;
        const desirabilityColumnNames = Object.keys(profile.properties)
          .map((prop) => `${mapping[prop] ?? prop} (${profileName} desirability)`);
        this.addRadarInCellViewer(desirabilityColumnNames);
      }
    } catch (e) {
      grok.shell.error(`Failed to run MPO calculation: ${e instanceof Error ? e.message : e}`);
    }
  }

  private addRadarInCellViewer(desirabilityColumnNames: string[]): void {
    const view = grok.shell.getTableView(this.dataFrame.name);
    if (!checkPackage('PowerGrid', 'radarCellRenderer')) {
      grok.shell.warning('PowerGrid package is not installed');
      return;
    }
    const gc = view.grid.columns.add({
      gridColumnName: `MPO (${this.currentProfile!.name})`,
      cellType: 'radar',
    });
    gc.settings = {columnNames: desirabilityColumnNames};
  }

  private getProfileControls(): HTMLElement {
    return ui.divH([this.profileInput.root, this.pmpoSettingsIcon, this.manageButton, this.saveButton],
      'chem-mpo-profile-actions');
  }

  private getMpoControls(): HTMLElement {
    return ui.divV([
      this.pmpoSettingsContainer,
      this.mpoProfileEditor.aggregationInput.root,
      this.designModeInput.root,
      this.mpoProfileEditor.root,
      this.addParetoFront.root,
      this.addRadarInCell.root,
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
    this.mpoContextPanel?.release();
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }
}
