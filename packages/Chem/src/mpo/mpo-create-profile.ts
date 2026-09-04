import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subscription} from 'rxjs';

import {
  DesirabilityProfile,
  getMpoScoreColumnName,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MPO_SCORE_CHANGED_EVENT} from '@datagrok-libraries/statistics/src/mpo/utils';

import {attachMpoProfileAi} from '../ai-tools/mpo';
import {MpoContextPanel} from './mpo-context-panel';
import {
  MpoMethod,
  MpoPathMode,
  MPO_PROFILE_DELETED_EVENT,
  updateMpoPath,
  setupMpoBreadcrumbs,
  createDefaultProfile,
  createProfileForDf,
  mergeProfileWithDf,
  isEdaPackageInstalled,
  MpoProfileInfo,
  MpoProfileRef,
  UNTITLED_PROFILE,
} from './utils';
import {mpoProfileStore} from './mpo-profile-store';
import {saveProfileInteractive} from './mpo-profile-actions';

const FIELD_DESCRIPTIONS: Record<string, string> = {
  'Method': 'Manual desirability curve editing or data-driven MPO trained from labeled data',
  'Dataset': 'Load data to preview desirability scores as you edit the profile',
  'Aggregation': 'How individual property scores combine into the final MPO score',
};

export class MpoProfileCreateView {
  static readonly PROFILE_ID_TAG = 'chem-mpo-profile-id';

  static focusOpenEditor(profileId: string): boolean {
    const isEditor = (v: DG.ViewBase | null) => v?.temp?.[MpoProfileCreateView.PROFILE_ID_TAG] === profileId;
    for (const v of grok.shell.views) {
      if (isEditor(v)) {
        grok.shell.v = v;
        return true;
      }
    }
    if (isEditor(grok.shell.preview as DG.View | null)) {
      grok.shell.windows.showBrowse = true;
      return true;
    }
    return false;
  }

  readonly view: DG.View;
  readonly showMethod: boolean;
  readonly isEditMode: boolean;

  df: DG.DataFrame | null = null;
  profile: DesirabilityProfile;
  editor: MpoProfileEditor;
  mpoContextPanel: MpoContextPanel | null = null;

  profileEditorContainer: HTMLDivElement;
  profileViewContainer!: HTMLDivElement;
  methodInput?: DG.ChoiceInput<string | null>;
  datasetInput?: DG.InputBase;
  saved: MpoProfileRef | null = null;
  saveButton: HTMLElement | null = null;
  resetButton: HTMLElement | null = null;

  private headerEl!: HTMLElement;
  private nameErrorEl!: HTMLElement;
  private descEl!: HTMLElement;
  private toolbarEl!: HTMLElement;
  private aggregationField!: HTMLElement;
  private contentEl!: HTMLDivElement;

  tableView: DG.TableView;
  private tableViewVisible: boolean = false;
  private subs: Subscription[] = [];

  private originalProfile: DesirabilityProfile;
  private profileModified = false;
  private updatingLayout = false;
  private suppressInputHandlers = false;
  private stashedManualProfile: { profile: DesirabilityProfile; modified: boolean } | null = null;

  private pMpoDockedItems: {
    statsGrid?: DG.DockNode;
    rocCurve?: DG.DockNode;
    confusionMatrix?: DG.DockNode;
    controls?: {
      form: DG.DockNode;
      saveBtn?: HTMLElement;
    };
  } | null = null;

  constructor(
    existingProfile?: DesirabilityProfile | MpoProfileInfo,
    showMethod: boolean = true,
  ) {
    this.view = DG.View.create();
    this.showMethod = showMethod;
    this.isEditMode = !!existingProfile;
    if (existingProfile && 'id' in existingProfile) {
      this.saved = existingProfile;
      this.activeView.temp[MpoProfileCreateView.PROFILE_ID_TAG] = this.saved.id;
    }

    this.profile = existingProfile ? structuredClone(existingProfile) : createDefaultProfile();
    this.originalProfile = structuredClone(this.profile);
    this.editor = new MpoProfileEditor(undefined, true);
    this.editor.setProfile(this.profile);
    this.profileEditorContainer = ui.divV([this.editor.root]);
    this.profileEditorContainer.classList.add('chem-profile-editor-container');

    this.tableView = DG.TableView.create(DG.DataFrame.create(0), false);
    this.tableView.name = this.view.name = this.displayName;
    this.dockTableView();

    updateMpoPath(
      this.activeView,
      this.isEditMode ? MpoPathMode.Edit : MpoPathMode.Create,
      this.profile.name,
    );

    this.initControls(showMethod);
    attachMpoProfileAi(this);
    this.attachLayout();
    this.listenForProfileChanges();
  }

  get activeView(): DG.View {
    return this.isEditMode ? this.view : this.tableView;
  }

  private get displayName(): string {
    return this.profile.name || UNTITLED_PROFILE;
  }

  setupBreadcrumbs(): void {
    setupMpoBreadcrumbs(this.activeView, this.displayName);
  }

  get isManualMode(): boolean {
    return !this.showMethod || this.methodInput?.value !== MpoMethod.DataDriven;
  }

  // --- Construction ---

  private dockTableView() {
    if (!this.isEditMode) {
      setTimeout(() => {
        this.tableView._onAdded();
        this.setTableViewVisible(false);
      }, 0);
    }
  }

  private initControls(showMethod: boolean) {
    if (showMethod) {
      this.methodInput = ui.input.choice('Method', {
        items: [MpoMethod.Manual, MpoMethod.DataDriven],
        value: MpoMethod.Manual,
        nullable: false,
        onValueChanged: () => this.onMethodChanged(),
      });
    }

    this.datasetInput = ui.input.table('Dataset', {
      nullable: true,
      onValueChanged: (df) => this.onDatasetChanged(df),
    });

    const field = (input: DG.InputBase) =>
      ui.divV([input.root, ui.divText(FIELD_DESCRIPTIONS[input.caption], 'chem-profile-field-desc')]);

    const controls: HTMLElement[] = [];
    if (this.methodInput)
      controls.push(field(this.methodInput));
    controls.push(field(this.datasetInput));
    this.aggregationField = field(this.editor.aggregationInput);
    controls.push(this.aggregationField);

    this.saveButton = ui.button('Save', () => this.save());
    this.resetButton = ui.button('Reset', () => this.resetProfile());
    this.nameErrorEl = ui.divText('', 'chem-profile-name-error');
    this.setModified(false);

    const editable = (el: HTMLElement, onChanged: () => void, singleLine = false) => {
      el.contentEditable = 'plaintext-only';
      el.addEventListener('input', () => {
        if (!this.updatingLayout) {
          onChanged();
          this.setModified(true);
        }
      });
      if (singleLine) {
        el.addEventListener('keydown', (e) => {
          if (e.key === 'Enter')
            e.preventDefault();
        });
      }
      return el;
    };

    this.headerEl = editable(ui.h1(this.displayName), () => this.applyProfileName(this.textOf(this.headerEl)), true);
    this.headerEl.classList.add('chem-profile-header');

    this.descEl = editable(ui.h3(this.profile.description || ''), () => {
      this.profile.description = this.textOf(this.descEl);
      this.updateDescPlaceholder();
    });
    this.descEl.classList.add('chem-profile-description');
    this.updateDescPlaceholder();

    this.toolbarEl = ui.divV([ui.divV(controls)], 'chem-profile-toolbar-wrap');

    this.contentEl = ui.divV([], 'chem-profile-content');
    this.profileViewContainer = ui.divV([this.headerEl, this.nameErrorEl, this.descEl, this.toolbarEl, this.contentEl]);
    this.profileViewContainer.classList.add('chem-profile-view');

    this.view.root.append(this.profileViewContainer);
    this.setupRibbon();
  }

  private setupRibbon(): void {
    this.activeView.setRibbonPanels([[this.saveButton!, this.resetButton!]]);
  }

  // --- Event handlers ---

  private async onMethodChanged(): Promise<void> {
    if (this.suppressInputHandlers)
      return;
    this.aggregationField.classList.toggle('chem-mpo-d-none', !this.isManualMode);

    if (this.methodInput!.value === MpoMethod.DataDriven) {
      this.stashedManualProfile = {
        profile: structuredClone(this.profile),
        modified: this.profileModified,
      };
      this.clearPreviousLayout();
      this.closeContextPanel();
      if (!this.df) {
        this.showError('Data-driven MPO requires a dataset. Please select a dataset first.');
        return;
      }
      this.tableView.dataFrame = this.df;
      await this.runDataDrivenMpo();
      return;
    }

    if (this.stashedManualProfile) {
      this.profile = this.stashedManualProfile.profile;
      this.profileModified = this.stashedManualProfile.modified;
      if (!this.profileModified)
        this.originalProfile = structuredClone(this.profile);
      this.stashedManualProfile = null;
      this.prepareManualLayout();
      await this.attachLayout();
      return;
    }

    const keepChanges = this.profileModified ? await this.showKeepChangesDialog() : false;
    if (keepChanges === null) {
      this.suppressInputHandlers = true;
      this.methodInput!.value = MpoMethod.DataDriven;
      this.suppressInputHandlers = false;
      return;
    }
    this.profile = this.resolveProfileForTransition(keepChanges);
    this.prepareManualLayout();
    await this.attachLayout();
  }

  private async onDatasetChanged(df: DG.DataFrame | null): Promise<void> {
    if (this.suppressInputHandlers)
      return;
    const keepChanges = this.showMethod && this.isManualMode && this.profileModified ?
      await this.showKeepChangesDialog() : false;
    if (keepChanges === null) {
      this.suppressInputHandlers = true;
      this.datasetInput!.value = this.df;
      this.suppressInputHandlers = false;
      return;
    }

    this.df = df;
    this.closePMpoPanels();
    const indicatorRoot = this.tableViewVisible ? this.tableView.root : this.view.root;
    this.setLoading(indicatorRoot, true, 'Switching dataset...');
    await new Promise((r) => setTimeout(r, 0));

    try {
      if (!this.df) {
        this.closeContextPanel();
        this.setTableViewVisible(false);

        if (!this.isManualMode) {
          this.clearPreviousLayout();
          this.showError('Data-driven MPO requires a dataset. Please select a dataset first.');
        } else {
          if (this.showMethod)
            this.profile = this.resolveProfileForTransition(keepChanges);
          await this.attachLayout();
        }
        return;
      }

      await grok.data.detectSemanticTypes(this.df);
      this.mpoContextPanel?.updateDataFrame(this.df);

      if (!this.isManualMode) {
        this.tableView.dataFrame = this.df;
        this.clearPreviousLayout();
        await this.runDataDrivenMpo();
        return;
      }

      if (this.showMethod)
        this.profile = this.resolveProfileForTransition(keepChanges);
      this.setTableViewVisible(false);
      await this.attachLayout();
    } finally {
      this.setLoading(indicatorRoot, false);
    }
  }

  async save(): Promise<boolean> {
    const result = await saveProfileInteractive(this.profile, this.saved);
    if (result) {
      this.saved = result;
      this.activeView.temp[MpoProfileCreateView.PROFILE_ID_TAG] = result.id;
      this.originalProfile = structuredClone(this.profile);
      this.setModified(false);
      this.tableView.name = this.view.name = this.displayName;
      this.setupBreadcrumbs();
    }
    return result != null;
  }

  private applyProfileName(newName: string): void {
    const oldName = this.profile.name;
    if (newName && oldName !== newName && this.df) {
      const oldCol = this.df.col(getMpoScoreColumnName(oldName));
      if (oldCol && !this.df.col(getMpoScoreColumnName(newName)))
        oldCol.name = getMpoScoreColumnName(newName);
    }
    this.profile.name = newName;
  }

  get isModified(): boolean {return this.profileModified;}

  setDataset(df: DG.DataFrame | null): void {
    this.datasetInput!.value = df;
  }

  async setMethod(method: MpoMethod): Promise<void> {
    if (this.methodInput!.value === method)
      return;
    this.suppressInputHandlers = true;
    this.methodInput!.value = method;
    this.suppressInputHandlers = false;
    await this.onMethodChanged();
  }

  setProfileName(name: string): void {
    this.applyProfileName(name);
    this.headerEl.innerText = this.displayName;
    this.setModified(true);
  }

  setProfileDescription(text: string): void {
    this.profile.description = text;
    this.descEl.innerText = text;
    this.updateDescPlaceholder();
    this.setModified(true);
  }

  private textOf(el: HTMLElement): string {
    return el.innerText?.trim() ?? '';
  }

  private updateDescPlaceholder(): void {
    this.descEl.classList.toggle('chem-placeholder', !this.profile.description);
  }

  private setModified(modified: boolean): void {
    this.profileModified = modified;
    this.saveButton!.classList.toggle('d4-disabled', (!modified && this.saved != null) || !this.validateName());
    this.resetButton!.classList.toggle('d4-disabled', !modified);
  }

  private validateName(): boolean {
    const name = this.profile.name?.trim() ?? '';
    const taken = mpoProfileStore.items.some((p) => p.name === name && p.id !== this.saved?.id);
    this.nameErrorEl.textContent = taken ? 'A profile with this name already exists' : '';
    this.nameErrorEl.classList.toggle('chem-mpo-d-none', !taken);
    return !taken;
  }

  async resetProfile(): Promise<void> {
    this.profile = structuredClone(this.originalProfile);
    this.setModified(false);
    this.updatingLayout = true;
    try {
      this.headerEl.innerText = this.displayName;
      this.descEl.innerText = this.profile.description;
      this.updateDescPlaceholder();
      this.editor.setProfile(this.profile);
    } finally {
      this.updatingLayout = false;
    }
    this.renderContextPanel();
  }

  // --- Layout ---

  private async attachLayout(): Promise<void> {
    this.setLoading(this.view.root, true, 'Updating layout...');
    try {
      if (this.isManualMode) {
        this.editor.design = true;
        this.editor.dataFrame = this.df ?? null as any;
      }

      this.clearPreviousLayout();
      this.updatingLayout = true;
      try {
        this.editor.setProfile(this.profile);
      } finally {
        this.updatingLayout = false;
      }

      this.setupRibbon();
      this.setModified(this.profileModified);

      if (!this.df) {
        this.contentEl.append(this.profileEditorContainer);
        return;
      }

      await this.setupGridAndContextPanel();
    } finally {
      this.setLoading(this.view.root, false);
    }
  }

  private clearPreviousLayout() {
    this.profileEditorContainer.remove();
    ui.empty(this.contentEl);
  }

  private async setupGridAndContextPanel() {
    if (!this.mpoContextPanel)
      this.mpoContextPanel = new MpoContextPanel(this.df!);

    const gridContainer = ui.div([]);
    gridContainer.classList.add('chem-data-grid-container');

    const gridPlot = this.df!.plot.grid();
    gridPlot.root.classList.add('chem-data-grid');
    gridContainer.append(gridPlot.root);

    const editorPanel = ui.box(this.profileEditorContainer);
    editorPanel.classList.add('chem-editor-panel');

    const gridPanel = ui.box(gridContainer);
    gridPanel.classList.add('chem-grid-panel');

    const split = ui.splitH([editorPanel, gridPanel], {}, true);
    split.classList.add('chem-view-split');

    this.contentEl.append(split);

    if (this.df!.currentRowIdx === -1 && this.df!.rowCount > 0)
      this.df!.currentCell = this.df!.cell(0, this.df!.columns.byIndex(0).name);

    this.renderContextPanel();
  }

  private setTableViewVisible(visible: boolean, ratio = 0.78): void {
    if (this.isEditMode)
      return;

    this.tableViewVisible = visible;
    this.tableView.grid.root.classList.toggle('chem-profile-view-loading', !visible);
    const viewNode = this.tableView.dockManager.findNode(this.view.root);

    if (visible) {
      if (!viewNode)
        this.tableView.dockManager.dock(this.view.root, DG.DOCK_TYPE.FILL, null, '');
      const gridNode = this.tableView.dockManager.findNode(this.tableView.grid.root);
      const vNode = this.tableView.dockManager.findNode(this.view.root);
      if (gridNode && vNode)
        this.tableView.dockManager.dock(this.tableView.grid.root, DG.DOCK_TYPE.DOWN, vNode, '', ratio);
    } else {
      this.tableView.dockManager.dock(this.tableView.grid.root, DG.DOCK_TYPE.FILL, null, '');
      this.tableView.dockManager.dock(this.view.root, DG.DOCK_TYPE.FILL, null, '');
    }
  }

  private renderContextPanel(): void {
    this.mpoContextPanel?.render(
      this.profile,
      this.editor.columnMapping,
      this.editor.aggregationInput.value ?? undefined,
    );
  }

  private setLoading(root: HTMLElement, loading: boolean, message?: string) {
    ui.setUpdateIndicator(root, loading, message);
    this.profileViewContainer.classList.toggle('chem-profile-view-loading', loading);
  }

  private showError(message: string) {
    this.clearPreviousLayout();
    const errorDiv = ui.divText(message);
    errorDiv.classList.add('chem-mpo-error-message');
    this.contentEl.append(errorDiv);
  }

  // --- Data-driven MPO mode ---

  private async runDataDrivenMpo(): Promise<void> {
    if (!this.df)
      return;

    this.setLoading(this.view.root, true, 'Running data-driven MPO...');

    if (!isEdaPackageInstalled()) {
      this.setLoading(this.view.root, false);
      return;
    }

    try {
      this.closePMpoPanels();

      const pMpoAppItems = await grok.functions.call('EDA:getPmpoAppItems', {view: this.tableView});
      if (!pMpoAppItems) {
        this.setTableViewVisible(false);
        this.showError('Data-driven MPO is not applicable for this dataset.');
        return;
      }

      this.setTableViewVisible(true);

      const dockMng = this.tableView.dockManager;
      const gridNode = dockMng.findNode(this.tableView.grid.root);
      if (!gridNode)
        throw new Error('Failed to train data-driven MPO: missing a grid in the table view.');

      const controlsNode = dockMng.dock(pMpoAppItems.controls.form, DG.DOCK_TYPE.LEFT, gridNode, undefined, 0.1);
      const statGridNode = dockMng.dock(pMpoAppItems.statsGrid, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);
      const rocNode = dockMng.dock(pMpoAppItems.rocCurve, DG.DOCK_TYPE.RIGHT, statGridNode, undefined, 0.3);
      dockMng.dock(this.tableView.grid.root, DG.DOCK_TYPE.DOWN, rocNode, 'Training Data', 0.5);
      const gridTabNode = dockMng.findNode(this.tableView.grid.root);
      const confusionNode = dockMng.dock(pMpoAppItems.confusionMatrix, DG.DOCK_TYPE.FILL, gridTabNode);

      if (pMpoAppItems.controls.saveBtn)
        this.tableView.setRibbonPanels([[pMpoAppItems.controls.saveBtn]]);

      this.pMpoDockedItems = {
        statsGrid: statGridNode,
        rocCurve: rocNode,
        confusionMatrix: confusionNode,
        controls: {
          form: controlsNode,
          saveBtn: pMpoAppItems.controls.saveBtn,
        },
      };
    } finally {
      this.setLoading(this.view.root, false);
    }
  }

  private closePMpoPanels(): void {
    if (!this.pMpoDockedItems)
      return;

    const dockMng = this.tableView.dockManager;
    const {statsGrid, rocCurve, confusionMatrix, controls} = this.pMpoDockedItems;

    if (statsGrid)
      dockMng.close(statsGrid);
    if (rocCurve)
      dockMng.close(rocCurve);
    if (confusionMatrix)
      dockMng.close(confusionMatrix);

    if (controls) {
      if (controls.form)
        dockMng.close(controls.form);
      if (controls.saveBtn)
        controls.saveBtn.remove();
    }

    this.pMpoDockedItems = null;
  }

  private prepareManualLayout(): void {
    this.clearPreviousLayout();
    this.closePMpoPanels();
    this.setTableViewVisible(false);
  }

  // --- Profile state ---

  private resolveProfileForTransition(keepChanges: boolean): DesirabilityProfile {
    if (!keepChanges) {
      this.profileModified = false;
      this.stashedManualProfile = null;
      const profile = this.df ? createProfileForDf(this.df) : createDefaultProfile();
      profile.name = this.profile.name;
      profile.description = this.profile.description;
      this.originalProfile = structuredClone(profile);
      return profile;
    }

    if (this.df)
      return mergeProfileWithDf(this.profile, this.df);

    return this.profile;
  }

  private showKeepChangesDialog(): Promise<boolean | null> {
    return new Promise((resolve) => {
      let resolved = false;
      const safeResolve = (value: boolean | null) => {
        if (resolved)
          return;
        resolved = true;
        dlg.close();
        resolve(value);
      };

      const dlg = ui.dialog('Keep profile changes?')
        .addButton('Keep', () => safeResolve(true))
        .addButton('Discard', () => safeResolve(false))
        .onCancel(() => safeResolve(null))
        .show({center: true});
    });
  }

  // --- Lifecycle ---

  private listenForProfileChanges() {
    const isOwnView = (v: DG.View | null) => v && (v.id === this.view.id || v.id === this.tableView.id);

    this.subs.push(grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      if (!this.updatingLayout)
        this.setModified(true);
      if (!this.df || !this.profile || !this.mpoContextPanel || !isOwnView(grok.shell.v as DG.View))
        return;
      this.renderContextPanel();
    }));

    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_DELETED_EVENT).subscribe((data) => {
      if (data?.id === this.saved?.id)
        this.closeView();
    }));

    this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
      if (!isOwnView(grok.shell.v as DG.View))
        this.mpoContextPanel?.release();
    }));
    this.subs.push(grok.events.onViewRemoving.subscribe((data: DG.EventData) => {
      if (isOwnView(data.args?.view))
        this.detach();
    }));
  }

  private closeContextPanel(): void {
    this.mpoContextPanel?.release();
    this.mpoContextPanel = null;
  }

  private detach(): void {
    this.closeContextPanel();
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }

  private closeView(): void {
    this.detach();
    this.activeView.close();
  }
}
