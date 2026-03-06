import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subscription} from 'rxjs';

import {
  DesirabilityProfile,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MPO_SCORE_CHANGED_EVENT} from '@datagrok-libraries/statistics/src/mpo/utils';
import {discoverComputeFunctions} from '@datagrok-libraries/statistics/src/compute-functions/discovery';
import {chemFunctionsDialog} from '@datagrok-libraries/statistics/src/compute-functions/dialog';
import {IComputeDialogResult, IFunctionArgs} from '@datagrok-libraries/statistics/src/compute-functions/types';

import {MpoContextPanel} from './mpo-context-panel';
import {
  MpoPathMode,
  MPO_PROFILE_DELETED_EVENT,
  updateMpoPath,
  createDefaultProfile,
  createProfileForDf,
  mergeProfileWithDf,
  templateFromCallString,
} from './utils';
import {MpoProfileManager} from './mpo-profile-manager';

const METHOD_MANUAL = 'Manual';
const METHOD_PROBABILISTIC = 'Probabilistic';

function buildCallString(funcKey: string, args: IFunctionArgs): string {
  const vals = Object.values(args);
  return vals.length ? `${funcKey}(${vals.map((v) => JSON.stringify(v)).join(', ')})` : funcKey;
}

export class MpoProfileCreateView {
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
  fileName?: string | null = null;
  saveButton: HTMLElement | null = null;

  private headerEl!: HTMLElement;
  private toolbarEl!: HTMLElement;

  tableView: DG.TableView;
  private tableViewVisible: boolean = false;
  private subs: Subscription[] = [];

  private profileModified = false;
  private updatingLayout = false;
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
    existingProfile?: DesirabilityProfile,
    showMethod: boolean = true,
    fileName?: string,
  ) {
    this.view = DG.View.create();
    this.showMethod = showMethod;
    this.isEditMode = !!existingProfile;
    this.fileName = fileName;

    this.profile = existingProfile ?? createDefaultProfile();
    this.editor = new MpoProfileEditor(undefined, true);
    this.editor.setProfile(this.profile);
    this.editor.onComputeRequested.subscribe((propName) => this.openComputeDialog(propName));
    this.profileEditorContainer = ui.divV([this.editor.root]);
    this.profileEditorContainer.classList.add('chem-profile-editor-container');

    this.tableView = DG.TableView.create(DG.DataFrame.create(0), false);
    const tabName = this.isEditMode ?
      `Edit ${this.profile.name || 'MPO'}` :
      'Create MPO';
    this.tableView.name = this.view.name = tabName;
    this.dockTableView();

    updateMpoPath(
      this.activeView,
      this.isEditMode ? MpoPathMode.Edit : MpoPathMode.Create,
      this.profile.name,
    );

    this.initControls(showMethod);
    this.attachLayout();
    this.listenForProfileChanges();
  }

  private get activeView(): DG.View {
    return this.isEditMode ? this.view : this.tableView;
  }

  private get isManualMode(): boolean {
    return !this.showMethod || this.methodInput?.value !== METHOD_PROBABILISTIC;
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
        items: [METHOD_MANUAL, METHOD_PROBABILISTIC],
        value: METHOD_MANUAL,
        nullable: false,
        onValueChanged: () => this.onMethodChanged(),
      });
      this.methodInput.addPostfix('Manual desirability curve editing or probabilistic MPO trained from labeled data');
    }

    this.datasetInput = ui.input.table('Dataset', {
      nullable: true,
      onValueChanged: (df) => this.onDatasetChanged(df),
    });
    this.datasetInput.addOptions(
      ui.divText('Load data to preview desirability scores as you edit the profile', 'ui-input-description'));

    this.editor.aggregationInput.addPostfix('How individual property scores combine into the final MPO score');

    const controls: DG.InputBase[] = [];
    if (this.methodInput)
      controls.push(this.methodInput);
    controls.push(this.datasetInput);
    controls.push(this.editor.aggregationInput);

    this.saveButton = ui.button('Save', () => this.showSaveDialog());
    this.saveButton.classList.add('d4-disabled');

    this.headerEl = ui.h1(this.isEditMode ? `Edit ${this.profile.name || 'MPO'}` : 'Create MPO Profile');
    this.headerEl.classList.add('chem-profile-header');

    this.toolbarEl = ui.divV([ui.form(controls), this.saveButton], 'chem-profile-toolbar-wrap');

    this.profileViewContainer = ui.divV([this.headerEl, this.toolbarEl]);
    this.profileViewContainer.classList.add('chem-profile-view');

    this.view.root.append(this.profileViewContainer);
  }

  // --- Event handlers ---

  private async onMethodChanged(): Promise<void> {
    this.editor.aggregationInput.root.classList.toggle('chem-mpo-d-none', !this.isManualMode);

    if (this.methodInput!.value === METHOD_PROBABILISTIC) {
      this.stashedManualProfile = {
        profile: structuredClone(this.profile),
        modified: this.profileModified,
      };
      this.clearPreviousLayout();
      this.closeContextPanel();
      if (!this.df) {
        this.showError('Probabilistic MPO requires a dataset. Please select a dataset first.');
        return;
      }
      this.tableView.dataFrame = this.df;
      await this.runProbabilisticMpo();
      return;
    }

    if (this.stashedManualProfile) {
      this.profile = this.stashedManualProfile.profile;
      this.profileModified = this.stashedManualProfile.modified;
      this.stashedManualProfile = null;
      this.prepareManualLayout();
      await this.attachLayout();
      return;
    }

    const keepChanges = this.profileModified ? await this.showKeepChangesDialog() : false;
    this.profile = this.resolveProfileForTransition(keepChanges);
    this.prepareManualLayout();
    await this.attachLayout();
  }

  private async onDatasetChanged(df: DG.DataFrame | null): Promise<void> {
    this.df = df;
    const keepChanges = this.showMethod && this.isManualMode && this.profileModified ?
      await this.showKeepChangesDialog() : false;

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
          this.showError('Probabilistic MPO requires a dataset. Please select a dataset first.');
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
        await this.runProbabilisticMpo();
        return;
      }

      if (this.showMethod) this.profile = this.resolveProfileForTransition(keepChanges);
      this.setTableViewVisible(false);
      await this.attachLayout();
    } finally {
      this.setLoading(indicatorRoot, false);
    }
  }

  private async showSaveDialog(): Promise<void> {
    await MpoProfileManager.ensureLoaded();

    const nameInput = ui.input.string('Name', {value: this.profile.name ?? '', nullable: false});
    const descInput = ui.input.textArea('Description', {value: this.profile.description ?? ''});

    const dlg = ui.dialog({title: 'Save MPO Profile'})
      .add(ui.divV([nameInput, descInput]))
      .onOK(async () => {
        this.profile.name = nameInput.value || '';
        this.profile.description = descInput.value || '';

        const fileName = this.isEditMode ?
          this.fileName! :
          MpoProfileManager.generateFileName(nameInput.value!.trim());
        const saved = await MpoProfileManager.save(this.profile, fileName);
        if (saved) {
          this.fileName = fileName;
          this.saveButton!.classList.add('d4-disabled');
          this.profileModified = false;
        }
      })
      .show();

    const okButton = dlg.getButton('OK');
    okButton.disabled = !nameInput.validate();
    nameInput.onInput.subscribe(() => okButton.disabled = !nameInput.validate());
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

      if (this.profileModified && this.saveButton)
        this.saveButton.classList.remove('d4-disabled');

      if (!this.df) {
        this.profileViewContainer.append(this.profileEditorContainer);
        return;
      }

      await this.setupGridAndContextPanel();
    } finally {
      this.setLoading(this.view.root, false);
    }
  }

  private clearPreviousLayout() {
    ui.empty(this.profileViewContainer);
    this.profileViewContainer.append(this.headerEl, this.toolbarEl);
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

    this.profileViewContainer.append(split);

    if (this.df!.currentRowIdx === -1 && this.df!.rowCount > 0)
      this.df!.currentCell = this.df!.cell(0, this.df!.columns.byIndex(0).name);

    await this.renderContextPanel();
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
    }
    else {
      this.tableView.dockManager.dock(this.tableView.grid.root, DG.DOCK_TYPE.FILL, null, '');
      this.tableView.dockManager.dock(this.view.root, DG.DOCK_TYPE.FILL, null, '');
    }
  }

  private async renderContextPanel(): Promise<void> {
    await this.mpoContextPanel?.render(
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
    this.profileViewContainer.append(errorDiv);
  }

  // --- pMPO mode ---

  private async runProbabilisticMpo(): Promise<void> {
    if (!this.df)
      return;

    this.setLoading(this.view.root, true, 'Running probabilistic MPO...');

    try {
      this.closePMpoPanels();

      const pMpoAppItems = await grok.functions.call('EDA:getPmpoAppItems', {view: this.tableView});
      if (!pMpoAppItems) {
        this.setTableViewVisible(false);
        this.showError('pMPO is not applicable for this dataset.');
        return;
      }

      this.setTableViewVisible(true);

      const dockMng = this.tableView.dockManager;
      const gridNode = dockMng.findNode(this.tableView.grid.root);
      if (!gridNode)
        throw new Error('Failed to train pMPO: missing a grid in the table view.');

      const controlsNode = dockMng.dock(pMpoAppItems.controls.form, DG.DOCK_TYPE.LEFT, gridNode, undefined, 0.1);
      const statGridNode = dockMng.dock(pMpoAppItems.statsGrid, DG.DOCK_TYPE.DOWN, gridNode, undefined, 0.5);
      const rocNode = dockMng.dock(pMpoAppItems.rocCurve, DG.DOCK_TYPE.RIGHT, statGridNode, undefined, 0.3);
      const confusionNode = dockMng.dock(pMpoAppItems.confusionMatrix, DG.DOCK_TYPE.RIGHT, rocNode, undefined, 0.2);

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
      return this.df ? createProfileForDf(this.df) : createDefaultProfile();
    }

    if (this.df)
      return mergeProfileWithDf(this.profile, this.df);

    return this.profile;
  }

  private showKeepChangesDialog(): Promise<boolean> {
    return new Promise((resolve) => {
      let resolved = false;
      const safeResolve = (value: boolean) => {
        if (resolved)
          return;
        resolved = true;
        dlg.close();
        resolve(value);
      };

      const dlg = ui.dialog('Keep profile changes?')
        .addButton('Keep', () => safeResolve(true))
        .addButton('Discard', () => safeResolve(false))
        .onCancel(() => safeResolve(false))
        .show({center: true});
    });
  }

  private async openComputeDialog(propName: string): Promise<void> {
    const prop = this.profile.properties[propName];
    const template = prop?.function ? templateFromCallString(prop.function) : undefined;

    const computeFunctions = discoverComputeFunctions('HitTriageFunction');
    let dialogResult: IComputeDialogResult | null = null;

    const {root, okProxy} = await chemFunctionsDialog(
      computeFunctions, (res) => { dialogResult = res; }, () => {},
      template, false,
    );

    ui.dialog('Compute Properties')
      .add(root)
      .onOK(() => {
        okProxy();
        if (!dialogResult)
          return;
        const entry = Object.entries(dialogResult.externals)[0];
        const prop = this.profile.properties[propName];
        if (entry && prop)
          prop.function = buildCallString(entry[0], entry[1]);
      })
      .show({resizable: true});
  }

  // --- Lifecycle ---

  private listenForProfileChanges() {
    const isOwnView = (v: DG.View | null) => v && (v.id === this.view.id || v.id === this.tableView.id);

    this.subs.push(grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      if (!this.updatingLayout) {
        this.saveButton!.classList.remove('d4-disabled');
        this.profileModified = true;
      }
      if (!this.df || !this.profile || !this.mpoContextPanel)
        return;
      await this.renderContextPanel();
    }));

    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_DELETED_EVENT).subscribe((data) => {
      if (data?.fileName === this.fileName)
        this.closeView();
    }));

    this.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
      if (!isOwnView(grok.shell.v as DG.View))
        this.closeContextPanel();
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
