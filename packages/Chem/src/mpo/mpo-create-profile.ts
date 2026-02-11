import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {Subscription} from 'rxjs';

import {
  DesirabilityProfile,
  PropertyDesirability,
  WEIGHTED_AGGREGATIONS_LIST,
  WeightedAggregation,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {FormulaEditor} from './formula-editor';
import {MpoContextPanel} from './mpo-context-panel';
import {MpoPathMode, MPO_PROFILE_DELETED_EVENT, updateMpoPath} from './utils';
import {MpoProfileManager} from './mpo-profile-manager';

const METHOD_MANUAL = 'Manual';
const METHOD_PROBABILISTIC = 'Probabilistic';

export class MpoProfileCreateView {
  readonly view: DG.View;
  readonly showMethod: boolean;
  readonly isEditMode: boolean;

  df: DG.DataFrame | null = null;
  profile: DesirabilityProfile;
  editor: MpoProfileEditor;
  mpoContextPanel: MpoContextPanel | null = null;

  profileEditorContainer: HTMLDivElement;
  profileViewContainer: HTMLDivElement = ui.div();
  methodInput?: DG.ChoiceInput<string | null>;
  aggregationInput?: DG.ChoiceInput<WeightedAggregation | null>;
  fileName?: string | null = null;
  saveButton: HTMLElement | null = null;

  tableView: DG.TableView;
  private tableViewVisible: boolean = false;
  private subs: Subscription[] = [];
  private formulaEditor = new FormulaEditor();

  private pMpoDockedItems: {
    statsGrid?: DG.DockNode;
    rocCurve?: DG.DockNode;
    confusionMatrix?: DG.DockNode;
    controls?: {
      form: DG.DockNode;
      saveBtn: HTMLElement;
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

    this.profile = existingProfile ?? this.createDefaultProfile();
    this.editor = new MpoProfileEditor(undefined, true);
    this.editor.setProfile(this.profile);
    this.profileEditorContainer = ui.divV([this.editor.root]);
    this.profileEditorContainer.classList.add('chem-profile-editor-container');

    this.tableView = DG.TableView.create(DG.DataFrame.create(0), false);
    this.tableView.name = this.view.name = this.isEditMode ?
      `Edit ${this.profile.name ?? 'MPO'} Profile` :
      'Create MPO Profile';
    this.dockTableView();

    updateMpoPath(
      this.isEditMode ? this.view : this.tableView,
      this.isEditMode ? MpoPathMode.Edit : MpoPathMode.Create,
      this.profile.name,
    );

    this.initControls(showMethod);
    this.initSaveButton();
    this.attachLayout();
    this.listenForProfileChanges();
  }

  private dockTableView() {
    if (!this.isEditMode) {
      setTimeout(() => {
        this.tableView._onAdded();
        this.tableView.grid.root.style.visibility = 'hidden';
        this.tableView.dockManager.dock(this.view.root, DG.DOCK_TYPE.TOP, null, '', 0.99);
      }, 0);
    }
  }

  private setTableViewVisible(visible: boolean, ratio = 0.25): void {
    if (this.isEditMode)
      return;

    this.tableViewVisible = visible;
    this.tableView.grid.root.style.visibility = visible ? 'visible' : 'hidden';
    this.tableView.dockManager.dock(
      this.view.root,
      DG.DOCK_TYPE.TOP,
      null,
      '',
      visible ? ratio : 0.99,
    );
  }

  private initSaveButton() {
    this.saveButton = ui.bigButton('Save', async () => {
      await MpoProfileManager.ensureLoaded();
      const existingFileNames = MpoProfileManager.existingFileNames;
      const fileNameInput = ui.input.string('File name', {value: this.fileName ?? 'mpo-profile.json', nullable: false});
      const normalizeFileName = (name: string) =>
        name.endsWith('.json') ? name : `${name}.json`;

      fileNameInput.addValidator((value) => {
        const normalized = normalizeFileName(value);
        const canOverwrite = this.isEditMode && normalized === this.fileName;
        if (!canOverwrite && existingFileNames?.has(normalized))
          return `File "${normalized}" already exists`;
        return null;
      });

      const nameInput = ui.input.string('Name', {value: this.profile.name ?? ''});
      const descInput = ui.input.string('Description', {value: this.profile.description ?? ''});

      const dlg = ui.dialog({title: 'Save MPO Profile'})
        .add(ui.divV([fileNameInput, nameInput, descInput]))
        .onOK(async () => {
          this.profile.name = nameInput.value || '';
          this.profile.description = descInput.value || '';

          const fileName = normalizeFileName(fileNameInput.value!.trim());
          const saved = await MpoProfileManager.save(this.profile, fileName);
          if (saved) {
            this.fileName = fileName;
            this.saveButton!.style.display = 'none';
          }
        })
        .show();

      const okButton = dlg.getButton('OK');
      fileNameInput.onInput.subscribe(() => okButton.disabled = !fileNameInput.validate());
    });

    this.saveButton.style.display = 'none';
    (this.isEditMode ? this.view : this.tableView).setRibbonPanels([[this.saveButton]]);
  }

  private initControls(showMethod: boolean) {
    const controls: DG.InputBase[] = [];

    if (showMethod) {
      this.methodInput = ui.input.choice('Method', {
        items: [METHOD_MANUAL, METHOD_PROBABILISTIC],
        value: METHOD_MANUAL,
        nullable: false,
        onValueChanged: async () => {
          this.clearPreviousLayout();
          if (this.methodInput!.value === METHOD_PROBABILISTIC) {
            this.closeContextPanel();
            if (!this.df) {
              this.showError('Probabilistic MPO requires a dataset. Please select a dataset first.');
              return;
            }
            this.tableView.dataFrame = this.df;
            await this.runProbabilisticMpo();
            return;
          }
          this.setTableViewVisible(false);
          this.attachLayout();
        },
      });
      controls.push(this.methodInput);
    }

    const datasetInput = ui.input.table('Dataset', {
      nullable: true,
      onValueChanged: async (df) => {
        this.closePMpoPanels();
        const indicatorRoot = this.tableViewVisible ? this.tableView.root : this.view.root;
        ui.setUpdateIndicator(indicatorRoot, true, 'Switching dataset...');
        await new Promise((r) => setTimeout(r, 0));

        try {
          this.df = df;

          if (!this.df) {
            this.closeContextPanel();
            this.setTableViewVisible(false);

            if (this.methodInput?.value === METHOD_PROBABILISTIC) {
              this.clearPreviousLayout();
              this.showError('Probabilistic MPO requires a dataset. Please select a dataset first.');
            } else
              this.attachLayout();
            return;
          }

          await grok.data.detectSemanticTypes(this.df);
          this.mpoContextPanel?.updateDataFrame(this.df);

          if (this.methodInput?.value === METHOD_PROBABILISTIC) {
            this.tableView.dataFrame = this.df;
            this.clearPreviousLayout();
            await this.runProbabilisticMpo();
            return;
          }

          this.setTableViewVisible(false);
          this.attachLayout();
        } finally {
          ui.setUpdateIndicator(indicatorRoot, false);
        }
      },
    });
    controls.push(datasetInput);

    this.aggregationInput = ui.input.choice('Aggregation', {
      items: WEIGHTED_AGGREGATIONS_LIST,
      nullable: false,
      onValueChanged: () => {
        const isFormula = this.aggregationInput!.value === 'Formula';
        this.formulaEditor.visible = isFormula;
        if (isFormula && this.profile && this.df)
          this.formulaEditor.build(this.profile, this.df);
        grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {});
      },
    });
    this.aggregationInput.enabled = false;
    controls.push(this.aggregationInput);

    const headerDiv = ui.divV([ui.h1(this.view.name), ui.form(controls), this.formulaEditor.root]);
    headerDiv.classList.add('chem-profile-header');

    this.profileViewContainer = ui.divV([headerDiv]);
    this.profileViewContainer.classList.add('chem-profile-view');

    this.view.root.append(this.profileViewContainer);
  }

  private async attachLayout() {
    ui.setUpdateIndicator(this.view.root, true, 'Updating layout...');
    try {
      this.clearPreviousLayout();
      this.ensureProfile();

      if (this.methodInput?.value === METHOD_MANUAL || !this.showMethod) {
        this.editor.design = true;
        this.editor.dataFrame = this.df ?? null as any;
        if (this.aggregationInput)
          this.aggregationInput.enabled = !!this.df;

        if (this.showMethod) {
          this.profile = this.df ?
            this.createProfileForDf() :
            this.createDefaultProfile();
        }
      }

      this.editor.setProfile(this.profile);

      if (!this.df) {
        this.profileViewContainer.append(this.profileEditorContainer);
        return;
      }

      await this.setupGridAndContextPanel();
    } finally {
      ui.setUpdateIndicator(this.view.root, false);
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

  private async runProbabilisticMpo(): Promise<void> {
    if (!this.df)
      return;

    ui.setUpdateIndicator(this.view.root, true, 'Running probabilistic MPO...');

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
      ui.setUpdateIndicator(this.view.root, false);
    }
  }

  private clearPreviousLayout() {
    const header = this.profileViewContainer.children[0];
    ui.empty(this.profileViewContainer);
    if (header)
      this.profileViewContainer.append(header);
  }

  private showError(message: string) {
    this.clearPreviousLayout();
    const errorDiv = ui.divText(message);
    errorDiv.classList.add('chem-mpo-error-message');
    this.profileViewContainer.append(errorDiv);
  }

  private ensureProfile() {
    if (!this.profile)
      this.profile = this.createDefaultProfile();
  }

  private createDefaultProfile(): DesirabilityProfile {
    return {
      name: '',
      description: '',
      properties: {
        'Property 1': {functionType: 'numerical', weight: 1, min: 0, max: 1, line: []},
        'Property 2': {functionType: 'numerical', weight: 1, min: 0, max: 1, line: []},
        'Property 3': {functionType: 'numerical', weight: 1, min: 0, max: 1, line: []},
      },
    };
  }

  private createProfileForDf(): DesirabilityProfile {
    const props: {[key: string]: PropertyDesirability} = {};
    for (const col of this.df!.columns.numerical)
      props[col.name] = {functionType: 'numerical', weight: 1, min: col.min, max: col.max, line: []};

    return {name: '', description: '', properties: props};
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

    await this.mpoContextPanel!.render(this.profile, this.editor.columnMapping, 'Average');
  }

  private listenForProfileChanges() {
    this.subs.push(grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      this.saveButton!.style.display = 'initial';
      if (!this.df || !this.profile || !this.mpoContextPanel) return;

      const agg = this.aggregationInput?.value ?? 'Average';
      await this.mpoContextPanel.render(this.profile, this.editor.columnMapping, agg, this.formulaEditor.getExpression());
    }));

    this.subs.push(grok.events.onCustomEvent(MPO_PROFILE_DELETED_EVENT).subscribe((data) => {
      if (data?.fileName === this.fileName)
        this.closeView();
    }));

    const viewHandler = (eventData: DG.EventData) => {
      const eventView = eventData.args?.view;
      if (eventView && (eventView.id === this.view.id || eventView.id === this.tableView.id))
        this.detach();
    };

    this.subs.push(grok.events.onViewChanging.subscribe(viewHandler));
    this.subs.push(grok.events.onViewRemoving.subscribe(viewHandler));
  }

  private closeContextPanel(): void {
    this.mpoContextPanel?.close();
    this.mpoContextPanel = null;
  }

  private detach(): void {
    this.closeContextPanel();
    this.formulaEditor.disconnect();
    this.subs.forEach((sub) => sub.unsubscribe());
    this.subs = [];
  }

  private closeView(): void {
    this.detach();
    (this.isEditMode ? this.view : this.tableView).close();
  }
}
