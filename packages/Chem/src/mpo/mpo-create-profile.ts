import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';

import {
  DesirabilityProfile,
  PropertyDesirability,
  WEIGHTED_AGGREGATIONS_LIST,
  WeightedAggregation,
} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';

import {MpoContextPanel} from './mpo-context-panel';
import {MPO_TEMPLATE_PATH} from './utils';

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

  tableView: DG.TableView | null = null;

  constructor(existingProfile?: DesirabilityProfile, showMethod: boolean = true, fileName?: string) {
    this.view = DG.View.create();
    this.showMethod = showMethod;
    this.isEditMode = !!existingProfile;
    this.fileName = fileName;

    this.profile = existingProfile ?? this.createDefaultProfile();
    this.editor = new MpoProfileEditor(undefined, true);
    this.editor.setProfile(this.profile);
    this.profileEditorContainer = ui.divV([this.editor.root]);
    this.profileEditorContainer.classList.add('chem-profile-editor-container');

    if (!this.isEditMode)
      this.initTableView();

    this.view.name = this.isEditMode ? 'Edit MPO Profile' : 'Create MPO Profile';
    this.updateViewPaths();

    this.initControls(showMethod);
    this.initSaveButton();
    this.attachLayout();
    this.listenForProfileChanges();
  }

  private updateViewPaths() {
    const path = this.isEditMode && this.profile.name ? this.profile.name : 'Create MPO Profile';
    // this.view.path = this.tableView!.path = /*createPath(path);*/path;
  }

  private initTableView(): void {
    const emptyDf = DG.DataFrame.create(1);
    this.tableView = DG.TableView.create(emptyDf, false);

    setTimeout(() => {
      this.tableView!._onAdded();
      this.tableView!.grid.root.style.visibility = 'hidden';
      this.tableView!.dockManager.dock(
        this.view.root,
        DG.DOCK_TYPE.TOP,
        null,
        '',
        0.99,
      );
    }, 0);
  }

  private setTableViewVisible(visible: boolean, ratio = 0.7): void {
    if (this.isEditMode || !this.tableView)
      return;

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
      const fileNameInput = ui.input.string('File name', {value: this.fileName ?? 'mpo-profile.json'});
      const nameInput = ui.input.string('Name', {value: this.profile.name ?? ''});
      const descInput = ui.input.string('Description', {value: this.profile.description ?? ''});

      ui.dialog({title: 'Save MPO Profile'})
        .add(ui.divV([fileNameInput, nameInput, descInput]))
        .onOK(() => {
          this.profile.name = nameInput.value || '';
          this.profile.description = descInput.value || '';

          const fileName = this.fileName ?? fileNameInput.value!.trim();
          grok.dapi.files.writeAsText(`${MPO_TEMPLATE_PATH}/${fileName}`, JSON.stringify(this.profile));
          this.saveButton!.style.display = 'none';
          grok.shell.info(`Profile "${this.profile.name}" saved.`);
        })
        .show();
    });

    this.saveButton.style.display = 'none';
    this.view.setRibbonPanels([[this.saveButton]]);
    this.tableView?.setRibbonPanels([[this.saveButton]]);
  }

  private initControls(showMethod: boolean) {
    const controls: DG.InputBase[] = [];

    if (showMethod) {
      this.methodInput = ui.input.choice('Method', {
        items: [METHOD_MANUAL, METHOD_PROBABILISTIC],
        value: METHOD_MANUAL,
        nullable: false,
        onValueChanged: async () => {
          if (this.methodInput!.value === METHOD_PROBABILISTIC) {
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
      onValueChanged: (df) => {
        this.df = df;

        if (this.tableView && df)
          this.tableView.dataFrame = df;

        this.attachLayout();
      },
    });
    controls.push(datasetInput);

    this.aggregationInput = ui.input.choice('Aggregation', {
      items: WEIGHTED_AGGREGATIONS_LIST,
      nullable: false,
      onValueChanged: () => grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {}),
    });
    this.aggregationInput.enabled = false;
    controls.push(this.aggregationInput);

    const headerDiv = ui.divV([ui.h1(this.view.name), ui.form(controls)]);
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

  private async runProbabilisticMpo(): Promise<void> {
    if (!this.df || !this.tableView)
      return;

    ui.setUpdateIndicator(this.view.root, true, 'Running probabilistic MPO...');

    this.setTableViewVisible(true);

    const prev = grok.shell.v;
    grok.shell.v = this.tableView;

    try {
      await grok.functions.call('Eda:trainPmpo', {});
    } finally {
      grok.shell.v = prev;

      ui.setUpdateIndicator(this.view.root, false);
    }
  }

  private clearPreviousLayout() {
    while (this.profileViewContainer.children.length > 1)
      this.profileViewContainer.removeChild(this.profileViewContainer.lastChild!);
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
        'Property 1': {weight: 1, min: 0, max: 1, line: []},
        'Property 2': {weight: 1, min: 0, max: 1, line: []},
        'Property 3': {weight: 1, min: 0, max: 1, line: []},
      },
    };
  }

  private createProfileForDf(): DesirabilityProfile {
    const props: {[key: string]: PropertyDesirability} = {};
    for (const col of this.df!.columns.numerical)
      props[col.name] = {weight: 1, min: col.min, max: col.max, line: []};

    return {name: '', description: '', properties: props};
  }

  private async setupGridAndContextPanel() {
    if (!this.mpoContextPanel)
      this.mpoContextPanel = new MpoContextPanel(this.df!);

    await grok.data.detectSemanticTypes(this.df!);

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

  private async listenForProfileChanges() {
    grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      this.saveButton!.style.display = 'initial';
      if (!this.df || !this.profile || !this.mpoContextPanel) return;

      const agg = this.aggregationInput?.value ?? 'Average';
      await this.mpoContextPanel.render(this.profile, this.editor.columnMapping, agg);
    });
  }
}
