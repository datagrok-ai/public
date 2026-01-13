/* eslint-disable guard-for-in */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DesirabilityProfile, PropertyDesirability, WeightedAggregation, WEIGHTED_AGGREGATIONS} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MPO_SCORE_CHANGED_EVENT, MpoProfileEditor} from '@datagrok-libraries/statistics/src/mpo/mpo-profile-editor';
import {MpoScoreViewer} from '../apps/mpo-scores';
import {MpoContextPanel} from '../apps/mpo-context-panel';

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

  // private mpoHistogramHost: HTMLDivElement = ui.div();
  private mpoContextPanel?: MpoContextPanel;
  // private mpoHistogram?: DG.Viewer;
  // private mpoBestScoreViewer?: MpoScoreViewer;
  // private mpoWorstScoreViewer?: MpoScoreViewer;

  constructor(dataFrame?: DG.DataFrame) {
    this.dataFrame = dataFrame ?? grok.shell.t;
    this.mpoProfileEditor = new MpoProfileEditor(this.dataFrame);

    this.aggregationInput = ui.input.choice('Aggregation', {
      items: [...WEIGHTED_AGGREGATIONS], nullable: false,
      onValueChanged: (v) => grok.events.fireCustomEvent(MPO_SCORE_CHANGED_EVENT, {}),
    });

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

    this.listenForProfileChanges();
    // this.mpoContextPanel = this.createMpoContextPanel();
    this.mpoContextPanel = new MpoContextPanel(this.dataFrame);
  }

  private async calculateMpoScores(): Promise<string[]> {
    if (!this.currentProfile)
      return [];

    const mapping = this.mpoProfileEditor.columnMapping;

    const mappedProperties: Record<string, PropertyDesirability> = {};
    for (const [propName, prop] of Object.entries(this.currentProfile.properties)) {
      const columnName = mapping[propName] ?? propName;
      mappedProperties[columnName] = prop;
    }

    const [func] = await DG.Func.find({name: 'mpoTransformFunction'});
    const funcCall = await func.prepare({
      df: this.dataFrame,
      profileName: this.currentProfile.name ?? 'MPO',
      currentProperties: mappedProperties,
      aggregation: this.aggregationInput.value,
    }).call(undefined, undefined, {processed: false});

    return funcCall.getOutputParamValue();
  }

  private createMpoContextPanel(): DG.Accordion {
    const acc = ui.accordion();

    const icon = ui.element('i');
    icon.className = 'grok-icon svg-icon svg-histogram';

    const statusLabel = ui.label('MPO not calculated yet');
    acc.addTitle(ui.span([icon, statusLabel]));

    grok.shell.o = acc.root;
    return acc;
  }


  // private createScoresPanes(mpoColName: string): void {
  //   if (!this.mpoBestScoreViewer) {
  //     this.mpoBestScoreViewer = new MpoScoreViewer(this.dataFrame, mpoColName);
  //     this.mpoBestScoreViewer.dataFrame = this.dataFrame;
  //     this.mpoBestScoreViewer.onTableAttached();
  //     this.mpoContextPanel?.addPane(
  //       'Best scores',
  //       () => this.mpoBestScoreViewer!.root,
  //       true,
  //     );
  //   }

  //   if (!this.mpoWorstScoreViewer) {
  //     this.mpoWorstScoreViewer = new MpoScoreViewer(this.dataFrame, mpoColName, 'worst');
  //     this.mpoWorstScoreViewer.dataFrame = this.dataFrame;
  //     this.mpoWorstScoreViewer.onTableAttached();
  //     this.mpoContextPanel?.addPane(
  //       'Worst scores',
  //       () => this.mpoWorstScoreViewer!.root,
  //       true,
  //     );
  //   }
  // }

  // async updateMpoScoresHistogram(): Promise<void> {
  //   if (!this.mpoContextPanel) return;

  //   try {
  //     const columnNames = await this.calculateMpoScores();
  //     if (!columnNames.length)
  //       return;

  //     const titleSpan = this.mpoContextPanel.root.querySelector('span');
  //     if (titleSpan)
  //       titleSpan.innerText = 'MPO Context';

  //     if (!this.mpoHistogram) {
  //       this.mpoHistogram = DG.Viewer.histogram(this.dataFrame);
  //       ui.empty(this.mpoHistogramHost);
  //       this.mpoHistogramHost.appendChild(this.mpoHistogram.root);

  //       this.mpoContextPanel.addPane(
  //         'Score distribution',
  //         () => this.mpoHistogramHost,
  //         true,
  //       );
  //     }

  //     this.mpoHistogram.setOptions({
  //       valueColumnName: columnNames[0],
  //     });

  //     this.createScoresPanes(columnNames[0]);

  //     if (this.mpoBestScoreViewer)
  //       this.mpoBestScoreViewer.render();
  //     if (this.mpoWorstScoreViewer)
  //       this.mpoWorstScoreViewer.render();

  //     grok.shell.o = this.mpoContextPanel.root;
  //   } catch (e) {
  //     grok.shell.error(`Failed to update MPO histogram: ${e instanceof Error ? e.message : e}`);
  //   }
  // }

  private async listenForProfileChanges(): Promise<void> {
    grok.events.onCustomEvent(MPO_SCORE_CHANGED_EVENT).subscribe(async () => {
      // this.updateMpoScoresHistogram();
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
      const columnNames = await this.calculateMpoScores();
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
