import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DesirabilityProfile, WeightedAggregation} from '@datagrok-libraries/statistics/src/mpo/mpo';

import {MpoScoreViewer} from './mpo-scores-viewer';
import {computeMpo} from './utils';
import {Subscription} from 'rxjs';


export class MpoContextPanel {
  root: HTMLElement;
  df: DG.DataFrame;

  private panel: DG.Accordion = ui.accordion();
  private histogramHost: HTMLDivElement = ui.div();
  private histogram?: DG.Viewer;
  private bestScoreViewer?: MpoScoreViewer;
  private worstScoreViewer?: MpoScoreViewer;
  private currentObjectChangingSub: Subscription | null = null;

  // Uncomment to disable interactivity
  // showCurrentRow: false,
  // showMouseOverRow: false,
  // showMouseOverRowGroup: false,
  private static HISTOGRAM_OPTIONS = {
    showColumnSelector: false,
    showSplitSelector: false,
    showRangeSlider: false,
    showXAxis: true,
    allowColumnSelection: false,
    showBinSelector: false,
  };

  constructor(df: DG.DataFrame) {
    grok.shell.windows.showHelp = false;
    this.df = df;
    this.root = this.panel.root;

    const icon = ui.element('i');
    icon.className = 'grok-icon svg-icon svg-histogram';
    this.panel.addTitle(ui.span([icon, ui.label('MPO not calculated yet')]));

    this.histogramHost.classList.add('chem-mpo-histogram-compact');
  }

  show(): void {
    grok.shell.o = this.root;
    this.attachCurrentObjectChanging();
  }

  private updateTitle(title: string) {
    const titleSpan = this.panel.root.querySelector('span');
    if (titleSpan)
      titleSpan.innerText = title;
  }

  private createHistogram(columnName: string): DG.Viewer {
    const hist = DG.Viewer.histogram(this.df, {
      ...MpoContextPanel.HISTOGRAM_OPTIONS,
      valueColumnName: columnName,
    });

    // Uncomment to disable interactivity
    // ['mousedown', 'click', 'dblclick'].forEach((ev) =>
    //   hist.root.addEventListener(ev, (e) => e.stopPropagation(), true),
    // );
    return hist;
  }

  private createScoreViewer(columnName: string, type: 'best' | 'worst' = 'best'): MpoScoreViewer {
    const viewer = new MpoScoreViewer(this.df, columnName, type === 'worst' ? 'worst' : undefined);
    viewer.dataFrame = this.df;
    viewer.onTableAttached();
    return viewer;
  }

  async render(
    profile: DesirabilityProfile,
    columnMapping: Record<string, string | null>,
    aggregation?: WeightedAggregation,
    formula?: string,
  ) {
    if (!profile)
      return;

    if (grok.shell.o !== this.root)
      this.show();

    const columnNames = await computeMpo(this.df, profile, columnMapping, aggregation, formula);
    if (!columnNames.length)
      return;

    this.updateTitle('MPO Profile');

    if (!this.histogram) {
      this.histogram = this.createHistogram(columnNames[0]);
      this.histogramHost.appendChild(this.histogram.root);
      this.panel.addPane('Scores', () => this.histogramHost, true);
    }
    this.histogram.apply({valueColumnName: columnNames[0]});

    if (!this.bestScoreViewer) {
      this.bestScoreViewer = this.createScoreViewer(columnNames[0], 'best');
      this.panel.addPane('Best scores', () => this.bestScoreViewer!.root, true);
    }

    if (!this.worstScoreViewer) {
      this.worstScoreViewer = this.createScoreViewer(columnNames[0], 'worst');
      this.panel.addPane('Worst scores', () => this.worstScoreViewer!.root, true);
    }

    this.bestScoreViewer.render();
    this.worstScoreViewer.render();
  }

  private attachCurrentObjectChanging(): void {
    if (this.currentObjectChangingSub)
      return;

    this.currentObjectChangingSub = grok.events.onEvent('d4-current-object-changing').subscribe((e) => {
      if (e.newObject !== this.root)
        e.preventDefault();
    });
  }

  close(): void {
    this.currentObjectChangingSub?.unsubscribe();
    this.currentObjectChangingSub = null;
    this.resetViewers();

    if (grok.shell.o === this.root)
      grok.shell.o = null;
  }

  updateDataFrame(df: DG.DataFrame): void {
    this.df = df;
    this.resetViewers();
  }

  private resetViewers(): void {
    this.histogram = undefined;
    ui.empty(this.histogramHost);
    this.bestScoreViewer = undefined;
    this.worstScoreViewer = undefined;

    while (this.panel.panes.length > 0)
      this.panel.removePane(this.panel.panes[0]);
  }
}
