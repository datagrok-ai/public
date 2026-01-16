import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DesirabilityProfile, PropertyDesirability, WeightedAggregation} from '@datagrok-libraries/statistics/src/mpo/mpo';
import {MpoScoreViewer} from '../apps/mpo-scores';
import {value} from '../widgets/scaffold-tree';

export class MpoContextPanel {
  root: HTMLElement;
  df: DG.DataFrame;
  private panel: DG.Accordion;
  private histogramHost: HTMLDivElement = ui.div();
  private histogram?: DG.Viewer;
  private bestScoreViewer?: MpoScoreViewer;
  private worstScoreViewer?: MpoScoreViewer;

  constructor(df: DG.DataFrame) {
    grok.shell.windows.showHelp = false;
    this.df = df;
    this.panel = ui.accordion();
    const icon = ui.element('i');
    icon.className = 'grok-icon svg-icon svg-histogram';
    this.panel.addTitle(ui.span([icon, ui.label('MPO not calculated yet')]));
    this.root = this.panel.root;
    this.histogramHost.classList.add('chem-mpo-histogram-compact');
    grok.shell.o = this.root;
  }

  async render(
    profile: DesirabilityProfile, columnMapping: Record<string, string | null>, aggregation?: WeightedAggregation,
  ) {
    if (!profile) return;

    if (grok.shell.o !== this.root)
      grok.shell.o = this.root;

    // Prepare mapped properties
    const mappedProperties: Record<string, PropertyDesirability> = {};
    for (const [propName, prop] of Object.entries(profile.properties)) {
      const columnName = columnMapping[propName] ?? propName;
      mappedProperties[columnName] = prop;
    }

    // Call MPO calculation function
    const [func] = await DG.Func.find({name: 'mpoTransformFunction'});
    const funcCall = await func.prepare({
      df: this.df,
      profileName: profile.name ?? 'MPO',
      currentProperties: mappedProperties,
      aggregation: aggregation ?? null,
    }).call(undefined, undefined, {processed: false});

    const columnNames: string[] = funcCall.getOutputParamValue();
    if (!columnNames.length) return;

    // Update title
    const titleSpan = this.panel.root.querySelector('span');
    if (titleSpan)
      titleSpan.innerText = 'MPO Profile';

    // Histogram
    if (!this.histogram) {
      this.histogram = DG.Viewer.histogram(this.df, {
        showColumnSelector: false,
        showSplitSelector: false,
        showRangeSlider: false,
        showXAxis: true,
        allowColumnSelection: false,
        showBinSelector: false,
        showCurrentRow: false,
        showMouseOverRow: false,
        showMouseOverRowGroup: false,
        valueColumnName: columnNames[0],
      });

      this.histogram.root.addEventListener('mousedown', (e) => e.stopPropagation(), true);
      this.histogram.root.addEventListener('click', (e) => e.stopPropagation(), true);
      this.histogram.root.addEventListener('dblclick', (e) => e.stopPropagation(), true);

      this.histogramHost.appendChild(this.histogram.root);
      this.panel.addPane('Scores', () => this.histogramHost, true);
    }

    this.histogram.apply({valueColumnName: columnNames[0]});

    // Best/Worst Scores
    if (!this.bestScoreViewer) {
      this.bestScoreViewer = new MpoScoreViewer(this.df, columnNames[0]);
      this.bestScoreViewer.dataFrame = this.df;
      this.bestScoreViewer.onTableAttached();
      this.panel.addPane('Best scores', () => this.bestScoreViewer!.root, true);
    }
    if (!this.worstScoreViewer) {
      this.worstScoreViewer = new MpoScoreViewer(this.df, columnNames[0], 'worst');
      this.worstScoreViewer.dataFrame = this.df;
      this.worstScoreViewer.onTableAttached();
      this.panel.addPane('Worst scores', () => this.worstScoreViewer!.root, true);
    }

    this.bestScoreViewer?.render();
    this.worstScoreViewer?.render();

    grok.shell.o = this.root;
  }
}
