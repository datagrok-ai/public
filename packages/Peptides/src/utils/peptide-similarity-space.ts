/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
// eslint-disable-next-line no-unused-vars
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getSequenceMolecularWeight} from './molecular-measure';
import {AlignedSequenceEncoder} from '@datagrok-libraries/utils/src/sequence-encoder';
import {DimensionalityReducer} from '@datagrok-libraries/utils/src/reduce-dimensionality';
import {Measurer} from '@datagrok-libraries/utils/src/string-measure';
import {Coordinates} from '@datagrok-libraries/utils/src/type-declarations';

function createDimensinalityReducingWorker(
  columnData: any[],
  method: string,
  measure: string,
  cyclesCount: number,
): Promise<unknown> {
  return new Promise(function(resolve) {
    const worker = new Worker(new URL('../workers/dimensionality-reducer.ts', import.meta.url));
    worker.postMessage({
      columnData: columnData,
      method: method,
      measure: measure,
      cyclesCount: cyclesCount,
    });
    worker.onmessage = ({data: {embedding}}) => {
      resolve(embedding);
    };
  });
}

/**
 * Finds a column with an activity.
 *
 * @param {DG.DataFrame} table The data frame to search for.
 * @return {(string | null)} Column name or null if not found.
 */
function inferActivityColumnsName(table: DG.DataFrame): string | null {
  const re = /activity|ic50/i;
  for (const name of table.columns.names()) {
    if (name.match(re)) {
      console.log(`${name} found.`);
      return name;
    }
  }
  return null;
}

/**
 * Creates scatter plot with sequences embeded.
 *
 * @export
 * @param {DG.DataFrame} table The table containing samples.
 * @param {DG.Column} alignedSequencesColumn Samples column.
 * @param {string} method Embedding method to apply.
 * @param {string} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @param {(DG.TableView | null)} view View to add scatter plot to
 * @param {(string | null)} [activityColumnName] Activity containing column to assign it to points radius.
 * @param {boolean} [zoom=false] Whether to fit view.
 * @return {Promise<DG.ScatterPlotViewer>} A viewer.
 */
export async function createPeptideSimilaritySpaceViewer(
  table: DG.DataFrame,
  alignedSequencesColumn: DG.Column,
  method: string,
  measure: string,
  cyclesCount: number,
  view: DG.TableView | null,
  activityColumnName?: string | null,
  zoom: boolean = false,
): Promise<DG.ScatterPlotViewer> {
  const pi = DG.TaskBarProgressIndicator.create('Creating embedding.');

  activityColumnName = activityColumnName ?? inferActivityColumnsName(table);

  const axesNames = ['~X', '~Y', '~MW'];
  const columnData = alignedSequencesColumn.toList().map((v, _) => AlignedSequenceEncoder.clean(v));

  const embcols = await createDimensinalityReducingWorker(columnData, method, measure, cyclesCount);

  const columns = Array.from(
    embcols as Coordinates,
    (v: Float32Array, k) => (DG.Column.fromFloat32Array(axesNames[k], v)),
  );

  function _getMW(sequences = columnData) {
    const mw: Float32Array = new Float32Array(sequences.length).fill(0);
    let currentSequence;

    for (let i = 0; i < sequences.length; ++i) {
      currentSequence = sequences[i];
      mw[i] = currentSequence == null ? 0 : getSequenceMolecularWeight(currentSequence);
    }
    return mw;
  }

  columns.push(DG.Column.fromFloat32Array('~MW', _getMW()));

  const edf = DG.DataFrame.fromColumns(columns);

  // Add new axes.
  for (const axis of axesNames) {
    const col = table.col(axis);

    if (col == null) {
      table.columns.insert(edf.getCol(axis));
    } else {
      table.columns.replace(col, edf.getCol(axis));
    }
  }

  // const viewer = DG.Viewer.scatterPlot(table, {x: '~X', y: '~Y', color: activityColumnName ?? '~MW', size: '~MW'});
  const viewerOptions = {x: '~X', y: '~Y', color: activityColumnName ?? '~MW', size: '~MW'};
  const viewer = view !== null ?
    view.addViewer(DG.VIEWER.SCATTER_PLOT, viewerOptions) : DG.Viewer.scatterPlot(table, viewerOptions);
  // Fit view if needed.
  /*if (zoom) {
    viewer.zoom(
      table.getCol('~X').min,
      table.getCol('~Y').min,
      table.getCol('~X').max,
      table.getCol('~Y').max,
    );
  }*/
  pi.close();
  return (viewer as DG.ScatterPlotViewer);
}

/**
 * Controls creation of the peptide similarity space viewer.
 *
 * @export
 * @class PeptideSimilaritySpaceWidget
 */
export class PeptideSimilaritySpaceWidget {
  protected method: string;
  protected metrics: string;
  protected cycles: number = 100;
  protected currentDf: DG.DataFrame;
  protected alignedSequencesColumn: DG.Column;
  protected availableMethods: string[];
  protected availableMetrics: string[];
  protected viewer: HTMLElement;
  view: DG.TableView;

  /**
   * Creates an instance of PeptideSimilaritySpaceWidget.
   * @param {DG.Column} alignedSequencesColumn The column to get amino acid sequences from.
   * @param {DG.TableView} view Current view
   * @memberof PeptideSimilaritySpaceWidget
   */
  constructor(alignedSequencesColumn: DG.Column, view: DG.TableView) {
    this.availableMethods = DimensionalityReducer.availableMethods;
    this.availableMetrics = Measurer.availableMeasures;
    this.method = this.availableMethods[0];
    this.metrics = this.availableMetrics[0];
    this.currentDf = alignedSequencesColumn.dataFrame;
    this.alignedSequencesColumn = alignedSequencesColumn;
    this.viewer = ui.div([]);
    this.view = view;
  }

  /**
   * Creates viewer itself.
   *
   * @return {Promise<DG.Viewer>} the viewer.
   * @memberof PeptideSimilaritySpaceWidget
   */
  public async drawViewer(): Promise<DG.Viewer> {
    const viewer = await createPeptideSimilaritySpaceViewer(
      this.currentDf,
      this.alignedSequencesColumn,
      this.method,
      this.metrics,
      this.cycles,
      null,
      null,
      true,
    );
    viewer.root.style.width = 'auto';
    return viewer;
  }

  /**
   * Updates the viewer on options changes.
   *
   * @protected
   * @memberof PeptideSimilaritySpaceWidget
   */
  protected async updateViewer() {
    this.viewer.lastChild?.remove();
    const viewer = await this.drawViewer();
    this.viewer.appendChild(viewer.root);
  }

  /**
   * Adds input controls to manage the viewer's parameters.
   *
   * @protected
   * @return {Promise<HTMLElement>} Bunch of control elements.
   * @memberof PeptideSimilaritySpaceWidget
   */
  protected async drawInputs(): Promise<HTMLElement> {
    const methodsList = ui.choiceInput('Embedding method', this.method, this.availableMethods,
      async (currentMethod: string) => {
        this.method = currentMethod;
        await this.updateViewer();
      },
    );
    methodsList.setTooltip('Embedding method to apply to the dataset.');

    const metricsList = ui.choiceInput('Distance metric', this.metrics, this.availableMetrics,
      async (currentMetrics: string) => {
        this.metrics = currentMetrics;
        await this.updateViewer();
      },
    );
    metricsList.setTooltip('Custom distance metric to pass to the embedding procedure.');

    const cyclesSlider = ui.intInput('Cycles count', this.cycles,
      async (currentCycles: number) => {
        this.cycles = currentCycles;
        await this.updateViewer();
      },
    );
    cyclesSlider.setTooltip('Number of cycles affects the embedding quality.');

    return ui.inputs([methodsList, metricsList, cyclesSlider]);
  }

  /**
   * Draws a viewer on property panel.
   *
   * @return {Promise<DG.Widget>} The corresponding widget.
   * @memberof PeptideSimilaritySpaceWidget
   */
  public async draw(): Promise<DG.Widget> {
    return new DG.Widget(ui.divV([(await this.drawViewer()).root, await this.drawInputs()]));
  }
}
