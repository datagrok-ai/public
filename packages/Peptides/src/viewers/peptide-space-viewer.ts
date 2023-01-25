import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {getSequenceMolecularWeight} from '../utils/molecular-measure';
import {AlignedSequenceEncoder} from '@datagrok-libraries/bio/src/sequence-encoder';
import {createDimensinalityReducingWorker, IReduceDimensionalityResult,
} from '@datagrok-libraries/ml/src/workers/dimensionality-reducing-worker-creator';
import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import * as C from '../utils/constants';
import {PeptidesModel} from '../model';

export class PeptideSpaceViewer extends DG.JsViewer {
  method: string;
  measure: string;
  cyclesCount: number;
  customProperties = new Set(['method', 'measure', 'cyclesCount']);
  isEmbeddingCreating: boolean = false;
  model!: PeptidesModel;
  //FIXME: even if the property stays the same, for some reason it's still triggering prop change
  prevProps!: {'method': string, 'measure': string, 'cyclesCount': number};

  constructor() {
    super();

    const methodChoices = ['UMAP', 't-SNE', 'SPE', 'pSPE', 'OriginalSPE'];
    this.method = this.addProperty('method', DG.TYPE.STRING, 'UMAP', {choices: methodChoices});
    const measureChoices = ['Levenshtein', 'Jaro-Winkler'];
    this.measure = this.addProperty('measure', DG.TYPE.STRING, 'Levenshtein', {choices: measureChoices});
    this.cyclesCount = this.addProperty('cyclesCount', DG.TYPE.INT, 100);
    this.updatePrevProperties();
  }

  updatePrevProperties(): void {
    this.prevProps = {'method': this.method, 'measure': this.measure, 'cyclesCount': this.cyclesCount};
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();

    this.model = PeptidesModel.getInstance(this.dataFrame);

    await this.render(!this.dataFrame.temp[C.EMBEDDING_STATUS]);
  }

  async onPropertyChanged(property: DG.Property | null): Promise<void> {
    super.onPropertyChanged(property);
    if (this.prevProps[property?.name as 'method' | 'measure' | 'cyclesCount' ?? ''] == property?.get(this))
      return;

    if (this.model)
      await this.render(this.customProperties.has(property?.name ?? '') || !this.dataFrame.temp[C.EMBEDDING_STATUS]);

    this.updatePrevProperties();
  }

  async render(computeData=false): Promise<void> {
    if (computeData && !this.isEmbeddingCreating && !this.model.isChangingEdfBitset) {
      this.isEmbeddingCreating = true;
      $(this.root).empty();
      const viewerHost = ui.waitBox(async () => {
        // const aligendSeqCol = this.dataFrame.columns.bySemType(C.SEM_TYPES.MACROMOLECULE)!;
        const alignedSeqCol = this.dataFrame.getCol(this.model.settings.sequenceColumnName!);
        const edf = await computeWeights(this.dataFrame, this.method, this.measure, this.cyclesCount, alignedSeqCol);
        this.dataFrame.temp[C.EMBEDDING_STATUS] = true;
        this.model.edf = edf;

        if (edf === null)
          return ui.label('Could not compute embeddings');

        const edfSelection = edf.selection;
        edfSelection.copyFrom(this.dataFrame.selection);
        edfSelection.onChanged.subscribe(() => {
          if (!this.model.isChangingEdfBitset)
            this.model.fireBitsetChanged(true);
        });

        const colorCol = this.dataFrame.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
        edf.columns.add(colorCol);

        const viewerOptions = {
          x: '~X', y: '~Y', color: colorCol.name ?? '~MW', size: '~MW', title: 'Peptide Space',
          showYSelector: false, showXSelector: false, showColorSelector: false, showSizeSelector: false,
          zoomAndFilter: 'no action', axesFollowFilter: false,
        };
        const scatterPlot = edf.plot.scatter(viewerOptions);
        const viewerRoot = scatterPlot.root;

        viewerRoot.addEventListener('mousemove', (ev) => {
          const idx = scatterPlot.hitTest(ev.offsetX, ev.offsetY);
          if (idx != -1) {
            const table = ui.tableFromMap({
              'Activity': colorCol.get(idx),
              'Sequence': alignedSeqCol.get(idx),
              'Row ID': idx,
            });
            ui.tooltip.show(table, ev.clientX, ev.clientY);
          }
        });
        viewerRoot.style.width = 'auto';
        this.isEmbeddingCreating = false;
        viewerHost.style.paddingLeft = 'unset';
        return viewerRoot;
      }) as HTMLDivElement;
      viewerHost.style.paddingLeft = '45%';
      this.root.appendChild(viewerHost);
    }
  }
}

//Do not accept table, only column
export async function computeWeights(
  table: DG.DataFrame, method: string, measure: string, cyclesCount: number, col: DG.Column,
): Promise<DG.DataFrame | null> {
  const pi = DG.TaskBarProgressIndicator.create('Creating embedding...');
  let edf: DG.DataFrame | null = null;
  try {
    const axesNames = ['~X', '~Y', '~MW'];
    // col ??= table.columns.bySemType(C.SEM_TYPES.MACROMOLECULE)!;
    // col ??= table.getCol(this.model.settings.sequenceColumnName!);
    const columnData = col.toList().map((v) => AlignedSequenceEncoder.clean(v));

    const reduceDimRes: IReduceDimensionalityResult = await createDimensinalityReducingWorker(
      {data: columnData, metric: measure as StringMetrics}, method, {cycles: cyclesCount});
    const embcols = reduceDimRes.embedding;

    const columns = Array.from(embcols, (v: Float32Array, k) => DG.Column.fromFloat32Array(axesNames[k], v));

    function _getMW(sequences: string[]): Float32Array {
      const mw: Float32Array = new Float32Array(sequences.length);

      mw.map((_, index) => getSequenceMolecularWeight(sequences[index] ?? ''));

      return mw;
    }

    columns.push(DG.Column.fromFloat32Array('~MW', _getMW(columnData)));

    edf = DG.DataFrame.fromColumns(columns);
  } catch (error) {
    grok.shell.error('Could not compute embeddings. See console for details.');
    console.error(error);
  } finally {
    pi.close();
    return edf;
  }
}
