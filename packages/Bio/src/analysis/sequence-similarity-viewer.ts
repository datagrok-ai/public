import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import $ from 'cash-dom';
import { SequenceSearchBaseViewer } from './sequence-search-base-viewer';
import { getMonomericMols } from '../calculations/monomerLevelMols';

export class SequenceSimilarityViewer extends SequenceSearchBaseViewer {
  hotSearch: boolean;
  sketchedMolecule: string = '';
  curIdx: number = 0;
  molCol: DG.Column | null = null;
  idxs: DG.Column | null = null;
  scores: DG.Column | null = null;
  cutoff: number;
  gridSelect: boolean = false;
  targetMoleculeIdx: number = 0;

  constructor() {
    super('similarity');
    this.cutoff = this.float('cutoff', 0.01, {min: 0, max: 1});
    this.hotSearch = this.bool('hotSearch', true);
    this.updateMetricsLink(this.metricsDiv, this, {fontSize: '10px', fontWeight: 'normal', height: '10px'});
  }

  init(): void {
    this.hotSearch = true;
    this.initialized = true;
  }

  async render(computeData = true): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.moleculeColumn) {
      this.curIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
      if (computeData && !this.gridSelect) {
        this.targetMoleculeIdx = this.dataFrame!.currentRowIdx == -1 ? 0 : this.dataFrame!.currentRowIdx;
        const monomericMols = await getMonomericMols(this.moleculeColumn);
        const df = await grok.functions.call('Chem:callChemSimilaritySearch', {
            df: this.dataFrame,
            col: monomericMols,
            molecule: monomericMols.get(this.targetMoleculeIdx),
            metricName: this.distanceMetric,
            limit: this.limit,
            minScore: this.cutoff,
            fingerprint: this.fingerprint
        });
        this.molCol = df.getCol('smiles');
        this.idxs = df.getCol('indexes');
        this.scores = df.getCol('score');
        this.root.append(df.plot.grid().root);
      } else if (this.gridSelect)
        this.gridSelect = false;
    }
  }
}
