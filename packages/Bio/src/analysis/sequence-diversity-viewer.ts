import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {similarityMetric, getDiverseSubset} from '@datagrok-libraries/utils/src/similarity-metrics';
import $ from 'cash-dom';
import {ArrayUtils} from '@datagrok-libraries/utils/src/array-utils';
import { SequenceSearchBaseViewer } from './sequence-search-base-viewer';

export class SequenceDiversityViewer extends SequenceSearchBaseViewer {
  renderMolIds: number[];
  columnNames = [];

  constructor() {
    super('diversity');
    this.renderMolIds = [];
    this.updateMetricsLink(this.metricsDiv, this, {fontSize: '10px', fontWeight: 'normal', paddingBottom: '15px'});
  }


  async render(computeData = true): Promise<void> {
    if (!this.beforeRender())
      return;
    if (this.dataFrame) {
      if (computeData) {
        this.renderMolIds =
        await grok.functions.call('Chem:chemDiversitySearch', {
            col: this.moleculeColumn,
            metricName: this.distanceMetric,
            limit: this.limit,
            fingerprint: this.fingerprint
        });
      }
    }
  }
}
