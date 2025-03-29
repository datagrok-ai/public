import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {CHEM_SIMILARITY_METRICS} from '@datagrok-libraries/ml/src/distance-metrics-methods';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SearchBaseViewer} from '@datagrok-libraries/ml/src/viewers/search-base-viewer';

const MAX_ROWS_FOR_DISTANCE_MATRIX = 22000;

export class SequenceSearchBaseViewer extends SearchBaseViewer {
  distanceMetric: string;
  fingerprint: string;
  metricsProperties = ['distanceMetric', 'fingerprint'];
  fingerprintChoices = ['Morgan', 'Pattern'];
  tags = [DG.TAGS.UNITS, bioTAGS.aligned, bioTAGS.separator, bioTAGS.alphabet];
  preComputeDistanceMatrix: boolean = false;

  constructor(name: string, semType: string) {
    super(name, semType);
    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {choices: this.fingerprintChoices});
    this.distanceMetric = this.string('distanceMetric', CHEM_SIMILARITY_METRICS[0], {choices: CHEM_SIMILARITY_METRICS});
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();

    if (this.dataFrame)
      this.preComputeDistanceMatrix = this.dataFrame.rowCount <= MAX_ROWS_FOR_DISTANCE_MATRIX;
  }
}
