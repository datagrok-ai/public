import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SearchBaseViewer} from '@datagrok-libraries/ml/src/viewers/search-base-viewer';

const MAX_ROWS_FOR_DISTANCE_MATRIX = 10000;

export class SequenceSearchBaseViewer extends SearchBaseViewer {
  distanceMetric: string;
  fingerprint: string;
  gapOpen: number;
  gapExtend: number;

  metricsProperties = ['distanceMetric', 'fingerprint', 'gapOpen', 'gapExtend'];
  fingerprintChoices = ['Morgan', 'RDKit', 'Pattern', 'AtomPair', 'MACCS', 'TopologicalTorsion'];
  distanceFunctionChoices = [
    MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH,
    MmDistanceFunctionsNames.HAMMING,
    MmDistanceFunctionsNames.LEVENSHTEIN,
    MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE
  ];

  tags = [DG.TAGS.UNITS, bioTAGS.aligned, bioTAGS.separator, bioTAGS.alphabet, 'cell.renderer'];
  preComputeDistanceMatrix: boolean = false;
  requiresSampling: boolean = false;

  constructor(name: string, semType: string) {
    super(name, semType);

    this.distanceMetric = this.string('distanceMetric', MmDistanceFunctionsNames.HAMMING, {
      choices: this.distanceFunctionChoices
    });

    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {
      choices: this.fingerprintChoices
    });

    this.gapOpen = this.float('gapOpen', 1);
    this.gapExtend = this.float('gapExtend', 0.6);
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();

    if (this.dataFrame) {
      const rowCount = this.dataFrame.rowCount;
      this.preComputeDistanceMatrix = rowCount <= MAX_ROWS_FOR_DISTANCE_MATRIX;
      this.requiresSampling = rowCount > MAX_ROWS_FOR_DISTANCE_MATRIX;
    }
  }

  needsGapPenalties(): boolean {
    return this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  }

  needsFingerprint(): boolean {
    return this.distanceMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE ||
           this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  }
}
