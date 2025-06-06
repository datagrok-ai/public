import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {MmDistanceFunctionsNames} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';
import {TAGS as bioTAGS} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {SearchBaseViewer} from '@datagrok-libraries/ml/src/viewers/search-base-viewer';

const MAX_ROWS_FOR_DISTANCE_MATRIX = 10000; // Reduced from 22000 to be safer
const MAX_ROWS_FOR_FULL_CALCULATION = 15000; // Threshold for switching to sampling

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

    // Initialize with macromolecule distance functions instead of chemical similarity metrics
    this.distanceMetric = this.string('distanceMetric', this.distanceFunctionChoices[0], {
      choices: this.distanceFunctionChoices
    });

    this.fingerprint = this.string('fingerprint', this.fingerprintChoices[0], {
      choices: this.fingerprintChoices
    });

    // Gap penalty parameters for Needleman-Wunsch
    this.gapOpen = this.float('gapOpen', -10.0);
    this.gapExtend = this.float('gapExtend', -0.5);
  }

  async onTableAttached(): Promise<void> {
    super.onTableAttached();

    if (this.dataFrame) {
      const rowCount = this.dataFrame.rowCount;
      this.preComputeDistanceMatrix = rowCount <= MAX_ROWS_FOR_DISTANCE_MATRIX;
      this.requiresSampling = rowCount > MAX_ROWS_FOR_FULL_CALCULATION;

      if (this.requiresSampling)
        grok.shell.info(`Large dataset detected (${rowCount} rows). Will use sampling for diversity calculation.`);
    }
  }

  // Helper method to get distance function options for current selection
  getDistanceFunctionOptions(): {[key: string]: any} {
    const options: {[key: string]: any} = {};

    if (this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH) {
      options.gapOpen = this.gapOpen;
      options.gapExtend = this.gapExtend;
    }

    if (this.distanceMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE ||
        this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH)
      options.fingerprintType = this.fingerprint;


    return options;
  }

  // Helper method to check if gap penalties are needed
  needsGapPenalties(): boolean {
    return this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  }

  // Helper method to check if fingerprint is needed
  needsFingerprint(): boolean {
    return this.distanceMetric === MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE ||
           this.distanceMetric === MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH;
  }
}
