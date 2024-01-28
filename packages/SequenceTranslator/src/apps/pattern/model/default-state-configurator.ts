import {PatternAppDataManager} from './external-data-manager';
import {STRANDS, TERMINI, DEFAULT_PHOSPHOROTHIOATE} from './const';
import {TerminalType, NucleotideSequences, StrandType, PhosphorothioateLinkageFlags, StrandTerminusModifications} from './types';
import {DEFAULT_PATTERN_CONFIG as DEFAULT} from './const';

export class DefaultStateConfigurator {
  constructor(
    private dataManager: PatternAppDataManager
  ) { }

  getPatternName(): string { return DEFAULT.PATTERN_NAME; }

  getAntiSenseStrandVisibilityFlag(): boolean { return DEFAULT.IS_ANTISENSE_STRAND_VISIBLE; }

  getNucleobases(): NucleotideSequences {
    const nucleobases = {} as NucleotideSequences;
    const defaultNucleobase = this.fetchDefaultNucleobase();
    STRANDS.forEach((strand) => {
      nucleobases[strand] = new Array(DEFAULT.SEQUENCE_LENGTH).fill(defaultNucleobase);
    });

    return nucleobases;
  }

  private fetchDefaultNucleobase(): string {
    return this.dataManager.fetchAvailableNucleotideBases()[0];
  }
  
  getPhosphorothioateLinkageFlags(): PhosphorothioateLinkageFlags {
    const phosphorothioateLinkageFlags = {} as PhosphorothioateLinkageFlags;
    STRANDS.forEach((strand) => {
      phosphorothioateLinkageFlags[strand] = new Array(DEFAULT.SEQUENCE_LENGTH).fill(DEFAULT.PHOSPHOROTHIOATE);
    });

    return phosphorothioateLinkageFlags;
  }
  
  getTerminusModifications(): StrandTerminusModifications {
    const terminusModifications = {} as StrandTerminusModifications;
    STRANDS.forEach((strand) => {
      terminusModifications[strand] = {} as Record<TerminalType, string>;
      TERMINI.forEach((terminus) => {
        terminusModifications[strand][terminus] = DEFAULT.TERMINUS_MODIFICATION;
      });
    });

    return terminusModifications;
  }

  getComment(): string { return DEFAULT.COMMENT; }
  
  getModificationsWithNumericLabels(): string[] { return DEFAULT.MODIFICATIONS_WITH_NUMERIC_LABELS; }
}
