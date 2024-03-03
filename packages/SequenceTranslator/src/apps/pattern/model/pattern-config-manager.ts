import {EventBus} from './event-bus';
import {PatternConfiguration} from './types';

/** A wrapper over EventBus encapsulating pattern config access  */
export class PatternConfigurationManager {
  constructor(private eventBus: EventBus) { }

  getConfig(): PatternConfiguration {
    const patternName = this.eventBus.getPatternName();
    const isAntisenseStrandIncluded = this.eventBus.isAntiSenseStrandVisible();
    const nucleotideSequences = this.eventBus.getNucleotideSequences();
    const phosphorothioateLinkageFlags = this.eventBus.getPhosphorothioateLinkageFlags();
    const strandTerminusModifications = this.eventBus.getTerminalModifications();
    const patternComment = this.eventBus.getComment();
    const nucleotidesWithNumericLabels =
      this.eventBus.getModificationsWithNumericLabels();

    return {
      patternName,
      isAntisenseStrandIncluded,
      nucleotideSequences,
      phosphorothioateLinkageFlags,
      strandTerminusModifications,
      patternComment,
      nucleotidesWithNumericLabels
    };
  }
}
