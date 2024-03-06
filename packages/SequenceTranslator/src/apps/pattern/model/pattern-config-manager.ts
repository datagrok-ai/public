import {EventBus} from './event-bus';
import {PatternConfiguration} from './types';

/** A wrapper over EventBus encapsulating pattern config access  */
export class PatternConfigurationManager {
  static getConfig(eventBus: EventBus): PatternConfiguration {
    const patternName = eventBus.getPatternName();
    const isAntisenseStrandIncluded = eventBus.isAntiSenseStrandVisible();
    const nucleotideSequences = eventBus.getNucleotideSequences();
    const phosphorothioateLinkageFlags = eventBus.getPhosphorothioateLinkageFlags();
    const strandTerminusModifications = eventBus.getTerminalModifications();
    const patternComment = eventBus.getComment();
    const nucleotidesWithNumericLabels =
      eventBus.getModificationsWithNumericLabels();

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
