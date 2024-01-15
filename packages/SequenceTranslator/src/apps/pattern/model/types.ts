/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PATTERN_KEY, TERMINI, STRANDS} from './const';

// todo: remove after full refactoring
export interface LegacyPatternConfig {
  [PATTERN_KEY.SS_BASES]: string[];
  [PATTERN_KEY.AS_BASES]: string[];
  [PATTERN_KEY.SS_PTO]:  boolean[];
  [PATTERN_KEY.AS_PTO]: boolean[];
  [PATTERN_KEY.SS_3]: string;
  [PATTERN_KEY.SS_5]: string;
  [PATTERN_KEY.AS_3]: string;
  [PATTERN_KEY.AS_5]: string;
  [PATTERN_KEY.COMMENT]: string;
}

export type StrandType = typeof STRANDS[number];
export type TerminalType = typeof TERMINI[number];

export type PatternConfiguration = {
  patternName: string,
  isAntisenseStrandIncluded: boolean,
  nucleotideSequences: Record<StrandType, string[]>,
  phosphorothioateLinkageFlags: Record<StrandType, boolean[]>,
  strandTerminusModifications: Record<StrandType, Record<TerminalType, string>>,
  patternComment: string,
  nucleotidesWithNumericLabels: string[],
}
