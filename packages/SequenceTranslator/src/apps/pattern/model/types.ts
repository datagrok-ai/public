/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  TERMINI, STRANDS,
  GRAPH_SETTINGS_KEYS as G, LEGEND_SETTINGS_KEYS as L, PATTERN_RECORD_KEYS as R
} from './const';

export type StrandType = typeof STRANDS[number];
export type TerminalType = typeof TERMINI[number];

export type NucleotideSequences = Record<StrandType, string[]>;
export type PhosphorothioateLinkageFlags = Record<StrandType, boolean[]>;
export type StrandTerminusModifications = Record<StrandType, Record<TerminalType, string>>;

export type PatternGraphSettings = {
  [G.IS_ANTISENSE_STRAND_INCLUDED]: boolean,
  [G.MODIFICATION_LABELS_VISIBLE]: boolean,
  [G.NUCLEOTIDE_SEQUENCES]: NucleotideSequences,
  [G.PHOSPHOROTHIOATE_LINKAGE_FLAGS]: PhosphorothioateLinkageFlags,
  [G.STRAND_TERMINUS_MODIFICATIONS]: StrandTerminusModifications,
}

export type PatternLegendSettings = {
  [L.PATTERN_NAME]: string,
  [L.PATTERN_COMMENT]: string,
  [L.NUCLEOTIDES_WITH_NUMERIC_LABELS]: string[],
}

export type PatternConfiguration = PatternGraphSettings & PatternLegendSettings;

export type PatternConfigRecord = {
  [R.PATTERN_CONFIG]: PatternConfiguration,
  [R.AUTHOR_ID]: string,
}

export class PatternNameExistsError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PatternNameExistsError';
  }
}

export class PatternExistsError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PatternExistsError';
  }
}

export type RawPatternRecords = {[patternName: string]: string};
