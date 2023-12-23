/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  DEFAULT_PHOSPHOROTHIOATE, DEFAULT_SEQUENCE_LENGTH, MAX_SEQUENCE_LENGTH, USER_STORAGE_KEY, SENSE_STRAND, ANTISENSE_STRAND, STRAND_LABEL, STRANDS, TERMINAL_KEYS, TERMINAL, THREE_PRIME_END, FIVE_PRIME_END, PATTERN_FIELD as PATTERN_KEY, StrandType, TerminalType
} from '../../../model/pattern-app/const';

export interface PatternConfiguration {
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

export type BooleanInput = DG.InputBase<boolean | null>;
export type StringInput = DG.InputBase<string | null>;
export type NumberInput = DG.InputBase<number | null>;
