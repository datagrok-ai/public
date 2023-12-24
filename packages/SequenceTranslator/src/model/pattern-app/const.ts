export const DEFAULT_PHOSPHOROTHIOATE: boolean = true;
export const DEFAULT_SEQUENCE_LENGTH: number = 23;
export const MAX_SEQUENCE_LENGTH: number = 35;
export const USER_STORAGE_KEY: string = 'SequenceTranslator';
export const EXAMPLE_MIN_WIDTH: string = '400px';

// todo: remove as legacy
export const enum PATTERN_KEY {
  SS_BASES = 'ssBases',
  AS_BASES = 'asBases',
  SS_PTO = 'ssPtoLinkages',
  AS_PTO = 'asPtoLinkages',
  SS_3 = 'ssThreeModification',
  SS_5 = 'ssFiveModification',
  AS_3 = 'asThreeModification',
  AS_5 = 'asFiveModification',
  COMMENT = 'comment',
};

export const SENSE_STRAND = 'SENSE_STRAND' as const;
export const ANTISENSE_STRAND = 'ANTISENSE_STRAND' as const;
export type StrandType = typeof SENSE_STRAND | typeof ANTISENSE_STRAND;

export const STRANDS: StrandType[] = [SENSE_STRAND, ANTISENSE_STRAND];
export const STRAND_LABEL: Record<StrandType, string> = {
  [SENSE_STRAND]: 'Sense strand',
  [ANTISENSE_STRAND]: 'Anti sense',
};

export const THREE_PRIME_END = 'THREE_PRIME_END' as const;
export const FIVE_PRIME_END = 'FIVE_PRIME_END' as const;
export type TerminalType = typeof THREE_PRIME_END | typeof FIVE_PRIME_END;

export const TERMINAL_KEYS: TerminalType[] = [THREE_PRIME_END, FIVE_PRIME_END];
export const TERMINAL: Record<TerminalType, number> = {
  [THREE_PRIME_END]: 3,
  [FIVE_PRIME_END]: 5,
};
