export const DEFAULT_PTO: boolean = true;
export const DEFAULT_SEQUENCE_LENGTH: number = 23;
export const MAX_SEQUENCE_LENGTH: number = 35;
export const USER_STORAGE_KEY: string = 'SequenceTranslator';
export const EXAMPLE_MIN_WIDTH: string = '400px';

export const enum JSON_FIELD {
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

export const SS = 'SS' as const;
export const AS = 'AS' as const;
export type StrandType = typeof SS | typeof AS;

export const STRANDS: StrandType[] = [SS, AS];
export const STRAND_NAME: Record<StrandType, string> = {
  [SS]: 'Sense strand',
  [AS]: 'Anti sense',
};

export const THREE_PRIME = 'THREE_PRIME' as const;
export const FIVE_PRIME = 'FIVE_PRIME' as const;
export const TERMINAL_KEYS = [THREE_PRIME, FIVE_PRIME];
export const TERMINAL = {
  [THREE_PRIME]: 3,
  [FIVE_PRIME]: 5,
}
