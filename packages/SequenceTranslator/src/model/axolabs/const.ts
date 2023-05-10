export const defaultPto: boolean = true;
export const defaultSequenceLength: number = 23;
export const maximalValidSequenceLength: number = 35;
export const userStorageKey: string = 'SequenceTranslator';
export const exampleMinWidth: string = '400px';

export const IDX = {
  SS: 0,
  AS: 1,
  THREE_PRIME: 0,
  FIVE_PRIME: 1,
};

export const strands = ['SS', 'AS'] as const;
export const strandLongNames = ['Sense Strand', 'Antisense Strand'] as const;
export const terminals = [3, 5] as const;
