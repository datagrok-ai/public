export enum CDDVaultSearchType {
    EXACT_MATCH = 'exact',
    SIMILARITY = 'similarity',
    SUBSTRUCTURE = 'substructure',
}

export enum ProtocolRunCond {
    ANY = '(any run)',
    SPECIFIC = 'specific run',
}

export const CDD_SEARCH_TYPES = [CDDVaultSearchType.SUBSTRUCTURE, CDDVaultSearchType.SIMILARITY,
  CDDVaultSearchType.EXACT_MATCH];

export const PROTOCOL_RUN_COND_CHOICES = [ProtocolRunCond.ANY, ProtocolRunCond.SPECIFIC];

export const ANY_RUN_CHOICE = '(any run)';

export const PROTOCOLS_TAB = 'Protocols';
export const COLLECTIONS_TAB = 'Collections';
export const SAVED_SEARCHES_TAB = 'Saved Searches';
export const MOLECULES_TAB = 'Molecules';
export const BATCHES_TAB = 'Batches';
export const SEARCH_TAB = 'Search';

export const ALL_TABS = [PROTOCOLS_TAB, COLLECTIONS_TAB, SAVED_SEARCHES_TAB, MOLECULES_TAB, BATCHES_TAB, SEARCH_TAB];
export const EXPANDABLE_TABS = [PROTOCOLS_TAB, COLLECTIONS_TAB, SAVED_SEARCHES_TAB];
