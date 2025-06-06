export enum CDDVaultSearchType {
    EXACT_MATCH = 'exact',
    SIMILARITY = 'similarity',
    SUBSTRUCTURE = 'substructure',
}

export enum CDDVaultInNotInCond {
    IN = 'in',
    NOT_IN = 'not in',
}

export enum ProtocolCond {
    PROTOCOL = 'Specific protocol',
    FIELD = 'Protocol field',
}

export enum ProtocolRunCond {
    ANY = '(any run)',
    SPECIFIC = 'specific run',
}

export const CDD_SEARCH_TYPES = [CDDVaultSearchType.SUBSTRUCTURE, CDDVaultSearchType.SIMILARITY, CDDVaultSearchType.EXACT_MATCH];

export const IN_NOT_IN_COND_CHOICES = [CDDVaultInNotInCond.IN, CDDVaultInNotInCond.NOT_IN];

export const PROTOCOL_COND_CHOICES = [ProtocolCond.PROTOCOL, ProtocolCond.FIELD];

export const PROTOCOL_RUN_COND_CHOICES = [ProtocolRunCond.ANY, ProtocolRunCond.SPECIFIC];

export const ANY_READOUT_CHOICE = '(any readout definition)';

export const ANY_RUN_CHOICE = '(any run)';

export const PROTOCOLS_TAB = 'Protocols';
export const COLLECTIONS_TAB = 'Collections';
export const SAVED_SEARCHES_TAB = 'Saved Searches';
export const MOLECULES_TAB = 'Molecules';
export const SEARCH_TAB = 'Search';

export const ALL_TABS = [PROTOCOLS_TAB, COLLECTIONS_TAB, SAVED_SEARCHES_TAB, MOLECULES_TAB, SEARCH_TAB];
export const EXPANDABLE_TABS = [PROTOCOLS_TAB, COLLECTIONS_TAB, SAVED_SEARCHES_TAB];