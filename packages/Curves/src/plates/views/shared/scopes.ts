export const MAPPING_SCOPES = {
  TEMPLATE: 'template',
  DRC: 'drc',
  DOSE_RATIO: 'doseRatio'
} as const;

export type MappingScope = typeof MAPPING_SCOPES[keyof typeof MAPPING_SCOPES];