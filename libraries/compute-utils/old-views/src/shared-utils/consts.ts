
export const RESTRICTED_PATH = 'restrictedValues' as const;
export const EDIT_STATE_PATH = 'editState' as const;
export enum DIRECTION {
  INPUT = 'Input',
  OUTPUT = 'Output'
};

export enum VISIBILITY_STATE {
  HIDDEN = 'hidden',
  VISIBLE = 'visible',
}

export enum ABILITY_STATE {
  ENABLED = 'enabled',
  DISABLED = 'disabled',
}

export type INPUT_STATE = 'disabled' | 'restricted' | 'restricted unlocked' | 'inconsistent' | 'user input';

export type VIEW_STATE = 'inconsistent' | 'consistent';

export enum SYNC_FIELD {
  INPUTS = 'inputs',
  OUTPUTS = 'outputs'
}

export interface ValidationRequestPayload {
  field?: string,
  isRevalidation: boolean,
  isNewOutput?: boolean,
  context?: any,
}

export type SyncFields = SYNC_FIELD.INPUTS | SYNC_FIELD.OUTPUTS;
export const syncParams = {
  [SYNC_FIELD.INPUTS]: 'inputParams',
  [SYNC_FIELD.OUTPUTS]: 'outputParams',
} as const;
