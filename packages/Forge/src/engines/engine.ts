import * as DG from 'datagrok-api/dg';

export const ENGINE_ROLES = ['train', 'apply', 'isApplicable', 'isInteractive', 'visualize'] as const;
export type EngineRole = typeof ENGINE_ROLES[number];
export type EngineKind = 'function' | 'script';
const TRAIN_DATA_INPUTS = ['df', 'predictColumn'] as const;

export interface Engine {
  name: string;
  namespace: string;
  kind: EngineKind;
  functions: Partial<Record<EngineRole, DG.Func>>;
  isLiveUpdate: boolean;
}

export function isComplete(engine: Engine): boolean {
  const {train, apply, isApplicable} = engine.functions;
  return train !== undefined && apply !== undefined && isApplicable !== undefined;
}

export function hyperparametersOf(engine: Engine): DG.Property[] {
  const inputs = engine.functions.train?.inputs ?? [];
  return inputs.filter((p) => !TRAIN_DATA_INPUTS.some((name) => name === p.name));
}

export function rolesOf(engine: Engine): EngineRole[] {
  return ENGINE_ROLES.filter((role) => engine.functions[role] !== undefined);
}
