import {Optimizer} from './optimizer';
import type {CommonSettings} from './types';

type OptimizerFactory = () => Optimizer<any>;

const registry = new Map<string, OptimizerFactory>();

/** Register an optimizer under a case-insensitive key. */
export function registerOptimizer(name: string, factory: OptimizerFactory): void {
  registry.set(name.toLowerCase(), factory);
}

/** Retrieve an optimizer by name. Throws if unknown. */
export function getOptimizer<S extends CommonSettings>(name: string): Optimizer<S> {
  const factory = registry.get(name.toLowerCase());
  if (!factory)
    throw new Error(`Unknown optimizer "${name}". Available: ${listOptimizers().join(', ')}`);
  return factory() as Optimizer<S>;
}

/** List registered optimizer names. */
export function listOptimizers(): string[] {
  return [...registry.keys()];
}
