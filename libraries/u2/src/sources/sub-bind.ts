/* A prop a source drives from one of the node's binds, resolved at start() — by then the input it
   names may be declared after the source (entity-ref's `id`, entity-source's `filter`). */
import type {Signal} from '../core/signals.js';
import type {ComponentEnv} from '../spec/registry.js';

export function subBind(env: ComponentEnv, key: string): Signal<unknown> | null {
  const path = env.subBinds[key];
  return path === undefined ? null : env.resolve(path);
}
