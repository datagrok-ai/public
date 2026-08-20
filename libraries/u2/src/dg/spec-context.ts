/* The platform SpecContext: the `cmd:` platform-function tier backed by grok.functions.call,
   with failures ballooned once at this boundary — the one handler that owns them (Q6). */
import * as grok from 'datagrok-api/grok';
import {SpecContext} from '../spec/spec.js';
import type {SpecContextOptions} from '../spec/spec.js';

export function platformContext(options: SpecContextOptions = {}): SpecContext {
  return new SpecContext({
    ...options,
    callFunction: options.callFunction ?? (async (name, args) => {
      try {
        return await grok.functions.call(name, args);
      } catch (e) {
        grok.shell.error(e instanceof Error ? e.message : String(e));
        return undefined;
      }
    }),
  });
}
