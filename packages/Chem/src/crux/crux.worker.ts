//@ts-ignore
import initCrux, {Collection, CollectionBuilder} from './crux_wasm.js';

const ctx: Worker = self as any;
const collections = new Map<string, Collection>();

/** Builds the collection of one index segment; returns the (segment-local) rows crux could not parse. */
function build(key: string, smiles: string[]): Uint32Array {
  const builder = new CollectionBuilder(true);
  const failed: number[] = [];
  // one molecule per call: the builder only counts failures, the NOT CONTAINS complement needs their rows
  for (let i = 0; i < smiles.length; i++) {
    if (builder.addMany([smiles[i]]) === 0)
      failed.push(i);
  }
  collections.set(key, builder.finish());
  return new Uint32Array(failed);
}

function drop(keys: string[]): void {
  for (const key of keys) {
    collections.get(key)?.free();
    collections.delete(key);
  }
}

ctx.addEventListener('message', async (e: MessageEvent) => {
  const {op, args} = e.data;
  const port = e.ports[0];
  try {
    let result;
    if (op === 'module::init')
      await initCrux({module_or_path: args[0]});
    else if (op === 'build')
      result = build(args[0], args[1]);
    else if (op === 'search')
      result = collections.get(args[0])!.substructureSearch(args[1], 0);
    else if (op === 'drop')
      drop(args[0]);
    else
      throw new Error(`Unknown operation: ${op}`);
    port.postMessage({op: op, retval: result});
  } catch (err: any) {
    port.postMessage({error: err instanceof Error ? err.message : String(err)});
  }
});
