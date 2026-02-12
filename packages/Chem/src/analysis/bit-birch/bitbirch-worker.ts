//@ts-ignore
import initBBWasm, {WasmBitBirch} from '../../bbwasm.js';

async function ensureBBWasmInit(webroot: string): Promise<void> {
  await initBBWasm(`${webroot}/dist/bbwasm_bg.wasm`);
}

onmessage = async (event: MessageEvent) => {
  const {threshold, branchingFactor, nFeatures, fps, fpCount, webRoot}: {
    threshold: number,
    branchingFactor: number,
    nFeatures: number,
    fps: Uint8Array,
    fpCount: number,
    webRoot: string,
    } = event.data;
  await ensureBBWasmInit(webRoot);

  const bb = new WasmBitBirch(threshold, branchingFactor, nFeatures);
  try {
    bb.fit(fps, fpCount);
    const assignments: Int32Array = bb.get_assignments();
    postMessage({assignments});
    // Map back to full column, leaving nulls for invalid rows
  } catch (e) {
    postMessage({error: e instanceof Error ? e.message : String(e)});
  } finally {
    bb.free();
  }
};
