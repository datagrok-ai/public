/*
 * WebGPU, opened in the background when PMX Studio opens, so that a calculation that wants it finds the device and
 * the shader compiler ready (the first pipeline costs ~0.13 s, the device ~15 ms). The default adapter: asking for
 * high performance would wake a laptop's discrete GPU with nothing to compute. Only in a secure context (https or
 * localhost) of a browser with WebGPU (Chrome / Edge 113+, Safari 26, Firefox 141+ on Windows); elsewhere null, and
 * the calculations stay on the CPU (wasm). A lost device (a driver reset, sleep) is forgotten: the next call opens
 * another. The digital twin asks it whether to run on WebGPU (Crux's module worker opens its own device; the
 * shader compiler this one warmed is the browser's).
 */
/* eslint-disable @typescript-eslint/no-explicit-any */

let device: Promise<any | null> | null = null;

/** The WebGPU device, or null where there is none; the first call opens it. */
export function gpu(): Promise<any | null> {
  return device ??= open();
}

async function open(): Promise<any | null> {
  const g = (navigator as any).gpu;
  if (!self.isSecureContext || !g)
    return null;
  try {
    // the first requestAdapter in a fresh headless Chrome can come back null: ask twice (as Crux's lab/gpu.ts)
    const adapter = (await g.requestAdapter()) ?? (await g.requestAdapter());
    if (!adapter)
      return null;
    // the adapter's largest buffers: the defaults (256 MB, 128 MB a binding) are below what a large grid needs
    const {maxBufferSize, maxStorageBufferBindingSize} = adapter.limits;
    const d = await adapter.requestDevice({requiredLimits: {maxBufferSize, maxStorageBufferBindingSize}});
    d.lost.then(() => device = null);
    // the first pipeline warms the shader compiler: 130 ms, the next ones about 12 (measured, Chrome on RDNA 3)
    await d.createComputePipelineAsync({layout: 'auto', compute: {entryPoint: 'main', module: d.createShaderModule({
      code: '@group(0) @binding(0) var<storage, read_write> x: array<f32>;\n' +
        '@compute @workgroup_size(64) fn main(@builtin(global_invocation_id) i: vec3u) { x[i.x] = exp(x[i.x]); }'})}});
    return d;
  } catch {
    return null;
  }
}
