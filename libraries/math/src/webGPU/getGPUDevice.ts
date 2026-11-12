let gpuAdapter: GPUAdapter | null = null;
let gpuDevice: GPUDevice | null = null;


export async function getGPUDevice() {
  if (!navigator.gpu) {
    console.error('WebGPU is not supported in this browser');
    return null;
  }
  if (!gpuAdapter) {
    //reason: only here we get the gpuAdapter
    // eslint-disable-next-line no-restricted-syntax
    gpuAdapter = await navigator.gpu.requestAdapter({powerPreference: 'high-performance'});
    if (gpuAdapter == null)
      return null;
  }

  if (!gpuDevice) {
    const requiredBufferSize = 1_000_000_000; // ~1000MB
    const adapterLimits = gpuAdapter.limits;
    const buffferSizeLimit = adapterLimits.maxBufferSize;
    const storageBufferSizeLimit = adapterLimits.maxStorageBufferBindingSize;
    try {
      gpuDevice = await gpuAdapter.requestDevice({requiredLimits: {
        maxBufferSize: Math.min(buffferSizeLimit, requiredBufferSize),
        maxStorageBufferBindingSize: Math.min(storageBufferSizeLimit, requiredBufferSize)
      }});
    } catch (e) {
      console.error('Failed to create device with required limits', e);
      gpuDevice = await gpuAdapter.requestDevice();
    }
    // a lost device is replaced on the next call; an adapter gives out one device only, so it goes too
    const device = gpuDevice;
    device.lost.then(() => {
      if (gpuDevice === device) {
        gpuDevice = null;
        gpuAdapter = null;
      }
    });
  }
  return gpuDevice;
}

export async function getGPUAdapterDescription() {
  if (!navigator.gpu) {
    console.error('WebGPU is not supported in this browser');
    return null;
  }
  // the adapter getGPUDevice computes on: requesting another here would leave getGPUDevice the wrong one
  if (!gpuAdapter && !(await getGPUDevice()))
    return null;

  let info: GPUAdapterInfo | null = null;
  if ('info' in (gpuAdapter as any))
    info = (gpuAdapter as any).info as GPUAdapterInfo;
  // this option is sort of deprecated but still available in every initial release
  // else if ('requestAdapterInfo' in gpuAdapter && typeof gpuAdapter.requestAdapterInfo === 'function')
  //   info = (await gpuAdapter.requestAdapterInfo()) as GPUAdapterInfo;
  if (!info) return 'No GPU description available';

  const outString = replaceEmptyString(
    info.description, replaceEmptyString(info.vendor, 'No GPU description available'));
  return outString;
}

function replaceEmptyString(str: string, replacement: string) {
  return !str || str == '' ? replacement : str;
}

/** The GPU the browser draws and computes on, as WebGL names it (Chrome on Windows:
 * "ANGLE (AMD, AMD Radeon(TM) 860M Graphics (0x00001114) Direct3D11 vs_5_0 ps_5_0, D3D11)"); null if hidden. */
export function getGPURenderer(): string | null {
  const gl = document.createElement('canvas').getContext('webgl');
  if (!gl)
    return null;
  const ext = gl.getExtension('WEBGL_debug_renderer_info');
  const renderer = ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) as string : null;
  gl.getExtension('WEBGL_lose_context')?.loseContext();
  return renderer;
}

/** Whether a renderer string (see getGPURenderer) names a GPU built into the CPU: Intel's but Arc A/B cards,
 * AMD's Radeon "Graphics" without RX, Pro or Instinct. */
// ponytail: name heuristic (no web API tells integrated from discrete); add names here as they show up
export function isIntegratedGPU(renderer: string): boolean {
  if (/intel/i.test(renderer))
    return !/arc\(tm\)\s*[ab]\d{3}/i.test(renderer);
  if (/amd|radeon/i.test(renderer))
    return /graphics/i.test(renderer) && !/\brx\b|\bpro\b|instinct/i.test(renderer);
  return false;
}

