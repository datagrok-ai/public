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
  if (!gpuAdapter) {
    // reason: only here we get the gpuAdapter
    // eslint-disable-next-line no-restricted-syntax
    gpuAdapter = await navigator.gpu.requestAdapter();
    if (gpuAdapter == null)
      return null;
  }

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

