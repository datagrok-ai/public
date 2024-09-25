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

  let isLost = false;
  if (gpuDevice) {
    gpuDevice.lost.then(() => {
      isLost = true;
    });
    await new Promise((r) => setTimeout(r, 10)); // wait to see if the device is lost
  }
  if (!gpuDevice || isLost) {
    const requiredBufferSize = 1_000_000_000; // ~1000MB
    const adapterLimits = gpuAdapter.limits;
    const buffferSizeLimit = adapterLimits.maxBufferSize;
    const storageBufferSizeLimit = adapterLimits.maxStorageBufferBindingSize;
    try {
      gpuDevice = await gpuAdapter.requestDevice({requiredLimits: {
        maxBufferSize: Math.min(buffferSizeLimit, requiredBufferSize),
        maxStorageBufferBindingSize: Math.min(storageBufferSizeLimit, requiredBufferSize)
      }});
      return gpuDevice;
    } catch (e) {
      console.error('Failed to create device with required limits', e);
      gpuDevice = await gpuAdapter.requestDevice();
      return gpuDevice;
    }
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
  const info = await gpuAdapter.requestAdapterInfo();
  if (!info) return null;

  const outString = replaceEmptyString(
    info.description, replaceEmptyString(info.vendor, 'No GPU description available'));
  return outString;
}

function replaceEmptyString(str: string, replacement: string) {
  return !str || str == '' ? replacement : str;
}

