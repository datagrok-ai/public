export async function getGPUDevice() {
  const adapter = await navigator.gpu.requestAdapter({powerPreference: 'high-performance'});
  if (adapter == null)
    return null;

  const requiredBufferSize = 1_000_000_000; // ~1000MB
  const adapterLimits = adapter.limits;
  const buffferSizeLimit = adapterLimits.maxBufferSize;
  const storageBufferSizeLimit = adapterLimits.maxStorageBufferBindingSize;
  try {
    const device = await adapter.requestDevice({requiredLimits: {
      maxBufferSize: Math.min(buffferSizeLimit, requiredBufferSize),
      maxStorageBufferBindingSize: Math.min(storageBufferSizeLimit, requiredBufferSize)
    }});
    return device;
  } catch (e) {
    console.error('Failed to create device with required limits', e);
    const device = await adapter.requestDevice();
    return device;
  }
}
