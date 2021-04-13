var grokCache = {}


async function grokCacheModuleInit() {
  const htmlMStatus = document.querySelector('#gk_module_status');

  if (window.createGrokCache) {
    grokCache = await createGrokCache();
  } else {
    let err = "GrokCache WASM Module is not available!";
    htmlMStatus.textContent = err;
    console.error(err);
  }
}

// Test entry point
//
function grokCacheTestStart() {
  grokCacheModuleInit().then(() => {
    {
      const htmlMStatus = document.querySelector('#gk_module_status');
      let greeting = grokCache.wasm_version();
      htmlMStatus.textContent = greeting;
    }
  });
}