# GrokCache

GrokCache is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.

GrokCache implements a RAM cache for resources based on WASM for compression 
and memory isolation safety.



## Stash

Things not needed but may be useful in other cases:

package.json
-------
"dependencies": {
    "fs": "^0.0.1-security", // <---  
    "worker-loader": "^3.0.8" // <---- not needed
}

"fs" is not needed when emscripten compiles: -s ENVIRONMENT="web"


webpack.config.js
------------

Useful for cases when WebAsm is supposed to run in a Worker?

experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true
}
