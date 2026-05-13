# WASM in Datagrok Plugins

For plugins with WASM compiled specifically for that plugin (C++ via Emscripten, Rust via wasm-pack).
**Reference implementations**: `Chem` (RDKit, C++/Emscripten), `EDA` (C++/Emscripten).

## Directory Layout

```
PluginName/
├── wasm/                        # Compiled WASM artifacts (check in to git)
│   ├── plugin.js                # Emscripten/wasm-pack JS glue
│   └── plugin_bg.wasm           # Binary WASM module
├── wasm-src/                    # WASM source (C++/Rust) — build output gitignored
│   ├── CMakeLists.txt or Cargo.toml
│   └── src/
├── src/
│   ├── package.ts
│   └── wasm-loader.ts           # Module wrapper with lazy init
├── webpack.config.js
└── package.json
```

## webpack.config.js

```js
module.exports = {
  experiments: {
    asyncWebAssembly: true,
    topLevelAwait: true,
  },
  resolve: {
    extensions: ['.wasm', '.ts', '.mjs', '.js', '.json'],
  },
  module: {
    rules: [
      {
        test: /\.wasm$/i,
        type: 'javascript/auto',
        loader: 'file-loader',
        options: {
          publicPath: 'dist/',
          name: '[name].[ext]',
        },
      },
    ],
  },
};
```

## package.json

```json
{
  "sources": [
    "wasm/plugin.js",
    "wasm/plugin_bg.wasm"
  ],
  "browserFeatures": ["wasm"]
}
```

`sources` tells the deploy tool to include these files alongside the webpack bundle.

## Lazy Module Init (src/wasm-loader.ts)

```ts
let _module: any = null;

export async function getModule() {
  if (_module) return _module;
  // file-loader rewrites this import to the dist/ URL at build time
  const init = await import('../wasm/plugin.js');
  _module = await init.default();
  return _module;
}
```

Use in package functions:

```ts
import { getModule } from './wasm-loader';

export async function myFunction(col: DG.Column): Promise<DG.DataFrame> {
  const wasm = await getModule();   // one-time init, cached thereafter
  const result = wasm.process(col.getRawData());
  return DG.DataFrame.fromColumns([...]);
}
```

## Build Pipeline

Add to `package.json` scripts:

```json
{
  "scripts": {
    "build:wasm": "emcc ... -o wasm/plugin.js",
    "build": "npm run build:wasm && grok api && grok check --soft && webpack"
  }
}
```

## Key Rules

- Check compiled `wasm/` artifacts into git — keeps deployment self-contained.
- Gitignore `wasm-src/` build output (object files, intermediate artifacts).
- Always lazy-init the module and cache the instance — WASM init is a one-time cost.
- Access the `.wasm` file at runtime via `_package.webRoot + 'dist/plugin_bg.wasm'` or via the JS glue import (file-loader handles the URL rewrite).
- For WASM from npm packages (e.g. `parquet-wasm`, `@biowasm/aioli`), no special setup is needed — see **Arrow** and **Bio** for examples.
