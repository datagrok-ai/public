# NodeJS API Demo

A minimal package demonstrating the **js-api under Node.js** (GROK-20443): a `nodejs`
script that uses `grok`/`DG` server-side and consumes **package functions** — no browser
anywhere in the chain, and **no bootstrap code in the script**: the platform injects the
js-api globals (`grok`, `DG`, `loadPackage`) with the calling user's token.

## What `NodeDemo` shows

`scripts/node-demo.js` (`NodejsApiDemo:NodeDemo`, `#language: nodejs`) runs in the
Jupyter Kernel Gateway's Node kernel and:

1. calls `grok.dapi.users.current()` — server API with the caller's identity;
2. clones the input `DG.DataFrame` and adds a computed column (`colName` doubled);
3. passes the dataframe to **this package's Python script** (`NodejsApiDemo:PyMean`)
   via `grok.functions.call` — cross-language package function call;
4. loads the **ApiTests package bundle** with `loadPackage()` and calls one of its
   **JS functions** in-process (`ApiTests:dummyPackageFunction`);
5. reads this package's AppData file (`files/hello.txt`) via `grok.dapi.files`.

## Run

```bash
grok publish <host>       # from this folder
```

then call it from anywhere — e.g. a standalone Node.js client:

```js
const {startDatagrok} = require('datagrok-api/datagrok');
await startDatagrok({apiUrl, apiToken, detached: true});
const df = DG.DataFrame.fromCsv('v\n1\n2\n3');
const r = await grok.functions.call('NodejsApiDemo:NodeDemo', {df, colName: 'v'});
// {user, enriched, mean, jsSum, fileText}
```

or from the browser console with the same two lines minus the bootstrap.

## Requirements

The JKG image must ship `datagrok-api` (built from repo source since GROK-20443;
`deploy/jupyter_kernel_gateway/Dockerfile`). Scripts run with `detached` js-api
init — nothing user-specific is cached across calls in the shared kernel.
