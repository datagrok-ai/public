# NodeJS API Demo

A minimal package demonstrating the **js-api under Node.js** (GROK-20443): a `nodejs`
script that uses `grok`/`DG` server-side and consumes **package functions** — no browser
anywhere in the chain.

## What it shows

`scripts/node-demo.js` (`NodejsApiDemo:NodeDemo`, `#language: nodejs`) runs in the
Jupyter Kernel Gateway's Node kernel and:

1. bootstraps the js-api once per kernel (`startDatagrok({..., detached: true})`) with the
   calling user's token re-bound per call;
2. calls `grok.dapi.users.current()` — server API from inside the script;
3. builds a `DG.DataFrame` and passes it to **this package's Python script**
   (`NodejsApiDemo:PyMean`) via `grok.functions.call` — a cross-language package
   function call with dataframe transfer;
4. reads this package's AppData file (`files/hello.txt`) via `grok.dapi.files`.

Expected outputs: `user` (caller's name), `cells = 4`, `mean = 2.5`,
`fileText = "hello from package files"`.

## Run

```bash
grok publish <host>       # from this folder
```

then call it from anywhere — for example from a standalone Node.js client:

```js
const {startDatagrok} = require('datagrok-api/datagrok');
await startDatagrok({apiUrl, apiToken, detached: true});
const r = await grok.functions.call('NodejsApiDemo:NodeDemo', {});
// {user, cells, mean, fileText}
```

## Notes

- The script `require`s `datagrok-api/datagrok` and falls back to `/tmp/datagrok-api` —
  the fallback covers JKG images that don't ship `datagrok-api` globally yet
  (the Dockerfile installs it behind the `DATAGROK_API_VERSION` build arg).
- `detached: true` is required in shared kernels: it skips per-user startup data, so
  nothing user-specific is cached across calls; functions resolve lazily over REST.
- Verified end-to-end on sandbox 2026-07-20 (see GROK-20443).
