import * as dayjs from "dayjs";
import * as wu from "wu";
import { AsyncLocalStorage } from 'node:async_hooks';
import { createWindowProxy, DGNotSupportedError } from './src/node/window-proxy';
import { installDomStub } from './src/node/dom-stub';

export { DGNotSupportedError };

(globalThis as any).self = globalThis;
// grok_* registrations from the dart2js bundle land on the Proxy target;
// unregistered names degrade gracefully instead of `undefined is not a function`.
(globalThis as any).window = createWindowProxy();

// Stub asset imports (.css etc.) for unbundled CommonJS consumption; the webpack
// bundle null-loads them already, tsx/ESM consumers use their own loader hooks.
try {
    const Module = require('node:module');
    for (const ext of ['.css', '.scss', '.svg', '.png', '.jpg', '.gif', '.woff', '.woff2', '.ttf', '.eot'])
        Module._extensions[ext] ??= (m: any) => { m.exports = {}; };
} catch (_) { /* ESM-only runtime — loader hooks handle assets there */ }

(globalThis as any).fetchWithCallback = function(url: any, params: any, callback: any, errorCallback: any) {
    fetch(url, params)
        .then(async (response: any) => {
            // Response.bytes() is unavailable before Node 22 (JKG runs Node 18)
            callback(new Uint8Array(await response.arrayBuffer()), response.status, response.statusText, response.headers);
        })
        .catch(err => errorCallback(err));
};

const WebSocket = require('ws');

(globalThis as any).createWebSocket = function(url: string, headers: any) {
    return new WebSocket(url, { headers: headers['o'] });
};

// ===== AsyncLocalStorage-based token context =====
export const datagrokAsyncLocalStorage = new AsyncLocalStorage();

// Express middleware
export function tokenContextMiddleware(req: any, res: any, next: any) {
    //@ts-ignore
    const token =
        req?.headers?.authorization ||
        req?.authorization ||
        (req?.cookies && req.cookies['auth']) ||
        req.headers["x-user-api-key"] ||
        null;

    const store = { token };
    datagrokAsyncLocalStorage.run(store, () => next());
}

// Generic wrapper for non-Express usage
export function withTokenContext(req: any, callback: any) {
    const cookies = parseCookies(req.headers?.cookie);
    const token =
        req?.headers?.authorization ||
        req?.authorization ||
        cookies["auth"] ||
        req.headers["x-user-api-key"] ||
        null;

    const store = { token };
    return datagrokAsyncLocalStorage.run(store, () => callback());
}

function parseCookies(cookieHeader: string | undefined): Record<string, string> {
    const cookies: Record<string, string> = {};
    if (!cookieHeader) return cookies;

    cookieHeader.split(";").forEach((cookie: string) => {
        const [name, ...rest] = cookie.trim().split("=");
        cookies[name] = decodeURIComponent(rest.join("="));
    });

    return cookies;
}


export function getContext() {
    return datagrokAsyncLocalStorage.getStore();
}

(globalThis as any).getTokenFromAsyncLocalStorage = function() {
    const store = datagrokAsyncLocalStorage.getStore() as { token?: string } | undefined;
    return store?.token ?? null;
}

/** Initializes the Datagrok environment under Node.js.
 *
 *  - Regular mode fetches startup data (functions registry, data sources, script
 *    handlers) for the token's user — full `grok.functions` support, single-user.
 *  - `detached: true` skips startup data: functions resolve lazily over REST and
 *    server scripts run remotely. Required in shared/multi-user processes
 *    (JKG kernels, per-request server contexts) so no per-user state is cached.
 *
 *  Idempotent: a second call only re-binds `grok.dapi.token` / `grok.dapi.root`. */
export async function startDatagrok(options: {apiUrl: string, apiToken?: string | undefined, detached?: boolean}): Promise<void> {
    const g = globalThis as any;
    if (g.grok && g.DG) {
        g.grok.dapi.token = options.apiToken;
        g.grok.dapi.root = options.apiUrl;
        return;
    }
    console.log('Initializing Datagrok environment');
    //@ts-ignore
    await import ('./src/datagrok/build/web/grok_shared.dart.js');

    // dart2js runs async main() in a microtask after the import; wait for the interop
    // surface ('in' bypasses the Proxy get-fallback and answers truthfully)
    for (let waited = 0; !('grok_Init' in g.window); waited += 10) {
        if (waited > 10000)
            throw new Error('Datagrok Dart runtime failed to initialize (grok_Init never registered)');
        await new Promise(resolve => setTimeout(resolve, 10));
    }
    installDomStub();
    let DG = await import('./node-api');
    // ui builders work against the DOM stub (inert elements); Dart-backed widgets
    // still throw DGNotSupportedError via the window Proxy
    let ui = await import('./ui');

    g.grok = DG.grok;
    g.DG = DG;
    g.ui = ui;
    const api: any = (typeof window !== 'undefined' ? window : global.window) as any;
    api.dayjs = dayjs;
    api.wu = wu;
    api.grok = DG.grok;
    api.DG = DG;
    g.wu = wu;
    g.dayjs = dayjs;


    DG.grok.dapi.token = options.apiToken;
    DG.grok.dapi.root = options.apiUrl;

    let status = await g.window.grok_Init(options.detached === true);
    if (status != 0)
        throw new Error(`Datagrok initialization failed (grok_Init returned ${status}, see log above)`);
    console.log('Datagrok environment is ready');
}
