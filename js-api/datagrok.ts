import * as dayjs from "dayjs";
import * as wu from "wu";
const { AsyncLocalStorage } = require('async_hooks');

(globalThis as any).self = globalThis;
(globalThis as any).window = {};

(globalThis as any).fetchWithCallback = function(url:any, params:any, callback:any, errorCallback:any) {
    fetch(url, params['o'])
        .then(async(response: any) => {callback(await response.bytes(), response.status, response.statusText, response.headers);})
        .catch(err => errorCallback(err));
};

const WebSocket = require('ws');

(globalThis as any).createWebSocket = function(url: string, headers: any) {
    return new WebSocket(url, { headers: headers['o'] });
};

// ===== AsyncLocalStorage-based token context =====
export const datagrokAsyncLocalStorage = new AsyncLocalStorage();

// Express middleware
export function tokenContextMiddleware(req, res, next) {
    //@ts-ignore
    const token =
        req?.headers?.authorization ||
        req?.authorization ||
        (req?.cookies && req.cookies['auth']) ||
        null;

    const store = { token };
    datagrokAsyncLocalStorage.run(store, () => next());
}

// Generic wrapper for non-Express usage
export function withTokenContext(req, callback) {
    const token =
        req?.headers?.authorization ||
        req?.authorization ||
        parseCookies(req.headers?.cookie)['auth'] ||
        null;

    const store = { token };
    return datagrokAsyncLocalStorage.run(store, () => callback());
}

function parseCookies(cookieHeader) {
    const cookies = {};
    if (!cookieHeader) return cookies;
    cookieHeader.split(';').forEach(cookie => {
        const [name, ...rest] = cookie.trim().split('=');
        cookies[name] = decodeURIComponent(rest.join('='));
    });
    return cookies;
}


export function getContext() {
    return datagrokAsyncLocalStorage.getStore();
}

(globalThis as any).getTokenFromAsyncLocalStorage = function() {
    const store = datagrokAsyncLocalStorage.getStore();
    return store ? store.token : null;
}

export async function startDatagrok(options: {apiUrl: string, apiToken?: string | undefined}): Promise<any> {
    //@ts-ignore
    await import ('./src/datagrok/web/grok_shared.dart.js');

    await new Promise(resolve => setTimeout(resolve, 100));
    (globalThis as any).document = {createElement: () => {return {bind: {}}}};
    let DG = await import('./node-api');

    (globalThis as any).grok = DG.grok;
    (globalThis as any).DG = DG;
    const api: any = (typeof window !== 'undefined' ? window : global.window) as any;
    api.dayjs = dayjs;
    api.wu = wu;
    api.grok = DG.grok;
    api.DG = DG;
    (globalThis as any).wu = wu;
    (globalThis as any).dayjs = dayjs;


    DG.grok.dapi.token = options.apiToken;
    DG.grok.dapi.root = options.apiUrl;

    let status = await (globalThis as any).window.grok_Init();
    console.log(status);
}