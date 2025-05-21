
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

export async function startDatagrok(): Promise<any> {
    await import ('./src/datagrok/web/grok_shared.dart.js');

    await new Promise(resolve => setTimeout(resolve, 100));
    (globalThis as any).document = {createElement: () => {return {bind: {}}}};


    await import('./node-api');
}


export type {DG as _DG, grok as _grok} from './node-api';

