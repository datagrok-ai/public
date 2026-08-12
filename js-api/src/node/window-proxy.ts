/**
 * Node.js `window` shim.
 *
 * The dart2js bundle (grok_shared.dart.js) registers its interop functions by
 * setting `grok_*` properties on `globalThis.window` — those land on the Proxy
 * target and win over the fallback. Reading a `grok_*` name the Dart runtime did
 * NOT register returns a stub that either maps to a console no-op (curated list)
 * or throws a descriptive {@link DGNotSupportedError} instead of the raw
 * `api.grok_X is not a function` TypeError.
 *
 * The `has` trap is left untouched on purpose: `'grok_X' in window` (and Dart-side
 * hasProperty checks) still answer truthfully whether a function is registered.
 */

export class DGNotSupportedError extends Error {
  constructor(name: string) {
    super(`${name} is not available under Node.js: it is implemented by the Datagrok browser client (UI/Dart). ` +
      'Use the server-side equivalents (grok.dapi.*, grok.functions.*, DG.DataFrame) or run this code in the browser.');
    this.name = 'DGNotSupportedError';
  }
}

function balloon(msg: any, type?: string): void {
  const log = type === 'error' ? console.error : type === 'warning' ? console.warn : console.log;
  log(`[shell.${type ?? 'info'}] ${typeof msg === 'string' ? msg : msg?.textContent ?? JSON.stringify(msg)}`);
}

/** UI notifications that degrade to console output / no-ops instead of throwing. */
const consoleStubs: Record<string, Function> = {
  grok_Balloon: balloon,
  grok_ShowHelp: () => {},
  grok_Utils_LoadJsCss: () => Promise.resolve([]),
};

/** Names that js-api feature-detects with `api.grok_X == null` and has a fallback
 * for — they must read as undefined when unregistered, not as a throwing stub. */
const nullCheckedNames = new Set(['grok_View_Set_Name', 'grok_UI_Render']);

export function createWindowProxy(): any {
  const target: any = {
    addEventListener: () => {},
    removeEventListener: () => {},
    dispatchEvent: () => true,
    navigator: {userAgent: `Node.js ${typeof process !== 'undefined' ? process.version : ''}`},
    location: {href: '', origin: '', pathname: '', search: '', hash: ''},
    /** Feature detection for API authors: true if the Dart runtime registered the function. */
    isDartRegistered: (name: string) => name in target,
  };
  return new Proxy(target, {
    get(t, prop, receiver) {
      if (prop in t)
        return Reflect.get(t, prop, receiver);
      if (typeof prop === 'string' && prop.startsWith('grok_') && !nullCheckedNames.has(prop)) {
        // grok_UI_* are element builders — per the hybrid degradation contract they
        // return inert stub elements instead of throwing, so shared layout code runs
        if (prop.startsWith('grok_UI_'))
          return (...args: any[]) => {
            const el: any = (globalThis as any).document?.createElement?.('div') ?? {};
            if (typeof args[0] === 'string') el.textContent = args[0];
            return el;
          };
        return consoleStubs[prop] ?? (() => { throw new DGNotSupportedError(prop as string); });
      }
      return undefined;
    },
  });
}
