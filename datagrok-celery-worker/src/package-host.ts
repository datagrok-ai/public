/** Lazy singleton host for the Datagrok Node runtime + the published package bundle.
 *
 *  'datagrok-api' is a peer dependency loaded at runtime via require() — no static
 *  imports, so this library compiles and unit-tests without it installed. */
import {FuncCall} from './func-call';
import {Settings} from './settings';
import {logInfo} from './logger';

export class PackageHost {
  private readonly settings: Settings;
  private loadPromise: Promise<any> | null = null;

  constructor(settings: Settings) {
    this.settings = settings;
  }

  /** The DG namespace, available after the first successful init. */
  get dg(): any {
    return (globalThis as any).DG;
  }

  private ensureLoaded(call: FuncCall): Promise<any> {
    if (this.loadPromise == null) {
      const promise = (async () => {
        // eslint-disable-next-line @typescript-eslint/no-var-requires
        const runtime = require('datagrok-api/datagrok');
        const apiUrl = call.apiUrl ?? this.settings.apiUrl;
        if (apiUrl == null) {
          throw new Error('Datagrok API url is not available: ' +
            'neither call.aux.DATAGROK_API_URL nor the DATAGROK_API_URL environment variable is set');
        }
        // detached: no per-user startup data is cached — required because the worker
        // serves calls from many users (see js-api/datagrok.ts startDatagrok docs)
        await runtime.startDatagrok({apiUrl: apiUrl, detached: true});
        logInfo(`Loading package ${this.settings.packageName}` +
          `${this.settings.packageVersion ? ` (expected version ${this.settings.packageVersion})` : ''}...`);
        const pkg = await runtime.loadPackage(this.settings.packageName);
        logInfo(`Loaded package ${pkg.name} v${pkg.version}, ${pkg.functions.length} functions registered`);
        return pkg;
      })();
      // clear the cached promise on failure so the next call retries the init
      promise.catch(() => {
        if (this.loadPromise === promise)
          this.loadPromise = null;
      });
      this.loadPromise = promise;
    }
    return this.loadPromise;
  }

  /** Resolves the raw JS implementation of call.func.name from the loaded package
   *  module (loadPackage returns {name, version, webRoot, module, functions} — see
   *  js-api/src/node/package-loader.ts). NEVER routed through grok.functions.call:
   *  that would resolve the server-side DockerFunc entity and recurse into this queue. */
  async resolve(call: FuncCall): Promise<(...args: any[]) => any> {
    const pkg = await this.ensureLoaded(call);
    const name = call.funcName;
    const module = pkg.module ?? {};
    let impl = module[name];
    if (typeof impl !== 'function' && name.length > 0)
      impl = module[name.charAt(0).toLowerCase() + name.slice(1)];
    if (typeof impl !== 'function')
      throw new Error(`Function "${name}" was not found in the module exports of package "${pkg.name}"`);
    return impl;
  }

  runWithToken<T>(token: string | null, fn: () => T): T {
    if (token == null)
      return fn();
    let storage: any;
    try {
      // eslint-disable-next-line @typescript-eslint/no-var-requires
      storage = require('datagrok-api/datagrok').datagrokAsyncLocalStorage;
    }
    catch (_) {
      return fn();
    }
    return storage.run({token: token}, fn);
  }
}
