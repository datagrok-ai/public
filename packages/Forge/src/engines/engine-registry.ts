import * as DG from 'datagrok-api/dg';
import {ENGINE_ROLES, Engine, EngineKind} from './engine';

export class EngineRegistry {
  static discover(): Engine[] {
    const engines = new Map<string, Engine>();
    // ENGINE_ROLES starts with 'train', so an engine is described by its train function when it has one.
    for (const role of ENGINE_ROLES) {
      for (const func of DG.Func.find({meta: {mlrole: role}})) {
        const name: unknown = func.options['mlname'];
        if (typeof name !== 'string' || name === '')
          continue;
        let engine = engines.get(name);
        if (engine === undefined) {
          engine = EngineRegistry.create(name, func);
          engines.set(name, engine);
        }
        engine.functions[role] = func;
        if (role === 'isInteractive')
          engine.isLiveUpdate = func.options['mlupdate'] !== 'false';
      }
    }
    return [...engines.values()].sort((a, b) =>
      a.kind === b.kind ? a.name.localeCompare(b.name) : a.kind === 'function' ? -1 : 1);
  }

  private static create(name: string, func: DG.Func): Engine {
    const kind: EngineKind = func instanceof DG.Script ? 'script' : 'function';
    const pkg: unknown = func.package;
    const namespace = kind === 'function' && pkg instanceof DG.Package ? pkg.name : func.nqName.split(':')[0];
    return {name, namespace, kind, functions: {}, isLiveUpdate: true};
  }
}
