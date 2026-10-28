/* Who owns what u2 builds: a scope, a component, or — given none — whoever is building. Every API with an `owner`
   resolves it here, so the rule and its error read the same everywhere. */
import {Scope, Usage} from './scope.js';
import {Component} from './component.js';

/** What an `owner` takes: a scope, or a component (its scope). */
export type Owner = Scope | Component;

export class Ownership {
  /** The scope of `owner`, else of whoever is building; throws, naming `what`, when there is none. */
  static scopeOf(owner: Owner | undefined, what: string): Scope {
    const scope = Ownership._resolve(owner);
    if (scope === undefined) {
      Usage.fail('owner-missing', `${what} needs an owner — pass one, or call it under \`owner.build(…)\` or ` +
        '`component((self) => …)`; after an `await`, inside `owner.task(async (t) => … t.build(…))`');
    }
    return scope;
  }

  /** Runs `dispose` when `owner` goes — else whoever is building; nobody is reported (`owner-missing`) — and answers
   * the release that takes it off again, for a thing disposed before its owner: a long-lived owner keeps no closure
   * per thing it outlived. */
  static adopt(owner: Owner | undefined, dispose: () => void): () => void {
    const scope = Ownership._resolve(owner);
    if (scope === undefined) {
      Usage.report('owner-missing', 'built with no `owner` while nothing was building, so nothing will dispose it — ' +
        'pass `owner`, or build it under `owner.build(…)` or `component((self) => …)`; after an `await`, inside ' +
        '`owner.task(async (t) => … t.build(…))`; a new root on purpose, under `Scope.root(\'<what>\')`');
      return () => {};
    }
    scope.own(dispose);
    return () => scope.disown(dispose);
  }

  private static _resolve(owner: Owner | undefined): Scope | undefined {
    return owner === undefined ? Scope.ambient : owner instanceof Scope ? owner : owner.scope;
  }
}
