import {Scope} from './scope.js';

/** Base of every u2 component: a visual root plus an effect/cleanup scope.
 * Platform-free; in Datagrok, `u2/dg`'s `host()` wires {@link dispose} into the existing
 * `DG.Widget` kill channel. Standalone hosts (gallery, tests) call {@link dispose} directly. */
export class Component {
  readonly root: HTMLElement;
  readonly scope: Scope;

  constructor(root?: HTMLElement, scope?: Scope) {
    this.root = root ?? document.createElement('div');
    this.root.classList.add('u2-component');
    this.scope = scope ?? new Scope();
    const owner = Scope.ambient;
    if (owner && owner !== this.scope)
      owner.own(() => this.dispose());
  }

  /** Builds a component out of plain elements: `builder` runs with the new component's scope
   * ambient, so signal bindings and nested components inside it are owned by the result.
   * The ambient scope only spans the synchronous part of `builder` — components constructed after
   * an `await` inside it are adopted by nobody and must be disposed by hand, or constructed before
   * the first await. */
  static build(builder: () => HTMLElement | (HTMLElement | Component | string)[]): Component {
    const scope = new Scope();
    const built = Scope.runWith(scope, builder);
    if (!Array.isArray(built))
      return new Component(built, scope);
    const root = document.createElement('div');
    for (const child of built)
      root.append(child instanceof Component ? child.root : child);
    return new Component(root, scope);
  }

  run<T>(fn: () => T): T {
    return Scope.runWith(this.scope, fn);
  }

  effect(fn: () => void): this {
    this.scope.effect(fn);
    return this;
  }

  own(dispose: () => void): this {
    this.scope.own(dispose);
    return this;
  }

  dispose(): void {
    this.scope.dispose();
  }
}
