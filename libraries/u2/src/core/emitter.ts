import type {ObservableLike} from './widget-like.js';

export class Emitter<T = void> implements ObservableLike<T> {
  private readonly _listeners = new Set<(x: T) => void>();

  subscribe(next: (x: T) => void): {unsubscribe(): void} {
    this._listeners.add(next);
    return {unsubscribe: () => this._listeners.delete(next)};
  }

  fire(x: T): void {
    for (const listener of [...this._listeners])
      listener(x);
  }

  get count(): number {
    return this._listeners.size;
  }

  clear(): void {
    this._listeners.clear();
  }
}
