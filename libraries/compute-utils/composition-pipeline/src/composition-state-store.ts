import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject, Observable, merge} from 'rxjs';
import {filter} from 'rxjs/operators';

export class CompositionStateStore {
  private state: {[name: string]: BehaviorSubject<any | undefined>} = {};
  private propagate = true;

  constructor(keys: string[] = []) {
    for (const key of keys)
      this.state[key] = new BehaviorSubject(undefined);
  }

  setItem(key: string, value: any, propagate = true) {
    if (this.state[key]) {
      try {
        this.propagate = propagate;
        this.state[key].next(value);
      } finally {
        this.propagate = true;
      }
    }
  }

  updateItemBy(key: string, fn: (item: any) => any, propagate = true) {
    if (this.state[key]) {
      try {
        this.propagate = propagate;
        const item = this.state[key].value;
        const updatedItem = fn(item);
        this.state[key].next(updatedItem);
      } finally {
        this.propagate = true;
      }
    }
  }

  getItem<T=any>(key: string): T | undefined {
    if (this.state[key])
      return this.state[key].value;
  }

  itemChanges$<T=any>(key: string): Observable<T> | undefined {
    if (this.state[key])
      return this.state[key].pipe(filter(() => this.propagate));
  }

  static mergeStores(stores: CompositionStateStore[]) {
    const mergedStore = new CompositionStateStore();
    const allKeys = new Set<string>();
    for (const store of stores) {
      for (const key of Object.keys(store))
        allKeys.add(key);
    }
    for (const key of allKeys) {
      const subjects = stores.map((s) => s.state[key]).filter((x) => x);
      mergedStore.setItem(key, merge(subjects));
    }
    return mergedStore;
  }

  narrowStore(keys: string[]) {
    const subStore = new CompositionStateStore();
    for (const key of keys) {
      if (this.state[key])
        subStore.state[key] = this.state[key];
    }
    return subStore;
  }
}
