import {EditorItem} from './types';

const HISTORY_LENGTH = 12;

/** Recently created items of one kind, kept in the local storage between the dialogs. */
export class ItemHistory<T extends EditorItem> {
  /** Items saved by the previous dialogs */
  public readonly items: T[];

  // Items created in the currently open dialog
  private justCreatedItems: T[] = [];

  constructor(private storageKey: string) {
    this.items = localStorage[storageKey] ? JSON.parse(localStorage[storageKey]) : [];
  }

  public add(item: T): void {
    this.justCreatedItems.unshift(item);
  }

  public save(): void {
    const compareItems = (a: T, b: T) => JSON.stringify(a) === JSON.stringify(b);
    /** Remove duplicates from just created items (object comparison via JSON.stringify) */
    const justCreatedItems = this.justCreatedItems.filter((val, ind, arr) =>
      arr.findIndex((t) => compareItems(t, val)) === ind);

    /** Remove identical older items from history */
    const olderItems = this.items.filter((old) => !justCreatedItems.find((val) => compareItems(val, old)));

    const newHistoryItems = justCreatedItems.concat(olderItems);
    newHistoryItems.splice(HISTORY_LENGTH);

    localStorage[this.storageKey] = JSON.stringify(newHistoryItems);
  }
}
