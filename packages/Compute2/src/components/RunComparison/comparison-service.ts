// The API of a mounted Run Comparison view for code outside it; the view gets it through provide/inject.

import type * as DG from 'datagrok-api/dg';
import type {InjectionKey} from 'vue';
import type {ComparisonEntry} from './types';

/** Runs to load and entries to add; the model joins the model list, and is selected on request. */
export interface ComparisonAddition {
  runIds?: string[];
  entries?: ComparisonEntry[];
  model?: DG.Func;
  selectModel?: boolean;
}

export interface ComparisonHandlers {
  add(addition: ComparisonAddition): Promise<void>;
  entryNames(): string[];
}

export class ComparisonService {
  private handlers?: ComparisonHandlers;

  /** Called by the view when it mounts. */
  connect(handlers: ComparisonHandlers) {
    this.handlers = handlers;
  }

  add(addition: ComparisonAddition) {
    return this.handlers!.add(addition);
  }

  get entryNames() {
    return this.handlers!.entryNames();
  }
}

export const COMPARISON_SERVICE: InjectionKey<ComparisonService> = Symbol('ComparisonService');
