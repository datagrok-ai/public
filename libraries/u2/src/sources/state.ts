/* The simplest v1 source (DD14): a value the form itself owns — a draft several controls share,
   a flag other nodes read. One writable signal, one `clear` function, no platform behind it. */
import {signal, Signal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {checkProp, specWarn} from '../spec/spec.js';
import type {BindProp, BindSource} from '../spec/bind-source.js';

export interface StateOptions {
  /** A platform TYPE string; what {@link StateSource.initial} is checked against. */
  type?: string;
  initial?: unknown;
}

export class StateSource extends Component implements BindSource {
  readonly value: Signal<unknown>;
  readonly type: string;
  readonly initial: unknown;

  constructor(options: StateOptions = {}) {
    super();
    this.type = options.type ?? 'string';
    const refusal = options.initial === undefined ? null :
      checkProp({name: 'initial', type: this.type}, options.initial);
    if (refusal !== null)
      specWarn(`u2-state: ${refusal} — starting empty`);
    this.initial = refusal === null && options.initial !== undefined ?
      options.initial : StateSource._empty(this.type);
    this.value = signal(this.initial);
    this.registerFunction({name: 'clear', description: 'Back to the initial value', inputs: [],
      apply: () => this.value.value = this.initial});
  }

  bindStep(name: string): Signal<unknown> | null {
    return name === 'value' || name === '' ? this.value : null;
  }

  bindProps(): BindProp[] {
    return [{name: 'value', type: this.type, writable: true, default: true}];
  }

  private static _empty(type: string): unknown {
    switch (type) {
      case 'int':
      case 'double':
        return 0;
      case 'bool':
        return false;
      case 'string_list':
        return [];
      case 'object':
        return null;
      default:
        return '';
    }
  }
}
