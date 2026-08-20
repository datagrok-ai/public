/* The status half the AsyncValue-backed sources (func-source, entity-ref) answer alike: the run's
   state projected over the design-time preview, and the `state`/`error` binding steps. Also owns
   the `DesignData` policy type the projection is keyed by (re-exported by func-source.ts) and
   `AsyncSteps.props`, the two status rows both sources' `bindProps()` open with. */
import {computed, Signal, ReadonlySignal} from '../core/signals.js';
import type {AsyncValue, AsyncValueState} from '../core/async-value.js';
import type {BindProp} from '../spec/bind-source.js';

/** What a source does while the designer is up (DD9): run for real, describe itself without
 * running, or preview the sample stored with the spec. */
export type DesignData = 'live' | 'schema' | 'sample';

export class AsyncSteps<T> {
  readonly state: ReadonlySignal<AsyncValueState<T>>;
  private readonly _stateStep: Signal<unknown>;
  private readonly _errorStep: Signal<unknown>;

  /** `preview` is read lazily, under the `sample` policy while nothing has run: the tag-level
   * designPreview is stamped onto `meta` only after the source's constructor returns. */
  constructor(async: AsyncValue<T>, live: boolean, designData: DesignData,
    preview: () => T | undefined) {
    // The engine exists whatever the policy, because Refresh is the one run design time always
    // allows (DD9): a policy that silently swallowed the user's own ask is what made a source with
    // no auto-kick look broken. What a run brought back outranks the preview it replaces.
    this.state = live ? async.state :
      computed((): AsyncValueState<T> => {
        const at = async.state.value;
        if (at.kind !== 'idle' || designData !== 'sample')
          return at;
        const sample = preview();
        return sample === undefined ? at : {kind: 'ready', value: sample};
      });
    this._stateStep = computed(() => this.state.value.kind) as unknown as Signal<unknown>;
    // the failure is bindable, so a form can show WHY it is empty — without it a run that threw
    // and a run that returned nothing look exactly alike anywhere in the designer
    this._errorStep = computed(() => {
      const at = this.state.value;
      return at.kind === 'error' ? at.message : undefined;
    }) as unknown as Signal<unknown>;
  }

  step(name: string): Signal<unknown> | null {
    return name === 'state' ? this._stateStep : name === 'error' ? this._errorStep : null;
  }

  static props(error: string): BindProp[] {
    return [
      {name: 'state', type: 'string', description: 'idle, loading, ready or error'},
      {name: 'error', type: 'string', description: error},
    ];
  }
}
