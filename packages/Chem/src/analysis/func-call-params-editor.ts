import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Observable, Subject} from 'rxjs';

/** The shape shared by the `@datagrok-libraries/ml` dialog editors (DimReductionBaseEditor,
 * ActivityCliffsEditor, ...): a plain widget builder exposing `getParams()` plus optional
 * history serialization. */
export interface ParamsEditorLike<TParams> {
  getEditor(): HTMLElement;
  getParams(): TParams;
  getStringInput?(): string;
  applyStringInput?(input: string): void | Promise<void>;
}

export interface FuncCallParamsEditorOptions<TParams> {
  /** The wrapped params-producing editor. */
  inner: ParamsEditorLike<TParams>;
  /** Maps `inner.getParams()` to a `{funcInputName: value}` record written into the funccall. */
  map: (params: TParams) => Record<string, unknown>;
  /** Returns whether the current selection is sufficient to enable OK / run the funccall. */
  isValid: (params: TParams) => boolean;
  /** Stable top-level inputs whose `onChanged` should trigger a re-sync. The library editors
   * regenerate some inputs on the fly; the DOM-level listener below is the catch-all for those. */
  stableInputs?: (DG.InputBase | null | undefined)[];
  /** Extra package-specific inputs appended below the wrapped editor (e.g. Chem's `clusterMCS`). */
  extraInputs?: {name: string, input: DG.InputBase}[];
  /** Optional `propertyName -> InputBase` mapping backing `FuncCallEditor.inputFor`. */
  inputFor?: Record<string, DG.InputBase>;
}

/** Adapts a `getParams()`-style dialog editor from `@datagrok-libraries/ml` to the canonical
 * `DG.FuncCallEditor` contract: it renders the wrapped editor, mirrors its current parameters into
 * the live funccall's inputs on every change, and lets the platform own the dialog + execution.
 *
 * Change detection is layered because the wrapped editors regenerate some of their inputs and their
 * gear-icon settings dialogs mutate options without firing any input event:
 *  1. subscriptions on the stable top-level inputs (primary signal),
 *  2. bubbling `input`/`change`/`click` DOM listeners on the editor root (catches regenerated
 *     inputs and settings gears),
 *  3. a re-sync inside `isValid`, which the platform evaluates right before enabling/running OK. */
export class FuncCallParamsEditor<TParams> extends DG.FuncCallEditor {
  private inputChangedSubject: Subject<any> = new Subject<any>();

  constructor(private funcCall: DG.FuncCall, private opts: FuncCallParamsEditorOptions<TParams>) {
    const root = ui.divV([]);
    super(root);
    const editorEl = opts.inner.getEditor();
    for (const extra of opts.extraInputs ?? [])
      editorEl.appendChild(extra.input.root);
    root.appendChild(editorEl);

    for (const input of opts.stableInputs ?? [])
      input?.onChanged.subscribe(() => this.syncCall());
    for (const extra of opts.extraInputs ?? [])
      extra.input.onChanged.subscribe(() => this.syncCall());
    const domSync = () => this.syncCall();
    editorEl.addEventListener('input', domSync, true);
    editorEl.addEventListener('change', domSync, true);
    // click covers settings-gear toggles and checkboxes whose value lands after the event
    editorEl.addEventListener('click', () => setTimeout(domSync), true);

    this.syncCall();
  }

  /** Writes the wrapped editor's current parameters into the live funccall and notifies the platform. */
  private syncCall(): void {
    const mapped = this.opts.map(this.opts.inner.getParams());
    for (const [key, value] of Object.entries(mapped))
      this.funcCall.inputs[key] = value;
    for (const extra of this.opts.extraInputs ?? [])
      this.funcCall.inputs[extra.name] = extra.input.value;
    this.inputChangedSubject.next(null);
  }

  get isValid(): boolean {
    // Pure read only: the platform re-evaluates isValid on every onInputChanged emission, so
    // emitting from here (via syncCall) would spin an infinite loop. Inputs are kept current by the
    // stable-input subscriptions and the bubbling DOM listeners, which fire before the platform
    // reads this. Reading getParams() here is side-effect free.
    return this.opts.isValid(this.opts.inner.getParams());
  }

  inputFor(propertyName: string): DG.InputBase {
    const input = this.opts.inputFor?.[propertyName];
    if (!input)
      throw new Error(`Unknown property name: ${propertyName}`);
    return input;
  }

  getHistoryString(): string {
    return this.opts.inner.getStringInput?.() ?? '';
  }

  loadHistoryString(historyString: string): void {
    if (!historyString || !this.opts.inner.applyStringInput)
      return;
    void Promise.resolve(this.opts.inner.applyStringInput(historyString)).then(() => this.syncCall());
  }

  get onInputChanged(): Observable<any> {
    return this.inputChangedSubject;
  }
}

/** A minimal invalid editor used when preconditions fail (e.g. no molecule column present):
 * shows a message and keeps OK disabled. */
export class MessageFuncCallEditor extends DG.FuncCallEditor {
  private inputChangedSubject: Subject<any> = new Subject<any>();
  constructor(message: string) {
    super(ui.divV([ui.divText(message)]));
  }
  get isValid(): boolean {
    return false;
  }
  get onInputChanged(): Observable<any> {
    return this.inputChangedSubject;
  }
}
