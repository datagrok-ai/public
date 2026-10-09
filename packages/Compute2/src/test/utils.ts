import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {after, awaitCheck} from '@datagrok-libraries/test/src/test';
import {loadInstanceState} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/runtime/funccall-utils';
import {
  isFuncCallSerializedState, PipelineSerializedState,
} from '@datagrok-libraries/compute-utils/reactive-tree-driver/src/config/PipelineInstance';

/** Wait for a selector to appear inside a root element, return the matched element. */
export async function awaitElement(
  root: HTMLElement, selector: string, msg: string, timeout = 5000,
): Promise<HTMLElement> {
  await awaitCheck(
    () => root.querySelector(selector) !== null,
    `${msg}: selector "${selector}" not found`, timeout,
  );
  return root.querySelector(selector) as HTMLElement;
}

/** Launch a Compute2 app function by name and poll until the view appears in shell. */
export async function launchApp(name: string, timeout = 10000): Promise<DG.ViewBase> {
  await grok.functions.call(`Compute2:${name}`);
  let view: DG.ViewBase | undefined;
  await awaitCheck(() => {
    const v = grok.shell.v;
    if (v && v.name === name) {
      view = v;
      return true;
    }
    return false;
  }, `View "${name}" did not appear in shell`, timeout);
  return view!;
}

/** Close a view, swallowing errors if already closed. */
export function closeView(view: DG.ViewBase): void {
  try { view.close(); } catch (_) { /* already closed */ }
}

/** Registers views a category's tests open, so every one of them is closed after the category. */
export function closeViewsAfter() {
  const views: DG.ViewBase[] = [];
  after(async () => views.splice(0).forEach(closeView));
  return <T extends DG.ViewBase>(view: T): T => {
    views.push(view);
    return view;
  };
}

const savedCallIds = (state: PipelineSerializedState): string[] => isFuncCallSerializedState(state) ?
  (state.funcCallId ? [state.funcCallId] : []) : state.steps.flatMap(savedCallIds);

/** Deletes a saved workflow run with its step calls; plain calls.delete() silently keeps them. */
export async function deleteSavedWorkflow(id: string) {
  const [state] = await loadInstanceState(id);
  const calls = grok.dapi.functions.calls.allPackageVersions();
  for (const callId of [id, ...savedCallIds(state)]) {
    const call = await calls.find(callId);
    if (call)
      await calls.delete(call);
  }
}

/** Wait for WebComponents custom elements to be registered. */
export async function awaitWebComponents(timeout = 10000): Promise<void> {
  await awaitCheck(
    () => customElements.get('dg-viewer') !== undefined,
    'WebComponents not initialized', timeout,
  );
}

/** Find a button inside a root element by partial text content (case-insensitive). */
export function findButton(root: HTMLElement, text: string): HTMLButtonElement {
  const btn = Array.from(root.querySelectorAll('button')).find(
    (b) => b.textContent?.toLowerCase().includes(text.toLowerCase()),
  );
  if (!btn)
    throw new Error(`Button containing "${text}" not found`);
  return btn as HTMLButtonElement;
}

/**
 * Click a button and poll until a condition becomes true.
 * Use this instead of delay() + awaitCheck() after interactions.
 */
export async function clickAndAwait(
  btn: HTMLElement, check: () => boolean, msg: string, timeout = 5000,
): Promise<void> {
  btn.click();
  await awaitCheck(check, msg, timeout);
}

/**
 * Snapshot the current innerHTML of a root, click, then poll until innerHTML differs.
 * Use when the click replaces content but the same selectors will appear again.
 */
export async function clickAndAwaitChange(
  btn: HTMLElement, root: HTMLElement, msg: string, timeout = 5000,
): Promise<void> {
  const before = root.innerHTML;
  btn.click();
  await awaitCheck(
    () => root.innerHTML !== before,
    msg, timeout,
  );
}
