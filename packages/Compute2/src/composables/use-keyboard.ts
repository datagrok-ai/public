import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {useEventListener} from '@vueuse/core';

/** Compute2's styles for platform controls apply under this class, see tailwind.css. */
export const SCOPE_CLASS = 'compute2-scope';

/** Marks the panels Tab goes round in a view: the form and the workflow. */
export const TAB_CYCLE_CLASS = 'compute2-tab-cycle';

const CANDIDATES = 'a[href], button, input, select, textarea, [tabindex]';

const tabStops = (containers: Element[]) => containers.flatMap((container) =>
  Array.from(container.querySelectorAll<HTMLElement>(CANDIDATES)).filter((el) => el.tabIndex >= 0 &&
    !el.matches(':disabled') && el.getClientRects().length > 0 && getComputedStyle(el).visibility !== 'hidden'));

/**
 * Tab and Shift+Tab go round the stops of the containers and never leave them. The focus is moved here rather than by
 * the browser, whose order jumps to the page start past an element with a positive tabindex, such as a grid's canvas.
 */
function tabWithin(e: KeyboardEvent, containers: Element[]) {
  if (e.key !== 'Tab' || e.defaultPrevented || e.ctrlKey || e.altKey || e.metaKey)
    return;
  const stops = tabStops(containers);
  if (!stops.length)
    return;
  e.preventDefault();
  const active = document.activeElement!;
  const index = stops.indexOf(active as HTMLElement);
  const step = e.shiftKey ? -1 : 1;
  if (index >= 0) {
    stops[(index + step + stops.length) % stops.length].focus();
    return;
  }
  // the focus is on no stop, like the dialog itself or an output grid: the nearest stop after it, or before it, else
  // round the other end
  const side = e.shiftKey ? Node.DOCUMENT_POSITION_PRECEDING : Node.DOCUMENT_POSITION_FOLLOWING;
  const near = stops.filter((el) => active.compareDocumentPosition(el) & side);
  (e.shiftKey ? near.pop() ?? stops[stops.length - 1] : near[0] ?? stops[0]).focus();
}

/**
 * While the view is shown, its keys work with the focus anywhere in the view or its ribbon, or lost to the page, which
 * a key on the view's root would miss. Tab goes round the form and workflow panels; from elsewhere in the view, such as
 * an output grid, it goes into them, and outputs, the run history and the ribbon are left to the pointer. The other
 * keys go to onKey. Ctrl/⌘+S never opens the browser's save, wherever the focus is.
 */
export function useViewKeys(view: DG.ViewBase, onKey: (e: KeyboardEvent) => void) {
  useEventListener(document, 'keydown', (e: KeyboardEvent) => {
    if ((e.key !== 'Tab' && !e.ctrlKey && !e.metaKey) || e.defaultPrevented || grok.shell.v?.root !== view.root)
      return;
    const active = document.activeElement;
    if (active === document.body || view.root.contains(active) ||
      view.getRibbonPanels().flat().some((el) => el.contains(active))) {
      if (e.key === 'Tab')
        tabWithin(e, Array.from(view.root.querySelectorAll(`.${TAB_CYCLE_CLASS}`)));
      else
        onKey(e);
    }
    if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 's')
      e.preventDefault();
  });
}

/**
 * Tab does nothing in the dialog: Enter confirms and Esc cancels wherever the focus is, so focusing a button misleads.
 */
export function ignoreTab(dialog: DG.Dialog) {
  dialog.root.addEventListener('keydown', (e) => {
    if (e.key === 'Tab')
      e.preventDefault();
  });
  return dialog;
}

const COMMON_SHORTCUTS = {
  'Enter, Space': 'Press the focused button or icon',
  'Enter, Esc': 'Dialogs: confirm, cancel',
};

const WORKFLOW_SHORTCUTS = {
  'Ctrl/⌘+Enter': 'Run the step focused in the tree, or the open step',
  'Ctrl/⌘+Shift+Enter': 'Update the focused or open step, or the substeps of a focused workflow',
  'Ctrl/⌘+S': 'Save the workflow',
  'Tab, Shift+Tab': 'Move between the steps tree and the open step',
  ...COMMON_SHORTCUTS,
  '↑, ↓, Home, End': 'Tree: go to another step',
  '←, →': 'Tree: collapse or expand a workflow, go to its parent or first step',
  'Enter': 'Tree: open the step',
  'Alt+↑, Alt+↓': 'Tree: move the step in its list',
  'Delete': 'Tree: remove the step',
};

const RFV_SHORTCUTS = {
  'Ctrl/⌘+Enter': 'Run',
  'Ctrl/⌘+S': 'Save the run',
  'Tab, Shift+Tab': 'Move between the inputs and buttons of the form',
  ...COMMON_SHORTCUTS,
};

export const shortcutsRibbonItem = (workflow: boolean) => ({
  icon: 'keyboard',
  tooltip: 'Keyboard shortcuts',
  onClick: () => ignoreTab(ui.dialog('Keyboard shortcuts')
    .add(ui.tableFromMap(workflow ? WORKFLOW_SHORTCUTS : RFV_SHORTCUTS)))
    .show({center: true, modal: true}),
});
