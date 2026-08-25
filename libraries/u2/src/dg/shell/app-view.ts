/* A u2 component tree as a platform view: the shell owns the chrome (ribbon, toolbox, status bar),
   so a u2 app never rebuilds it. Chrome lives outside the view root, where the Dart kill-walk never
   reaches, which is why components passed as chrome are disposed with the content instead. */
import * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {Signal, ReadonlySignal, rawEffect} from '../../core/signals.js';

type ChromeItem = Control | HTMLElement;

export interface AppViewOptions {
  name: string;
  /** The view content — disposed through the platform kill channel when the view closes. */
  content: Control;
  /** Ribbon panel groups → `view.setRibbonPanels`. Main view controls (filter, refresh, mode
   * switches, menu bar) belong here, not inside the content area. */
  ribbon?: ChromeItem[][];
  /** → the view's toolbox pane (`ViewBase.toolbox`): navigation and structure, not view controls. */
  toolbox?: ChromeItem;
  /** → the shell's per-view status bar. A signal renders as a live text panel; elements and
   * components are placed as given. */
  status?: ReadonlySignal<string> | ChromeItem | ChromeItem[];
}

export function appView(options: AppViewOptions): DG.ViewBase {
  const content = options.content;
  const view = DG.View.fromRoot(content.root);
  view.name = options.name;
  // core contract (2026-08): a docked pane's ✕ kills elements carrying this attribute
  content.root.setAttribute('data-kill-on-close', 'true');
  DG.Widget.registerCleanup(content.root, () => content.dispose());

  const place = (item: ChromeItem): HTMLElement => {
    if (!Control.is(item))
      return item;
    content.own(() => item.dispose());
    return item.root;
  };
  const panel = (item: ChromeItem): HTMLDivElement => {
    const el = place(item);
    if (el instanceof HTMLDivElement)
      return el;
    const wrap = document.createElement('div');
    wrap.append(el);
    return wrap;
  };

  if (options.ribbon)
    view.setRibbonPanels(options.ribbon.map((group) => group.map(place)));
  if (options.toolbox)
    view.toolbox = place(options.toolbox);
  const status = options.status;
  if (status instanceof Signal) {
    const el = document.createElement('div');
    content.own(rawEffect(() => {
      el.textContent = status.value;
    }));
    view.statusBarPanels = [el];
  } else if (status != null) {
    const items = (Array.isArray(status) ? status : [status]) as ChromeItem[];
    view.statusBarPanels = items.map(panel);
  }
  return view;
}
