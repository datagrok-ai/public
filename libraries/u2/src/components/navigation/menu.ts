/* Popup / context menu. Hand-written introspectable machine per D6: 'closed' | 'open' |
   'submenu-open' exposed as a signal, WAI-ARIA 1.2 menu pattern, focus parked on the panel with
   aria-activedescendant. The fluent surface mirrors DG.Menu so platform users read it unaided.
   Anchored and at a point alike, the menu is a popup of the one overlay stack: Escape and an
   outside press close it there, and the menu installs no page listener of its own. */

import {autoUpdate, computePosition, flip, offset, shift} from '@floating-ui/dom';
import {signal, batch, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {icon as iconElement} from '../display/icon.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface MenuItemOptions {
  icon?: string;
  shortcut?: string;
  enabled?: boolean;
  check?: boolean;
  /** Grey text after the label — what the item stands for ("hydroxyl · O[*:2]"); the details of one
   * menu start in one column. */
  detail?: string;
}

export type MenuState = 'closed' | 'open' | 'submenu-open';

interface MenuNode {
  kind: 'item' | 'separator' | 'group';
  label: string;
  onClick?: () => void;
  options?: MenuItemOptions;
  children?: MenuNode[];
}

interface MenuLevel {
  depth: number;
  panel: HTMLElement;
  entries: {node: MenuNode, el: HTMLElement, path: string}[];
  active: number;
  openIndex: number;
  scope: Scope;
}

const SUBMENU_DELAY = 200;

/** One open level of a menu as a control: what it shows, for tests and tools. It belongs to the
 * menu and closes with it, never to whatever scope is ambient when the menu opens. */
class MenuPanel extends Control {
  constructor(private readonly _level: MenuLevel) {
    super(_level.panel, _level.scope);
  }

  protected get ambientAdoption(): boolean {
    return false;
  }

  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), open: true};
  }

  /** The items, keyed by their path (`File > Open`); separators are not items. */
  getItems(): IItemsState {
    const level = this._level;
    const list = level.entries.filter((e) => e.node.kind !== 'separator').map((e, index) => ({
      key: e.path, label: e.node.label, index, disabled: e.node.options?.enabled === false,
      ...e.node.options?.check === undefined ? {} : {checked: e.node.options.check},
      ...e.node.kind === 'group' ? {expanded: level.entries[level.openIndex] === e} : {},
      rendered: true}));
    return {count: list.length, list, truncated: false};
  }
}

export class Menu {
  private static _seq = 0;

  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<MenuState>;

  private readonly _id = `u2-menu-${++Menu._seq}`;
  private readonly _nodes: MenuNode[] = [];
  private readonly _stack: MenuNode[][] = [];
  private readonly _levels: MenuLevel[] = [];
  private readonly _open = signal(false);
  private readonly _machine = signal<MenuState>('closed');
  private _scope: Scope | undefined;
  private _returnFocus: HTMLElement | undefined;
  private _hoverCancel: (() => void) | undefined;

  constructor() {
    this.isOpen = this._open;
    this.machineState = this._machine;
    this._stack.push(this._nodes);
  }

  static popup(): Menu {
    return new Menu();
  }

  /** `onClick` runs owned by whoever added the item (none once that owner is gone), after the menu closed. */
  item(label: string, onClick: () => void, options?: MenuItemOptions): Menu {
    this._current.push({kind: 'item', label, onClick: Menu._owned(onClick), options});
    return this;
  }

  separator(): Menu {
    this._current.push({kind: 'separator', label: ''});
    return this;
  }

  /** Starts a submenu; subsequent items go into it until {@link endGroup}. */
  group(label: string): Menu {
    const children: MenuNode[] = [];
    this._current.push({kind: 'group', label, children});
    this._stack.push(children);
    return this;
  }

  endGroup(): Menu {
    if (this._stack.length > 1)
      this._stack.pop();
    return this;
  }

  show(options: {anchor: HTMLElement} | {x: number, y: number}): void {
    this.close();
    const focused = document.activeElement;
    this._returnFocus = focused instanceof HTMLElement ? focused : undefined;
    const scope = new Scope(undefined, {label: 'menu'});
    this._scope = scope;
    const level = this._buildLevel(this._nodes, 0, '');
    this._levels.push(level);
    this._open.value = true;
    this._machine.value = 'open';
    const anchor = 'anchor' in options ? options.anchor : Menu._point(level.panel, options.x, options.y, scope);
    Overlay.show(anchor, level.panel, scope);
    level.scope.listen(level.panel, OVERLAY_CLOSE_EVENT, () => this.close());
    level.panel.focus({preventScroll: true});
  }

  close(): void {
    const scope = this._scope;
    if (!scope)
      return;
    this._scope = undefined;
    this._hoverCancel = undefined;
    while (this._levels.length > 0)
      Menu._disposeLevel(this._levels.pop()!);
    scope.dispose();
    this._open.value = false;
    this._machine.value = 'closed';
    const restore = this._returnFocus;
    this._returnFocus = undefined;
    if (restore && restore.isConnected)
      restore.focus();
  }

  /** Right-click on `el` opens the menu at the cursor; the listener is released with `scope`. */
  bindContext(el: HTMLElement, scope: Scope): void {
    scope.listen(el, 'contextmenu', (e) => {
      e.preventDefault();
      const me = e as MouseEvent;
      this.show({x: me.clientX, y: me.clientY});
    });
  }

  private get _current(): MenuNode[] {
    return this._stack[this._stack.length - 1];
  }

  /** A zero-size anchor at the point, moved left of it when the panel would run past the right
   * edge — the overlay places the panel under its anchor and flips it up when there is no room. */
  private static _point(panel: HTMLElement, x: number, y: number, scope: Scope): HTMLElement {
    panel.style.position = 'fixed';
    panel.style.visibility = 'hidden';
    Overlay.host.append(panel);
    const width = panel.offsetWidth;
    panel.style.visibility = '';
    const anchor = document.createElement('div');
    anchor.className = 'u2-menu-point';
    const left = x + width > window.innerWidth - 4 ? Math.max(4, x - width) : x;
    anchor.style.cssText = `position: fixed; left: ${left}px; top: ${y}px; width: 0; height: 0;`;
    Overlay.host.append(anchor);
    scope.own(() => anchor.remove());
    return anchor;
  }

  private _buildLevel(nodes: MenuNode[], depth: number, path: string): MenuLevel {
    const panel = document.createElement('div');
    panel.className = 'u2-menu';
    panel.dataset.u2 = 'menu';
    panel.tabIndex = -1;
    panel.setAttribute('role', 'menu');
    const level: MenuLevel = {depth, panel, entries: [], active: -1, openIndex: -1,
      scope: this._scope!.child('menu-level')};
    const detailed = nodes.some((n) => n.options?.detail !== undefined);
    panel.classList.toggle('u2-menu-detailed', detailed);
    for (const node of nodes) {
      const itemPath = path === '' ? node.label : `${path} > ${node.label}`;
      const el = this._renderNode(node, depth, level.entries.length, itemPath, detailed);
      level.entries.push({node, el, path: itemPath});
      panel.append(el);
    }
    new MenuPanel(level);
    const scope = level.scope;
    scope.listen(panel, 'keydown', (e) => this._onKeyDown(level, e as KeyboardEvent));
    scope.listen(panel, 'pointerover', (e) => this._onPointerOver(level, e));
    scope.listen(panel, 'pointerdown', (e) => e.preventDefault());
    scope.listen(panel, 'click', (e) => this._onClick(level, e));
    return level;
  }

  private _renderNode(node: MenuNode, depth: number, index: number, path: string, detailed: boolean): HTMLElement {
    if (node.kind === 'separator') {
      const sep = document.createElement('div');
      sep.className = 'u2-menu-separator';
      sep.setAttribute('role', 'separator');
      return sep;
    }
    const options = node.options ?? {};
    const isGroup = node.kind === 'group';
    const el = document.createElement('div');
    el.className = isGroup ? 'u2-menu-item u2-menu-item-sub' : 'u2-menu-item';
    el.id = `${this._id}-${depth}-${index}`;
    el.dataset.index = String(index);
    el.dataset.u2Item = path;
    el.setAttribute('role', options.check === undefined ? 'menuitem' : 'menuitemcheckbox');
    el.setAttribute('aria-disabled', String(options.enabled === false));
    if (options.check !== undefined)
      el.setAttribute('aria-checked', String(options.check));
    if (isGroup) {
      el.setAttribute('aria-haspopup', 'menu');
      el.setAttribute('aria-expanded', 'false');
    }
    const icon = document.createElement('span');
    icon.className = 'u2-menu-icon';
    if (options.check)
      icon.textContent = '✓';
    else if (options.icon)
      icon.append(iconElement(options.icon));
    const label = document.createElement('span');
    label.className = 'u2-menu-label';
    label.dataset.u2Part = 'item-label';
    label.textContent = node.label;
    const accel = document.createElement('span');
    accel.className = isGroup ? 'u2-menu-arrow' : 'u2-menu-shortcut';
    accel.textContent = isGroup ? '' : (options.shortcut ?? '');
    el.append(icon, label);
    // an item of a detailed menu with no detail of its own ("Edit definition…") lets its label run on
    // across the details' column, so a long label does not push every detail right
    label.classList.toggle('u2-menu-label-wide', detailed && options.detail === undefined);
    if (options.detail !== undefined) {
      const detail = document.createElement('span');
      detail.className = 'u2-menu-detail';
      detail.dataset.u2Part = 'item-detail';
      detail.textContent = options.detail;
      el.append(detail);
    }
    el.append(accel);
    return el;
  }

  private _openSubmenu(level: MenuLevel, index: number, focus: boolean): void {
    const entry = level.entries[index];
    if (!entry || entry.node.kind !== 'group')
      return;
    if (level.openIndex !== index) {
      this._closeFrom(level.depth + 1);
      const sub = this._buildLevel(entry.node.children!, level.depth + 1, entry.path);
      this._levels.push(sub);
      level.openIndex = index;
      entry.el.setAttribute('aria-expanded', 'true');
      entry.el.classList.add('u2-menu-item-open');
      sub.panel.style.position = 'fixed';
      sub.panel.style.left = '0';
      sub.panel.style.top = '0';
      // Nested in the owning item rather than in Overlay.host: the root overlay dismisses on any
      // pointerdown outside its content, and a detached submenu panel would count as outside.
      entry.el.append(sub.panel);
      sub.scope.own(autoUpdate(entry.el, sub.panel, () => {
        computePosition(entry.el, sub.panel, {
          strategy: 'fixed',
          placement: 'right-start',
          middleware: [
            offset({mainAxis: -1, crossAxis: -3}),
            flip({fallbackPlacements: ['left-start']}),
            shift({padding: 4}),
          ],
        }).then(({x, y}) => {
          sub.panel.style.left = `${x}px`;
          sub.panel.style.top = `${y}px`;
        });
      }));
    }
    const child = this._levels[level.depth + 1];
    this._machine.value = 'submenu-open';
    if (focus) {
      child.panel.focus({preventScroll: true});
      if (child.active < 0)
        Menu._setActive(child, Menu._step(child, 1));
    }
  }

  private _closeFrom(depth: number): void {
    while (this._levels.length > depth)
      Menu._disposeLevel(this._levels.pop()!);
    const parent = this._levels[this._levels.length - 1];
    if (parent && parent.openIndex >= 0) {
      const el = parent.entries[parent.openIndex].el;
      el.setAttribute('aria-expanded', 'false');
      el.classList.remove('u2-menu-item-open');
      parent.openIndex = -1;
    }
    this._machine.value = this._levels.length > 1 ? 'submenu-open' : 'open';
  }

  private _activate(level: MenuLevel, index: number): void {
    const entry = level.entries[index];
    if (!entry || !Menu._navigable(entry.node))
      return;
    if (entry.node.kind === 'group') {
      this._openSubmenu(level, index, true);
      return;
    }
    const onClick = entry.node.onClick;
    this.close();
    if (onClick)
      onClick();
  }

  private _onClick(level: MenuLevel, e: Event): void {
    const index = Menu._indexAt(level, e);
    if (index >= 0)
      this._activate(level, index);
  }

  private _onPointerOver(level: MenuLevel, e: Event): void {
    const index = Menu._indexAt(level, e);
    if (index < 0)
      return;
    // Deeper levels handle their own hover; letting it bubble would reset this level's intent timer.
    e.stopPropagation();
    const entry = level.entries[index];
    const isGroup = entry.node.kind === 'group';
    if (Menu._navigable(entry.node))
      Menu._setActive(level, index);
    this._hoverCancel?.();
    this._hoverCancel = undefined;
    if (isGroup ? level.openIndex === index : level.openIndex < 0)
      return;
    this._hoverCancel = this._scope?.setTimeout(() => {
      this._hoverCancel = undefined;
      if (isGroup)
        this._openSubmenu(level, index, false);
      else
        this._closeFrom(level.depth + 1);
    }, SUBMENU_DELAY);
  }

  private _onKeyDown(level: MenuLevel, e: KeyboardEvent): void {
    if (e.target !== level.panel)
      return;
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        Menu._setActive(level, Menu._step(level, 1));
        break;
      case 'ArrowUp':
        e.preventDefault();
        Menu._setActive(level, Menu._step(level, -1));
        break;
      case 'Home':
        e.preventDefault();
        Menu._setActive(level, Menu._edge(level, 1));
        break;
      case 'End':
        e.preventDefault();
        Menu._setActive(level, Menu._edge(level, -1));
        break;
      case 'ArrowRight':
        if (level.entries[level.active] && level.entries[level.active].node.kind === 'group') {
          e.preventDefault();
          this._openSubmenu(level, level.active, true);
        }
        break;
      case 'ArrowLeft':
        if (level.depth > 0) {
          e.preventDefault();
          const parent = this._levels[level.depth - 1];
          this._closeFrom(level.depth);
          parent.panel.focus({preventScroll: true});
        }
        break;
      case 'Enter':
        e.preventDefault();
        this._activate(level, level.active);
        break;
      case 'Escape':
        e.preventDefault();
        this.close();
        break;
      default:
        if (e.key.length === 1 && !e.ctrlKey && !e.altKey && !e.metaKey) {
          const index = Menu._jump(level, e.key.toLowerCase());
          if (index >= 0) {
            e.preventDefault();
            Menu._setActive(level, index);
          }
        }
    }
  }

  /** A handler run the way its registrant would run it: with the scope ambient where it was
   * registered (none once that scope is gone), in one batch. */
  private static _owned(fn: () => void): () => void {
    const owner = Scope.ambient;
    return () => batch(() => Scope.runWith(owner?.isDisposed ? undefined : owner, fn));
  }

  private static _disposeLevel(level: MenuLevel): void {
    level.scope.dispose();
    level.panel.remove();
  }

  private static _navigable(node: MenuNode): boolean {
    return node.kind !== 'separator' && (!node.options || node.options.enabled !== false);
  }

  private static _indexAt(level: MenuLevel, e: Event): number {
    const el = (e.target as Element).closest('.u2-menu-item') as HTMLElement | null;
    // Submenu panels live inside their item, so an inner hit resolves to the inner item.
    if (!el || el.parentElement !== level.panel)
      return -1;
    return Number(el.dataset.index);
  }

  private static _setActive(level: MenuLevel, index: number): void {
    level.active = index;
    for (let i = 0; i < level.entries.length; i++)
      level.entries[i].el.classList.toggle('u2-menu-item-active', i === index);
    if (index < 0)
      level.panel.removeAttribute('aria-activedescendant');
    else
      level.panel.setAttribute('aria-activedescendant', level.entries[index].el.id);
  }

  private static _step(level: MenuLevel, delta: number): number {
    const count = level.entries.length;
    if (count === 0)
      return -1;
    const from = level.active >= 0 ? level.active : (delta > 0 ? -1 : 0);
    for (let k = 1; k <= count; k++) {
      const i = ((from + delta * k) % count + count) % count;
      if (Menu._navigable(level.entries[i].node))
        return i;
    }
    return -1;
  }

  private static _edge(level: MenuLevel, delta: number): number {
    const count = level.entries.length;
    for (let i = delta > 0 ? 0 : count - 1; i >= 0 && i < count; i += delta) {
      if (Menu._navigable(level.entries[i].node))
        return i;
    }
    return -1;
  }

  private static _jump(level: MenuLevel, letter: string): number {
    const count = level.entries.length;
    for (let k = 1; k <= count; k++) {
      const i = (level.active + k + count) % count;
      const node = level.entries[i].node;
      if (Menu._navigable(node) && node.label.toLowerCase().startsWith(letter))
        return i;
    }
    return -1;
  }
}
