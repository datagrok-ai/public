/* Popup / context menu. Hand-written introspectable machine per D6: 'closed' | 'open' |
   'submenu-open' exposed as a signal, WAI-ARIA 1.2 menu pattern, focus parked on the panel with
   aria-activedescendant. The fluent surface mirrors DG.Menu so platform users read it unaided. */

import {autoUpdate, computePosition, flip, offset, shift} from '@floating-ui/dom';
import {signal, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {icon as iconElement} from '../display/icon.js';

export interface MenuItemOptions {
  icon?: string;
  shortcut?: string;
  enabled?: boolean;
  check?: boolean;
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
  entries: {node: MenuNode, el: HTMLElement}[];
  active: number;
  openIndex: number;
  disposers: (() => void)[];
}

const SUBMENU_DELAY = 200;

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
  private _hoverTimer = 0;

  constructor() {
    this.isOpen = this._open;
    this.machineState = this._machine;
    this._stack.push(this._nodes);
  }

  static popup(): Menu {
    return new Menu();
  }

  item(label: string, onClick: () => void, options?: MenuItemOptions): Menu {
    this._current.push({kind: 'item', label, onClick, options});
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
    const scope = new Scope();
    this._scope = scope;
    const level = this._buildLevel(this._nodes, 0);
    this._levels.push(level);
    this._open.value = true;
    this._machine.value = 'open';
    if ('anchor' in options) {
      Overlay.show(options.anchor, level.panel, scope);
      this._listen(level, level.panel, OVERLAY_CLOSE_EVENT, () => this.close());
    } else
      this._showAt(level.panel, options.x, options.y, scope);
    level.panel.focus({preventScroll: true});
  }

  close(): void {
    const scope = this._scope;
    if (!scope)
      return;
    this._scope = undefined;
    window.clearTimeout(this._hoverTimer);
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
    const handler = (e: Event) => {
      e.preventDefault();
      const me = e as MouseEvent;
      this.show({x: me.clientX, y: me.clientY});
    };
    el.addEventListener('contextmenu', handler);
    scope.own(() => el.removeEventListener('contextmenu', handler));
  }

  private get _current(): MenuNode[] {
    return this._stack[this._stack.length - 1];
  }

  private _showAt(panel: HTMLElement, x: number, y: number, scope: Scope): void {
    panel.style.position = 'fixed';
    panel.style.left = '0';
    panel.style.top = '0';
    Overlay.host.append(panel);

    const onPointerDown = (e: Event) => {
      if (!panel.contains(e.target as Node))
        this.close();
    };
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Escape')
        this.close();
    };
    document.addEventListener('pointerdown', onPointerDown, true);
    document.addEventListener('keydown', onKeyDown);
    scope.own(() => {
      document.removeEventListener('pointerdown', onPointerDown, true);
      document.removeEventListener('keydown', onKeyDown);
    });

    computePosition({getBoundingClientRect: () => new DOMRect(x, y, 0, 0)}, panel, {
      strategy: 'fixed',
      placement: 'right-start',
      middleware: [flip({fallbackPlacements: ['left-start', 'right-end', 'left-end']}), shift({padding: 4})],
    }).then(({x: left, y: top}) => {
      panel.style.left = `${left}px`;
      panel.style.top = `${top}px`;
    });
  }

  private _buildLevel(nodes: MenuNode[], depth: number): MenuLevel {
    const panel = document.createElement('div');
    panel.className = 'u2-menu';
    panel.dataset.u2 = 'menu';
    panel.tabIndex = -1;
    panel.setAttribute('role', 'menu');
    const level: MenuLevel = {depth, panel, entries: [], active: -1, openIndex: -1, disposers: []};
    for (const node of nodes) {
      const el = this._renderNode(node, depth, level.entries.length);
      level.entries.push({node, el});
      panel.append(el);
    }
    this._listen(level, panel, 'keydown', (e) => this._onKeyDown(level, e as KeyboardEvent));
    this._listen(level, panel, 'pointerover', (e) => this._onPointerOver(level, e));
    this._listen(level, panel, 'pointerdown', (e) => e.preventDefault());
    this._listen(level, panel, 'click', (e) => this._onClick(level, e));
    return level;
  }

  private _renderNode(node: MenuNode, depth: number, index: number): HTMLElement {
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
    label.textContent = node.label;
    const accel = document.createElement('span');
    accel.className = isGroup ? 'u2-menu-arrow' : 'u2-menu-shortcut';
    accel.textContent = isGroup ? '▸' : (options.shortcut ?? '');
    el.append(icon, label, accel);
    return el;
  }

  private _openSubmenu(level: MenuLevel, index: number, focus: boolean): void {
    const entry = level.entries[index];
    if (!entry || entry.node.kind !== 'group')
      return;
    if (level.openIndex !== index) {
      this._closeFrom(level.depth + 1);
      const sub = this._buildLevel(entry.node.children!, level.depth + 1);
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
      sub.disposers.push(autoUpdate(entry.el, sub.panel, () => {
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
      }, {animationFrame: true}));
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
    window.clearTimeout(this._hoverTimer);
    if (isGroup ? level.openIndex === index : level.openIndex < 0)
      return;
    this._hoverTimer = window.setTimeout(() => {
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

  private _listen(level: MenuLevel, target: EventTarget, type: string, handler: (e: Event) => void): void {
    target.addEventListener(type, handler);
    level.disposers.push(() => target.removeEventListener(type, handler));
  }

  private static _disposeLevel(level: MenuLevel): void {
    for (const d of level.disposers)
      d();
    level.disposers.length = 0;
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
