import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {signal, computed, untracked, Control, Scope, Splitter, VirtualTree, divV} from '@datagrok-libraries/u2';
import type {Signal, ReadonlySignal} from '@datagrok-libraries/u2';
import {controlAt} from '@datagrok-libraries/u2/dg';
import {DEMO_TREE, DemoContext, DemoLeaf, leafById, pageShown} from './nav';

const LAST_TAB = 'u2demo.tab';
/** The rail: the narrowest nav whose longest label ("Sections & wizard") still fits. */
const NAV_MIN = 148;
/** The nav is a width in pixels, not a share: a share grows with the window and leaves the
 * capped pages an empty strip on a wide screen. */
const NAV_WIDTH = 180;
const NAV_MAX = 0.4;
/** What the widest sub-demo form needs; the nav gives way before the content goes under it. */
const CONTENT_MIN = 480;
/** The pane width at which the nav folds away unless asked for (css/u2demo.css). */
const NAV_FOLD = 460;

export interface DemoOptions {
  /** The leaf to open on; falls back to `localStorage['u2demo.tab']`, then the first leaf. */
  initial?: string;
}

export interface DemoShell {
  readonly current: ReadonlySignal<DemoLeaf>;
  readonly inspect: Signal<boolean>;
  /** The pane is too narrow for the nav beside the demo, so the nav is away unless `navShown`. */
  readonly navFolded: ReadonlySignal<boolean>;
  readonly navShown: Signal<boolean>;
  navigate(id: string): void;
  rebuild(): void;
}

export interface DemoParts {
  content: Control;
  shell: DemoShell;
  status: ReadonlySignal<string>;
}

export function buildDemo(options?: DemoOptions): DemoParts {
  let shell!: DemoShell;
  let status!: ReadonlySignal<string>;
  const content = Control.build(() => {
    const root = Scope.ambient!;
    const last = signal('');
    // `navigate` is declared below and only ever called from a click, long after this runs
    const ctx: DemoContext = {status: last, navigate: (id) => navigate(id)};
    const inspect = signal(false);
    const navFolded = signal(false);
    const navShown = signal(false);

    const initial = leafById(options?.initial ?? '') ??
      leafById(localStorage.getItem(LAST_TAB) ?? '') ??
      DEMO_TREE[0].children[0];
    const current = signal(initial);

    const host = document.createElement('div');
    host.className = 'u2demo-content';

    let pageScope: Scope | undefined;
    root.own(() => pageScope?.dispose());

    const show = (leaf: DemoLeaf, force = false): void => {
      if (!force && current.peek() === leaf)
        return;
      pageScope?.dispose();
      pageScope = new Scope();
      last.value = '';
      ui.empty(host);
      host.replaceChildren(pageScope.build(() => leaf.build(ctx)));
      current.value = leaf;
      localStorage.setItem(LAST_TAB, leaf.id);
      pageShown(leaf);
    };

    const tree = new VirtualTree<DemoLeaf>();
    tree.root.classList.add('u2demo-nav');
    tree.root.setAttribute('aria-label', 'Demos');
    tree.setRoots(DEMO_TREE.map((g) => ({id: g.id, label: g.label,
      children: g.children.map((leaf) =>
        ({id: leaf.id, label: leaf.label, tooltip: leaf.description, data: leaf}))})));
    tree.expanded.value = new Set(DEMO_TREE.map((g) => g.id));

    root.effect(() => {
      const node = tree.selectedNode.value;
      if (node?.data)
        untracked(() => show(node.data!));
    });

    const navigate = (id: string): void => {
      const leaf = leafById(id);
      if (leaf)
        void tree.expandPath([leaf.group.id, leaf.id]);
    };

    // CAPTURE phase, writing WITHOUT freezing: `grok.shell.o` mutes the channel for a second and
    // drops muted writes, so a freezing write here would swallow the later bubble-phase write a
    // chip (entities/spaces pages) or the convergence '…' button makes for its own object.
    let panelOpened = false;
    host.addEventListener('click', (e) => {
      if (!inspect.peek())
        return;
      // bounded to the page: an unbounded walk resolves whitespace to the shell's own splitter.
      // A click that lands on a heading, a hint or a plain button hits no u2 control — leave the
      // panel showing whatever it holds rather than wiping it with an empty state
      const control = controlAt(e.target as Element, host);
      if (!control)
        return;
      if (!panelOpened) {
        // docking the panel re-parents the view, which blurs the field this click just focused:
        // hand the focus back, or the first keystrokes after an Inspect click go nowhere
        const focused = document.activeElement;
        grok.shell.windows.showContextPanel = true;
        panelOpened = true;                      // open once per view, then respect the user
        if (focused instanceof HTMLElement && focused !== document.activeElement && host.contains(focused))
          focused.focus();
      }
      grok.shell.setCurrentObject(control, false, false);
    }, true);

    // paging is how you read a long sub-demo, not how you leave it: the tree's own list handler
    // moves the selection a page at a time, which navigates. Arrows and Home/End still walk the
    // tree — only Page Up/Down are taken over, and they scroll the page you are reading.
    tree.root.addEventListener('keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key !== 'PageUp' && key !== 'PageDown')
        return;
      e.stopPropagation();
      e.preventDefault();
      host.scrollBy({top: (key === 'PageDown' ? 1 : -1) * Math.max(1, host.clientHeight - 40)});
    }, true);

    // a group row toggles its branch and leaves the selection alone: taking it would move the
    // highlight off the leaf the content still shows, so the tree would lie about where you are
    tree.root.addEventListener('click', (e) => {
      const node = tree.nodeForRow(e.target);
      if (!node || node.data)
        return;
      e.stopPropagation();
      const expanded = new Set(tree.expanded.peek());
      if (!expanded.delete(node.id))
        expanded.add(node.id);
      tree.expanded.value = expanded;
    }, true);

    // the two context-panel behaviours share one panel: turning Inspect off hands it back to the
    // source of the sub-demo that is open, so the toggle round-trips
    root.effect(() => {
      if (!inspect.value)
        untracked(() => pageShown(current.peek()));
    });

    // `show` clears the status tail; a page with no visible state gives no other sign it was rebuilt
    const rebuild = (): void => {
      show(current.peek(), true);
      last.value = 'rebuilt';
    };
    shell = {current, inspect, navFolded, navShown, navigate, rebuild};
    status = computed(() => {
      const leaf = current.value;
      const tail = last.value;
      return `${leaf.group.label} / ${leaf.label}${tail ? ` · ${tail}` : ''}`;
    });

    show(initial, true);
    navigate(initial.id);

    const split = new Splitter([tree.root, host],
      {direction: 'horizontal', sizes: [22, 78], minSize: NAV_MIN});
    split.root.classList.add('u2demo-shell');

    // the splitter holds one minimum for both panels and sizes them in shares; the nav keeps its
    // width in pixels across resizes instead, stops at NAV_MAX, and gives way whenever the content
    // would go under CONTENT_MIN — on a drag and on a resize alike
    let navPx = NAV_WIDTH;
    const fitNav = (): void => {
      const width = split.root.clientWidth;
      if (width <= 0)
        return;
      navFolded.value = width <= NAV_FOLD;
      const nav = Math.min(navPx, Math.max(NAV_MIN, Math.min(NAV_MAX * width, width - CONTENT_MIN))) / width;
      if (Math.abs(split.sizes.peek()[0] - nav) * width > 0.5)
        split.sizes.value = [nav, 1 - nav];
    };
    // only a drag or a sash key is the user's width: the splitter's own clamp on a resize is not
    const bySash = (): boolean => split.root.classList.contains('u2-splitter-dragging') ||
      document.activeElement?.parentElement === split.root;
    root.effect(() => {
      const nav = split.sizes.value[0];
      untracked(() => {
        if (bySash())
          navPx = nav * split.root.clientWidth;
        fitNav();
      });
    });
    root.listen(split.root, 'dblclick', (e) => {
      if ((e.target as Element).parentElement !== split.root)
        return;
      navPx = NAV_WIDTH;
      fitNav();
    });
    const resize = new ResizeObserver(fitNav);
    resize.observe(split.root);
    root.own(() => resize.disconnect());
    const shellRoot = divV([split], 'u2demo-root');
    root.effect(() => shellRoot.classList.toggle('u2demo-nav-shown', navShown.value));
    return shellRoot;
  });
  return {content, shell, status};
}
