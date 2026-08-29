import * as grok from 'datagrok-api/grok';
import {signal, computed, untracked, Control, Scope, Splitter, VirtualTree, divV} from '@datagrok-libraries/u2';
import type {Signal, ReadonlySignal} from '@datagrok-libraries/u2';
import {controlAt} from '@datagrok-libraries/u2/src/dg/index.js';
import {DEMO_TREE, DemoContext, DemoLeaf, leafById, pageShown} from './nav';

const LAST_TAB = 'u2demo.tab';

export interface DemoOptions {
  /** The leaf to open on; falls back to `localStorage['u2demo.tab']`, then the first leaf. */
  initial?: string;
}

export interface DemoShell {
  readonly current: ReadonlySignal<DemoLeaf>;
  readonly inspect: Signal<boolean>;
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
    const inspect = signal(true);

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
      host.replaceChildren(Scope.runWith(pageScope, () => leaf.build(ctx)));
      current.value = leaf;
      localStorage.setItem(LAST_TAB, leaf.id);
      pageShown(leaf);
    };

    const tree = new VirtualTree<DemoLeaf>();
    tree.root.classList.add('u2demo-nav');
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
        grok.shell.windows.showContextPanel = true;
        panelOpened = true;                      // open once per view, then respect the user
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
    shell = {current, inspect, navigate, rebuild};
    status = computed(() => {
      const leaf = current.value;
      const tail = last.value;
      return `${leaf.group.label} / ${leaf.label}${tail ? ` · ${tail}` : ''}`;
    });

    show(initial, true);
    navigate(initial.id);

    const split = new Splitter([tree.root, host],
      {direction: 'horizontal', sizes: [22, 78], minSize: 140});
    split.root.classList.add('u2demo-shell');
    return divV([split], 'u2demo-root');
  });
  return {content, shell, status};
}
