/* Schema-driven query builder over one immutable `FilterGroup`. Simple mode shows the root group
   as a list of condition rows joined by one `and`/`or` toggle; horizontal lays the same rows out
   inline with a connector chip between them; advanced mode renders sub-groups recursively (a
   header per group, a left rule per depth), adds `not`, `+ group` and the drag handles. Rows are
   keyed by node id and reordered in place under whichever group holds them, so an edit never
   detaches the editor being typed into.
   Automation: `data-u2-node="<id>"` on rows and groups; `data-u2-part` values `prop op prop-text
   op-text value value2 problem add add-group remove connector not handle rows mode none
   hint hint-mode flatten query`. */
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Tooltip} from '../../core/tooltip.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {div, span, link, button} from '../../core/elements.js';
import type {IWidgetStatus} from '../../core/widget-like.js';
import {Filters} from '../../core/filter/index.js';
import type {FilterCondition, FilterGroup, FilterNode, FilterProblem, FilterProperty, FilterSchema, FilterTarget,
  FilterTemplate, FilterValueEditorFactory, Lock} from '../../core/filter/index.js';
import {FilterRow, FilterNestedRow, FilterGroupRow} from './filter-row.js';
import type {FilterRowHost} from './filter-row.js';
import {FilterDragLayer} from './filter-dnd.js';
import type {FilterDropHit} from './filter-dnd.js';

export type FilterMode = 'simple' | 'advanced';

export interface FilterBuilderOptions extends InputOptions<FilterGroup> {
  schema: FilterSchema;
  /** 'simple' by default; a Signal is adopted as the control's own. */
  mode?: FilterMode | Signal<FilterMode>;
  /** 'vertical' (a header with the and/or toggle, one row per line) or 'horizontal' (rows
   * inline, wrapping, a connector chip between them). Advanced mode is always vertical. */
  orientation?: 'vertical' | 'horizontal';
  template?: FilterTemplate;
  /** What the tree must be expressible for; validation only. */
  target?: FilterTarget;
  /** Footer line with the canonical query string; off by default. */
  showQuery?: boolean;
  /** Value editors beyond the core kind defaults; {@link FilterBuilder.defaultEditors} otherwise. */
  editors?: FilterValueEditorFactory;
}

export interface FilterBuilderStatus extends IWidgetStatus {
  tree: FilterGroup;
  query: string;
  problems: FilterProblem[];
  mode: FilterMode;
}

type Row = FilterRow | FilterNestedRow | FilterGroupRow;
type RowKind = 'cond' | 'nested' | 'group';

export class FilterBuilder extends Input<FilterGroup, FilterBuilderOptions> {
  /** The platform layer assigns its factory here at import; null = core kind defaults only. */
  static defaultEditors: FilterValueEditorFactory | null = null;

  readonly query: ReadonlySignal<string>;
  readonly problems: ReadonlySignal<FilterProblem[]>;
  readonly mode: Signal<FilterMode>;

  // every field below is built by createEditor(), which the base constructor calls before
  // subclass initializers would run — none of them may carry an inline initializer
  private _mode!: Signal<FilterMode>;
  private _query!: ReadonlySignal<string>;
  private _problems!: ReadonlySignal<FilterProblem[]>;
  private _rows!: Map<string, Row>;
  /** Per node: owns the row and its join chip; released with the node, disowning itself. */
  private _scopes!: Map<string, Scope>;
  private _joins!: Map<string, HTMLElement>;
  private _host!: FilterRowHost;
  private _panel!: HTMLElement;
  private _rootRow!: FilterGroupRow;
  private _modeLink!: HTMLElement;
  private _none!: HTMLElement;
  private _hint!: HTMLElement;
  private _hintText!: HTMLElement;
  private _hintMode!: HTMLElement;
  private _flatten!: HTMLElement;
  private _queryEl: HTMLElement | undefined;
  private _drag!: FilterDragLayer;
  private _horizontal!: boolean;
  private _warned!: boolean;

  constructor(options: FilterBuilderOptions) {
    super(options, Filters.group('and'));
    this.query = this._query;
    this.problems = this._problems;
    this.mode = this._mode;
    this.root.classList.add('u2-filter-builder');
    this.root.dataset.u2 = 'filter-builder';
  }

  /** Appends a condition on the first offered property with its first operator; the new id. */
  addCondition(parentId?: string, index?: number): string {
    const root = this.value.peek();
    const prop = this._host.properties()[0];
    if (!prop)
      return '';
    const cond = Filters.cond(prop.name, this._host.operators(prop)[0]?.id ?? '=');
    this._append(root, cond, parentId, index);
    (this._rows.get(cond.id) as FilterRow | undefined)?.focus();
    return cond.id;
  }

  /** Appends an empty sub-group with the other connector; edited by advanced mode. */
  addGroup(parentId?: string, index?: number): string {
    const root = this.value.peek();
    const parent = (parentId && Filters.find(root, parentId)) || root;
    const group = Filters.group(Filters.isGroup(parent) && parent.op === 'and' ? 'or' : 'and');
    this._append(root, group, parentId, index);
    return group.id;
  }

  removeNode(id: string): void {
    this.value.value = Filters.remove(this.value.peek(), id);
  }

  /** False — with the hint shown — when simple mode cannot show the tree or the template keeps
   * the builder simple. */
  setMode(mode: FilterMode): boolean {
    if (mode === 'simple' && !Filters.isFlat(this.value.peek())) {
      this._showHint(this.value.peek());
      return false;
    }
    if (mode === 'advanced' && this.options.template?.allowAdvanced === false) {
      this._showHint(null);
      return false;
    }
    this._mode.value = mode;
    return true;
  }

  getWidgetStatus(): FilterBuilderStatus {
    const status = super.getWidgetStatus() as FilterBuilderStatus;
    status.parts.rows = this._rootRow.rows;
    status.parts.add = this._rootRow.add;
    status.parts.connector = this._rootRow.connector.root;
    status.parts.mode = this._modeLink;
    if (this._queryEl)
      status.parts.query = this._queryEl;
    status.tree = this.value.peek();
    status.query = this._query.peek();
    status.problems = this._problems.peek();
    status.mode = this._advanced() ? 'advanced' : 'simple';
    return status;
  }

  /** The tree with an id on every node: the same object when it has them, a fresh-id clone of a
   * literal (a spec's `value`, a template's `root`), an empty group for anything else. */
  static withIds(root: FilterGroup | undefined): FilterGroup {
    if (!root || !Array.isArray(root.nodes))
      return Filters.group('and');
    let ids = true;
    Filters.walk(root, (n) => {
      if (typeof n.id !== 'string')
        ids = false;
    });
    return ids ? root : Filters.clone(root, true);
  }

  protected createEditor(): HTMLElement {
    const o = this.options;
    this._mode = o.mode instanceof Signal ? o.mode : signal(o.mode ?? 'simple');
    this._query = computed(() => Filters.format(this.value.value));
    this._problems = computed(() => Filters.validate(this.value.value, o.schema, o.target, o.template));
    this._rows = new Map();
    this._scopes = new Map();
    this._joins = new Map();
    this._horizontal = o.orientation === 'horizontal';
    this._warned = false;
    this._host = {
      schema: o.schema,
      editors: o.editors ?? FilterBuilder.defaultEditors,
      allowAdd: o.template?.allowAdd !== false,
      properties: () => this._properties(),
      operators: (prop) => this._operators(prop),
      lockOf: (id) => this._lockOf(id),
      advanced: () => this._advanced(),
      change: (id, patch) => this._change(id, patch),
      add: (parentId, kind) => kind === 'group' ? this.addGroup(parentId) : this.addCondition(parentId),
      remove: (id) => this.removeNode(id),
    };

    // the first render adopts an id-less literal; until then the root row must not see it
    this._rootRow = new FilterGroupRow(FilterBuilder.withIds(this.value.peek()), this._host, true);
    this._modeLink = link(computed(() => this._mode.value === 'simple' ? 'Advanced' : 'Simple'),
      () => this.setMode(this._mode.peek() === 'simple' ? 'advanced' : 'simple'));
    this._modeLink.classList.add('u2-fb-mode');
    this._modeLink.dataset.u2Part = 'mode';
    this._modeLink.hidden = o.template?.allowAdvanced === false;
    this._rootRow.header.append(this._modeLink);
    this._none = this._host.allowAdd ? link('No conditions — add one', () => this.addCondition()) :
      span('No conditions');
    this._none.classList.add('u2-fb-none');
    this._none.dataset.u2Part = 'none';
    this._hintText = span('', 'u2-fb-hint-text');
    this._flatten = link('Flatten', () => {
      const flat = Filters.flatten(this.value.peek());
      if (flat) {
        this.value.value = flat;
        this._mode.value = 'simple';
      }
    });
    this._flatten.dataset.u2Part = 'flatten';
    this._hintMode = link('Advanced', () => this.setMode('advanced'));
    this._hintMode.dataset.u2Part = 'hint-mode';
    this._hint = div([this._hintText, this._flatten, this._hintMode], 'u2-fb-hint');
    this._hint.dataset.u2Part = 'hint';
    this._hint.setAttribute('role', 'status');
    this._hint.hidden = true;

    this._panel = div([this._rootRow, this._hint], 'u2-fb');
    if (o.showQuery) {
      this._queryEl = span(this._query, 'u2-fb-query');
      this._queryEl.dataset.u2Part = 'query';
      this._panel.append(this._queryEl);
    }
    this._drag = new FilterDragLayer({
      root: this._panel,
      tree: () => this.value.peek(),
      hits: () => this._hits(),
      accepts: (movingId, parentId) => this._accepts(movingId, parentId),
      drop: (id, target) => this.value.value = Filters.move(this.value.peek(), id, target.parentId, target.index),
    });
    this.own(() => this._drag.dispose());

    this.addValidator(() => this._problems.value[0]?.message ?? null);
    this.effect(() => this._render());
    return this._panel;
  }

  /** The mode as rendered: a template that keeps the builder simple wins over the signal. */
  private _advanced(): boolean {
    return this._mode.value === 'advanced' && this.options.template?.allowAdvanced !== false;
  }

  /** A patch key set to `undefined` leaves the node (`not`, a cleared value) rather than staying
   * on it as a key. */
  private _change(id: string, patch: Partial<FilterCondition> | Partial<FilterGroup>): void {
    const root = this.value.peek();
    const node: Record<string, unknown> = {...Filters.find(root, id), ...patch};
    for (const key of Object.keys(patch)) {
      if (node[key] === undefined)
        delete node[key];
    }
    this.value.value = Filters.replace(root, id, node as unknown as FilterNode);
  }

  private _accepts(movingId: string, parentId: string): boolean {
    if (this._lockOf(movingId) !== 'none' || this._lockOf(parentId) !== 'none')
      return false;
    return this._host.allowAdd || Filters.parentOf(this.value.peek(), movingId)?.id === parentId;
  }

  private _append(root: FilterGroup, node: FilterNode, parentId?: string, index?: number): void {
    const parent = (parentId && Filters.find(root, parentId)) || root;
    const into = Filters.isGroup(parent) ? parent : root;
    this.value.value = Filters.insert(root, into.id, index ?? into.nodes.length, node);
  }

  private _properties(): FilterProperty[] {
    const allowed = this.options.template?.allowedProperties;
    const props = this.options.schema.properties;
    return allowed ? props.filter((p) => allowed.includes(p.name)) : props;
  }

  private _operators(prop: FilterProperty) {
    const ops = Filters.operators.for(prop);
    const allowed = this.options.template?.allowedOperators?.[prop.name];
    return allowed ? ops.filter((o) => allowed.includes(o.id)) : ops;
  }

  /** The strongest lock on the node and its ancestors, in the value and in the template. */
  private _lockOf(id: string): Lock {
    const own = Filters.lockOf(this.value.peek(), id);
    const template = this.options.template?.root;
    const inherited = template ? Filters.lockOf(template, id) : 'none';
    return own === 'all' || inherited === 'all' ? 'all' : own === 'value' || inherited === 'value' ? 'value' : 'none';
  }

  private _hits(): FilterDropHit[] {
    const hits: FilterDropHit[] = Array.from(this._panel.querySelectorAll<HTMLElement>('[data-u2-node]'))
      .map((el) => ({id: el.dataset.u2Node!, rect: el.getBoundingClientRect()}));
    hits.push({id: this.value.peek().id, rect: this._rootRow.rows.getBoundingClientRect()});
    return hits;
  }

  /** The reason the mode cannot change: a nested tree (flatten offered when the connectors
   * allow, advanced mode when simple mode is showing it), or a template that keeps the builder
   * simple (`root` null). */
  private _showHint(root: FilterGroup | null): void {
    const canFlatten = root !== null && Filters.flatten(root) !== null;
    const canAdvance = root !== null && !this._advanced() && this.options.template?.allowAdvanced !== false;
    this._hintText.textContent = root === null ? 'This filter cannot use advanced mode.' :
      canFlatten ? 'This filter has nested groups. ' : 'This filter has nested groups that cannot be flattened. ';
    this._flatten.hidden = !canFlatten;
    this._hintMode.hidden = !canAdvance;
    this._hint.hidden = false;
  }

  private _render(): void {
    const value = this.value.value;
    // a bound literal without ids (a spec's `value`, a fresh state) is adopted with ids first
    const root = FilterBuilder.withIds(value);
    if (root !== value) {
      Input.system(() => this.value.value = root);
      return;
    }
    const advanced = this._advanced();
    const horizontal = this._horizontal && !advanced;
    if (advanced && this._horizontal && !this._warned) {
      this._warned = true;
      console.warn('u2: FilterBuilder ignores orientation in advanced mode');
    }
    if (!advanced && !Filters.isFlat(root))
      this._showHint(root);
    else
      this._hint.hidden = true;
    this._panel.classList.toggle('u2-fb-horizontal', horizontal);
    this._panel.classList.toggle('u2-fb-advanced', advanced);
    this._rootRow.header.hidden = horizontal;
    this._rootRow.update(root);
    if (!horizontal && this._rootRow.add.parentElement !== this._rootRow.header)
      this._rootRow.header.insertBefore(this._rootRow.add, this._rootRow.addGroup);

    const seen = new Set<string>();
    this._renderChildren(root, this._rootRow, advanced, horizontal, this._problems.value, seen);
    for (const id of this._rows.keys()) {
      if (!seen.has(id))
        this._releaseRow(id);
    }
    for (const join of this._joins.values())
      join.textContent = root.op;
  }

  private _renderChildren(group: FilterGroup, into: FilterGroupRow, advanced: boolean, horizontal: boolean,
    problems: FilterProblem[], seen: Set<string>): void {
    const order: HTMLElement[] = [];
    for (const node of group.nodes) {
      const kind: RowKind = !Filters.isGroup(node) ? 'cond' : advanced ? 'group' : 'nested';
      let row = this._rows.get(node.id);
      if (row && FilterBuilder._kindOf(row) !== kind) {
        this._releaseRow(node.id);
        row = undefined;
      }
      if (!row)
        row = this._mountRow(kind, node);
      else if (row instanceof FilterRow)
        row.update(node as FilterCondition);
      else
        row.update(node as FilterGroup);
      row.setProblem(problems.find((p) => p.nodeId === node.id)?.message ?? null);
      seen.add(node.id);
      if (horizontal && order.length > 0)
        order.push(this._join(node.id));
      order.push(row.root);
      if (row instanceof FilterGroupRow)
        this._renderChildren(node as FilterGroup, row, advanced, false, problems, seen);
    }
    if (order.length === 0 && into === this._rootRow)
      order.push(this._none);
    if (horizontal && !this._rootRow.add.hidden)
      order.push(this._rootRow.add);
    FilterBuilder._place(into.rows, order);
  }

  private _mountRow(kind: RowKind, node: FilterNode): Row {
    const scope = new Scope();
    const disposer = () => scope.dispose();
    this.own(disposer);
    scope.own(() => this.scope.disown(disposer));
    this._scopes.set(node.id, scope);
    const row = Scope.runWith(scope, () => kind === 'cond' ? new FilterRow(node as FilterCondition, this._host) :
      kind === 'group' ? new FilterGroupRow(node as FilterGroup, this._host) :
        new FilterNestedRow(node as FilterGroup, this._host));
    this._rows.set(node.id, row);
    return row;
  }

  private _releaseRow(id: string): void {
    this._scopes.get(id)?.dispose();
    this._scopes.delete(id);
    this._rows.delete(id);
    this._joins.delete(id);
  }

  private static _kindOf(row: Row): RowKind {
    return row instanceof FilterRow ? 'cond' : row instanceof FilterGroupRow ? 'group' : 'nested';
  }

  /** Puts `desired` under `host` in order, moving only what is out of place — an element already
   * at its index is never detached, so a focused editor keeps its focus. */
  private static _place(host: HTMLElement, desired: HTMLElement[]): void {
    for (let i = 0; i < desired.length; i++) {
      const el = desired[i];
      if (host.children[i] !== el)
        host.insertBefore(el, host.children[i] ?? null);
    }
    while (host.children.length > desired.length)
      host.children[host.children.length - 1].remove();
  }

  /** The connector chip shown before the row `id` in horizontal mode; it flips the group op. */
  private _join(id: string): HTMLElement {
    let join = this._joins.get(id);
    if (!join) {
      const scope = this._scopes.get(id)!;
      join = Scope.runWith(scope, () => button('', () => {
        const root = this.value.peek();
        this.value.value = Filters.update(root, root.id, {op: root.op === 'and' ? 'or' : 'and'});
      }));
      join.classList.add('u2-fb-join');
      join.setAttribute('aria-label', 'Switch and/or');
      Tooltip.bind(join, 'Switch and/or', scope);
      this._joins.set(id, join);
    }
    return join;
  }
}
