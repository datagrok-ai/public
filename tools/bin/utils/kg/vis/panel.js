// The right pane: one node with its edge groups and operations, one edge, or a Cypher answer.
// Every element is built with DOM calls; ids and values never go through innerHTML.

const HIDDEN = new Set(['_id', '_label', '_src', '_dst', 'id', 'name', 'type', 'types', 'batch']);
const PATH_KEYS = new Set(['home', 'path', 'file', 'evidence']);
const OPS = [['explain', 'Explain'], ['impact', 'Impact'], ['tests-for', 'Tests for']];

export function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') node.className = v;
    else if (k === 'on') for (const [ev, fn] of Object.entries(v)) node.addEventListener(ev, fn);
    else if (v !== undefined && v !== null) node.setAttribute(k, v);
  }
  for (const c of children.flat())
    if (c !== null && c !== undefined) node.append(c.nodeType ? c : document.createTextNode(String(c)));
  return node;
}

export class Panel {
  /** ctx: {ids, index(id) -> graph index or -1, badge(type), isVisible(i), select(i), pin(i), scope(id), api, repoRoot, showQuery} */
  constructor(container, empty, ctx) {
    this.container = container;
    this.empty = empty;
    this.ctx = ctx;
  }

  clear() {
    this.container.hidden = true;
    this.container.replaceChildren();
    this.empty.hidden = false;
  }

  show(...content) {
    this.empty.hidden = true;
    this.container.hidden = false;
    this.container.replaceChildren(...content);
  }

  async node(id) {
    this.show(el('div', {class: 'kgv-section-empty'}, `Loading ${id}…`));
    let answer;
    try {
      answer = await this.ctx.api.node(id);
    }
    catch (e) {
      this.show(el('div', {class: 'kgv-error'}, e.message));
      return;
    }
    const {node, edges} = answer;
    const index = this.ctx.index(id);
    const ops = el('div');
    const tabs = el('div', {class: 'kgv-tabs'},
      el('button', {class: 'kgv-on', on: {click: (ev) => this.tab(ev.target, ops, () => this.groups(node, edges))}}, 'Edges'),
      ...OPS.map(([op, label]) => el('button', {on: {click: (ev) => this.tab(ev.target, ops, () => this.operation(op, id))}}, label)));
    ops.append(this.groups(node, edges));
    this.show(
      el('h2', {}, node.name ?? id),
      el('div', {class: 'kgv-id-line'}, this.ctx.badge(node.type), el('code', {}, id),
        el('button', {title: 'Copy the id', on: {click: () => navigator.clipboard.writeText(id)}}, '⧉')),
      el('div', {class: 'kgv-meta'}, ...['_label', 'source_layer', 'visibility', 'status', 'provenance'].filter((k) => node[k])
        .map((k) => el('span', {class: 'kgv-badge', title: k}, String(node[k])))),
      el('div', {class: 'kgv-actions'},
        el('button', {disabled: index < 0 || !this.ctx.isVisible(index) ? '' : null, on: {click: () => this.ctx.zoom(index)}}, 'Zoom'),
        el('button', {disabled: index < 0 ? '' : null, title: 'Pull every neighbour into the view', on: {click: () => this.ctx.pin(index)}}, 'Expand'),
        el('button', {disabled: index < 0 ? '' : null, title: 'Start exploring here: this node and its neighbours, then click nodes to open theirs', on: {click: () => this.ctx.explore(index)}}, 'Explore'),
        el('button', {disabled: index < 0 ? '' : null, title: 'Show this node alone; Back restores the view', on: {click: () => this.ctx.isolate([index])}}, 'Isolate'),
        el('button', {title: 'Restrict the view to this node, what is under it and its neighbours', on: {click: () => this.ctx.scope(id)}}, 'Focus'),
        node.home ? this.link(node.home, 'Home') : null,
        node.path && node.path !== node.home ? this.link(node.path, 'Source') : null),
      this.properties(node),
      tabs, ops);
  }

  tab(button, host, render) {
    for (const b of button.parentElement.children) b.classList.toggle('kgv-on', b === button);
    host.replaceChildren(el('div', {class: 'kgv-section-empty'}, 'Loading…'));
    Promise.resolve(render()).then((content) => host.replaceChildren(content));
  }

  properties(node) {
    const rows = Object.entries(node).filter(([k, v]) => !HIDDEN.has(k) && v !== null && v !== undefined && v !== '' && !(Array.isArray(v) && !v.length));
    if (!rows.length) return null;
    return el('dl', {class: 'kgv-props'}, ...rows.flatMap(([k, v]) => [el('dt', {}, k), el('dd', {}, this.value(k, v))]));
  }

  value(key, v) {
    if (Array.isArray(v)) return PATH_KEYS.has(key) ? v.map((p, i) => [i ? ', ' : '', this.link(String(p))]) : v.join(', ');
    if (typeof v === 'object') return JSON.stringify(v);
    const text = String(v);
    if (PATH_KEYS.has(key)) return this.link(text);
    if (/^https?:\/\//.test(text)) return el('a', {href: text, target: '_blank', rel: 'noopener'}, text);
    if (this.ctx.index(text) >= 0) return this.idLink(text);
    return text;
  }

  /** A repo path opens in the editor; the text stays the path so it can be copied. */
  link(p, label) {
    const href = `vscode://file/${this.ctx.repoRoot}/${p}`;
    return el('a', {href, title: `Open ${p}`}, label ?? p);
  }

  idLink(id, name) {
    const index = this.ctx.index(id);
    return el('a', {href: '#', class: index >= 0 && !this.ctx.isVisible(index) ? 'kgv-hidden-node' : '',
      title: index >= 0 && !this.ctx.isVisible(index) ? `${id} (hidden by the filters; opens the pane, Expand pulls it in)` : id,
      on: {click: (ev) => { ev.preventDefault(); this.ctx.select(index >= 0 ? index : null, id); }}}, name ?? id);
  }

  groups(node, edges) {
    if (!edges.length) return el('div', {class: 'kgv-section-empty'}, 'No edges.');
    return el('div', {}, ...edges.map((g) => {
      const body = el('div', {class: 'kgv-group-body'}, ...g.targets.map((s) =>
        el('div', {class: 'kgv-neighbor' + (this.ctx.index(s.id) >= 0 && !this.ctx.isVisible(this.ctx.index(s.id)) ? ' kgv-hidden-node' : '')},
          this.ctx.badge(s.type), this.idLink(s.id, s.name ?? s.id),
          el('span', {class: 'kgv-edge-link', title: 'The edge itself', on: {click: () => this.edge(g.direction === 'out' ? node.id : s.id, g.direction === 'out' ? s.id : node.id, g.edge)}}, 'edge'))));
      if (g.count > g.targets.length)
        body.append(el('div', {class: 'kgv-more', on: {click: () => this.ctx.showQuery(allOf(node.id, g))}}, `… ${g.count - g.targets.length} more: run as a query`));
      body.hidden = g.count > 8;
      const arrow = el('span', {class: 'kgv-arrow'}, body.hidden ? '▸' : '▾');
      const head = el('div', {class: 'kgv-group-head', on: {click: () => { body.hidden = !body.hidden; arrow.textContent = body.hidden ? '▸' : '▾'; }}},
        arrow, el('span', {class: 'kgv-kind'}, g.direction === 'out' ? `${g.edge} →` : `← ${g.edge}`), el('span', {class: 'kgv-n'}, g.count),
        el('span', {class: 'kgv-how', title: `derived by ${g.derived_by.join(', ')}; confidence ${confidence(g)}`}, `${g.derived_by.join(', ')} ${confidence(g)}`));
      return el('div', {class: 'kgv-group'}, head, body);
    }));
  }

  async edge(from, to, kind) {
    let answer;
    try {
      answer = await this.ctx.api.edge(from, to, kind);
    }
    catch (e) {
      this.show(el('div', {class: 'kgv-error'}, e.message));
      return;
    }
    const e = answer.edge;
    const ends = [this.ctx.index(from), this.ctx.index(to)];
    this.show(
      el('h2', {}, kind),
      el('div', {class: 'kgv-id-line'}, this.idLink(from), el('span', {}, '→'), this.idLink(to)),
      el('div', {class: 'kgv-actions'},
        el('button', {disabled: ends.some((i) => i < 0) ? '' : null, title: 'Show the two ends and this edge alone; Back restores the view',
          on: {click: () => this.ctx.isolate(ends)}}, 'Isolate')),
      el('div', {class: 'kgv-meta'}, ...['derived_by', 'confidence'].filter((k) => e[k] !== null && e[k] !== undefined)
        .map((k) => el('span', {class: 'kgv-badge', title: k}, `${k} ${e[k]}`))),
      this.properties(Object.fromEntries(Object.entries(e).filter(([k]) => k !== 'derived_by' && k !== 'confidence'))) ?? el('div', {class: 'kgv-section-empty'}, 'No properties.'));
  }

  async operation(op, id) {
    let result;
    try {
      result = await this.ctx.api.op(op, id);
    }
    catch (e) {
      return el('div', {class: 'kgv-error'}, e.message);
    }
    return this.opsResult(result);
  }

  opsResult(result) {
    return el('div', {},
      result.notes?.length ? el('div', {class: 'kgv-notes'}, ...result.notes.map((n) => el('div', {}, n))) : null,
      ...result.sections.map((s) => el('div', {},
        el('div', {class: 'kgv-section-title'}, `${s.title}${s.total !== undefined ? ` (${s.rows.length < s.total ? `${s.rows.length} of ` : ''}${s.total})` : ''}`),
        s.rows.length ? this.table(Object.keys(s.rows[0]).filter((k) => k !== 'path'), s.rows) : el('div', {class: 'kgv-section-empty'}, s.empty ?? 'none'))));
  }

  /** A table whose id-shaped cells select the node. */
  table(columns, rows) {
    return el('div', {class: 'kgv-table'}, el('table', {},
      el('thead', {}, el('tr', {}, ...columns.map((c) => el('th', {}, c)))),
      el('tbody', {}, ...rows.map((r) => el('tr', {}, ...columns.map((c) => {
        const v = r[c];
        const text = v === null || v === undefined ? '' : typeof v === 'object' ? JSON.stringify(v) : String(v);
        const index = typeof v === 'string' ? this.ctx.index(v) : -1;
        return index >= 0
          ? el('td', {class: 'kgv-id', title: text, on: {click: () => this.ctx.select(index, v)}}, text)
          : el('td', {title: text}, text);
      }))))));
  }
}

/** `1`, `0.8–1`, or nothing when no edge of the group carries a confidence. */
function confidence(g) {
  return g.confidence === null ? '' : g.confidence[0] === g.confidence[1] ? String(g.confidence[0]) : `${g.confidence[0]}–${g.confidence[1]}`;
}

/** The Cypher that lists one edge group in full. */
function allOf(id, g) {
  const pattern = g.direction === 'out' ? `(n)-[e:\`${g.edge}\`]->(m)` : `(n)<-[e:\`${g.edge}\`]-(m)`;
  return `MATCH ${pattern} WHERE n.id = '${id.replace(/\\/g, '\\\\').replace(/'/g, "\\'")}'\nRETURN m.id AS id, m.type AS type, m.name AS name, e.confidence AS confidence, e.derived_by AS derived_by\nORDER BY id LIMIT 1000`;
}
