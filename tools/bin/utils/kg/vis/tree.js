// A tri-state checkbox tree over a type hierarchy: branches reflect their leaves, leaves carry the state,
// counts show what the generation holds and what the view shows.

/** The synthetic row above the roots. */
const ALL = '*';

export class Tree {
  /**
   * items: [{name, parent, abstract, color, line}] in display order (parents before children);
   * enabled: the Set of leaf names the tree edits in place.
   */
  constructor(container, items, enabled, onChange) {
    this.container = container;
    this.enabled = enabled;
    this.onChange = onChange;
    this.rows = new Map();
    this.items = new Map(items.map((i) => [i.name, i]));
    this.children = new Map();
    for (const item of items) {
      if (!this.children.has(item.parent)) this.children.set(item.parent, []);
      this.children.get(item.parent).push(item);
    }
    this.renderAll(container);
    for (const item of this.children.get(undefined) ?? []) this.render(item, container, 0);
  }

  /** One row above the roots that ticks or clears every leaf of the tree. */
  renderAll(parent) {
    const row = document.createElement('div');
    row.className = 'kgv-node kgv-branch kgv-all';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.title = 'Tick or clear everything in this tree';
    const name = document.createElement('span');
    name.className = 'kgv-name';
    name.textContent = 'all';
    const n = document.createElement('span');
    n.className = 'kgv-n';
    row.append(document.createElement('span'), box, document.createElement('span'), name, n);
    row.firstChild.className = 'kgv-toggle';
    row.children[2].className = 'kgv-swatch';
    row.children[2].style.visibility = 'hidden';
    parent.append(row);
    this.rows.set(ALL, {item: {name: ALL}, row, box, n, children: null});
    box.addEventListener('change', () => {
      this.set(ALL, box.checked);
      this.refresh();
      this.onChange();
    });
  }

  render(item, parent, depth) {
    const row = document.createElement('div');
    row.className = 'kgv-node' + (item.abstract ? ' kgv-branch' : '');
    const toggle = document.createElement('span');
    toggle.className = 'kgv-toggle';
    const kids = this.children.get(item.name) ?? [];
    toggle.textContent = kids.length ? '▾' : '';
    const box = document.createElement('input');
    box.type = 'checkbox';
    const swatch = document.createElement('span');
    swatch.className = 'kgv-swatch' + (item.line ? ' kgv-line' : '');
    swatch.style.background = item.color;
    swatch.title = `Only ${item.name}${this.children.get(item.name)?.length ? ' and what is under it' : ''}; click again to restore the previous selection`;
    swatch.addEventListener('click', () => this.solo(item.name));
    const name = document.createElement('span');
    name.className = 'kgv-name';
    name.textContent = item.name;
    name.title = item.description ?? item.name;
    const n = document.createElement('span');
    n.className = 'kgv-n';
    row.append(toggle, box, swatch, name, n);
    parent.append(row);
    let children = null;
    if (kids.length) {
      children = document.createElement('div');
      children.className = 'kgv-children';
      parent.append(children);
      for (const kid of kids) this.render(kid, children, depth + 1);
      if (depth >= 1) {
        children.hidden = true;
        toggle.textContent = '▸';
      }
      toggle.addEventListener('click', () => {
        children.hidden = !children.hidden;
        toggle.textContent = children.hidden ? '▸' : '▾';
      });
    }
    this.rows.set(item.name, {item, row, box, n, children});
    box.addEventListener('change', () => {
      this.set(item.name, box.checked);
      this.refresh();
      this.onChange();
    });
  }

  /** Only this subtree on; the second click on the same dot brings the previous selection back. */
  solo(name) {
    if (this.soloed?.name === name) {
      this.enabled.clear();
      for (const l of this.soloed.previous) this.enabled.add(l);
      this.soloed = null;
    }
    else {
      const previous = this.soloed?.previous ?? new Set(this.enabled);
      this.enabled.clear();
      for (const l of this.leaves(name)) this.enabled.add(l);
      this.soloed = {name, previous};
    }
    this.refresh();
    this.onChange('solo');
  }

  /** Enables or disables every leaf under [name], itself included when it is one. */
  set(name, on) {
    for (const leaf of this.leaves(name))
      if (on) this.enabled.add(leaf);
      else this.enabled.delete(leaf);
  }

  /** The togglable types under [name]: a concrete type counts as its own leaf even when it has subtypes, an
   * abstract one only through them, and a type with nothing in the generation is not counted at all. */
  leaves(name) {
    const item = this.items.get(name);
    if (item?.empty) return [];
    const kids = name === ALL ? this.children.get(undefined) : this.children.get(name);
    const own = item && !item.abstract ? [name] : [];
    return [...own, ...(kids ?? []).flatMap((k) => this.leaves(k.name))];
  }

  /** Reflects `enabled` and the counts: `total` per leaf from the generation, `visible` per leaf from the view. */
  refresh(total, visible) {
    if (total) this.total = total;
    if (visible) this.visible = visible;
    for (const [name, r] of this.rows) {
      const leaves = this.leaves(name);
      const on = leaves.filter((l) => this.enabled.has(l)).length;
      r.box.checked = on === leaves.length && on > 0;
      r.box.indeterminate = on > 0 && on < leaves.length;
      const all = leaves.reduce((s, l) => s + (this.total?.get(l) ?? 0), 0);
      const shown = leaves.reduce((s, l) => s + (this.visible?.get(l) ?? 0), 0);
      r.n.textContent = all ? (shown === all ? `${all}` : `${shown} / ${all}`) : '';
      r.row.classList.toggle('kgv-empty', !all);
      r.row.classList.toggle('kgv-solo', this.soloed?.name === name);
    }
  }
}
