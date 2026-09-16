// The page: loads the render tier, wires the three panes together and keeps the view in the URL hash.
import {loadAll, api} from './data.js';
import {Filter, PRESETS, compact} from './filter.js';
import {Tree} from './tree.js';
import {GraphView} from './graph.js';
import {Panel, el} from './panel.js';

const ROOT_COLORS = {feature: '#2f6fdb', concept: '#b3661a', component: '#2e7d32', artifact: '#7a4fb3', work: '#c2413b',
  actor: '#0e8a8a', infra: '#6b7280', type: '#b5179e'};
const FACETS = [['layer', 'layer'], ['visibility', 'visibility'], ['status', 'status'], ['provenance', 'provenance'], ['derivedBy', 'edge derived by']];
const $ = (id) => document.getElementById(id);

async function main() {
  const [data, questions] = await Promise.all([loadAll(), api.questions()]);
  const style = makeStyle(data);
  const filter = new Filter(data);
  const idIndex = new Map(data.ids.map((id, i) => [id, i]));
  const index = (id) => idIndex.get(String(id).replace(/^~/, '')) ?? -1;
  const positions = seedPositions(data, style);
  const lower = {ids: data.ids.map((s) => s.toLowerCase()), names: data.names.map((s) => s.toLowerCase())};
  const state = {preset: 'spine', selection: null, isolated: null};
  const viewStack = [];
  let masks = null;
  let view = null;
  let highlighted = null;

  const graph = new GraphView($('graph'), $('labels'), {
    names: data.names,
    degree: data.degree,
    isolated: () => state.isolated,
    onSimulation: () => syncSim(),
    onSelect: (i) => {
      if (i === null) {
        panel.clear();
        state.selection = null;
        writeHash();
        return;
      }
      state.selection = data.ids[i];
      if (filter.explore && !filter.explore.expanded.has(i)) {
        graph.savePositions(positions);
        seedAround(i);
        filter.explore.open.add(i);
        filter.explore.expanded.add(i);
        refresh({resimulate: true});
      }
      panel.node(data.ids[i]);
      writeHash();
    },
    onHover: (i, event) => tooltip(i, event),
    onExpand: (i) => expand(i),
    onIsolate: (i) => { if (state.isolated === i) popView(); else isolate([i]); },
  });

  const panel = new Panel($('pane'), $('pane-empty'), {
    api, index, repoRoot: data.manifest.repo_root,
    badge: (type) => badge(type, style),
    isVisible: (i) => graph.isVisible(i),
    select: (i, id) => {
      if (i !== null && i >= 0) {
        if (filter.explore) {
          graph.savePositions(positions);
          seedAround(i);
          filter.explore.open.add(i);
          filter.explore.expanded.add(i);
          refresh({resimulate: true});
        }
        graph.select(i);
        graph.zoomTo(i);
        state.selection = data.ids[i];
      }
      panel.node(id ?? data.ids[i]);
      writeHash();
    },
    zoom: (i) => graph.zoomTo(i),
    pin: (i) => expand(i),
    explore: (i) => explore(i),
    isolate: (indices) => isolate(indices),
    scope: (id) => setScope(id),
    showQuery: (cypher) => {
      $('drawer').hidden = false;
      $('cypher').value = cypher;
      runQuery();
    },
  });

  // header
  const m = data.manifest;
  $('batch').textContent = `${m.gen} · built ${new Date(m.built_at).toLocaleString()}`;
  $('batch').title = `revisions ${Object.entries(m.revisions).map(([k, v]) => `${k} ${String(v).slice(0, 10)}`).join(', ')}`;
  const problems = Object.entries(m.problems ?? {}).filter(([, n]) => n > 0);
  $('chips').replaceChildren(...(m.notes ?? []).map((n) => el('span', {class: 'kgv-chip', title: n}, n.split(':')[0])));
  if (problems.length)
    $('chips').append(el('span', {class: 'kgv-chip kgv-chip-info', title: problems.map(([k, n]) => `${k}: ${n}`).join('\n')}, `${problems.length} problem kinds`));
  if (data.dropped)
    $('chips').append(el('span', {class: 'kgv-chip', title: `${data.dropped} edges point at ids the generation has no node for`}, `${data.dropped} dangling`));

  // trees
  const totals = counts(data.nodeTypes, data.type);
  const kindTotals = counts(data.edgeKinds, data.kind);
  const typeTree = new Tree($('node-types'), nodeTypeItems(data, style), filter.nodeTypes, (reason) => {
    state.preset = 'custom';
    // the solo has already replaced the type set; the view to come back to is the one before it
    if (reason === 'solo' && leaveExplore() && typeTree.soloed) viewStack[viewStack.length - 1].filter.nodeTypes = new Set(typeTree.soloed.previous);
    syncPresets();
    refresh({resimulate: true});
  });
  const kindTree = new Tree($('edge-kinds'), edgeKindItems(data, style), filter.edgeKinds, () => refresh());
  $('node-type-count').textContent = `${data.nodeTypes.length} of ${data.schema.nodeTypes.filter((t) => !t.abstract).length} populated`;
  $('edge-kind-count').textContent = `${data.edgeKinds.length} kinds`;

  // presets
  for (const [key, p] of Object.entries(PRESETS))
    $('presets').append(el('button', {'data-preset': key, title: p.types ? p.types.join(', ') : 'every node of the generation', on: {click: () => applyPreset(key)}}, p.label));
  function applyPreset(key) {
    const p = PRESETS[key];
    if (p.confirm && data.nodes > p.confirm && !confirm(`${data.nodes.toLocaleString()} nodes and ${data.edges.toLocaleString()} edges: show everything?`)) return;
    filter.nodeTypes.clear();
    for (const t of p.types ?? data.nodeTypes) filter.nodeTypes.add(t);
    state.preset = key;
    leaveExplore();
    syncPresets();
    refresh({resimulate: true});
  }
  /** A preset or a solo asks for the whole type, not the explored corner of it; Back still returns there. */
  function leaveExplore() {
    if (!filter.explore) return false;
    pushView();
    filter.explore = null;
    state.isolated = null;
    return true;
  }
  function syncPresets() {
    for (const b of $('presets').children) b.classList.toggle('kgv-on', b.dataset.preset === state.preset);
  }

  // questions
  for (const q of questions)
    $('question').append(el('option', {value: q.id}, `${q.question}${q.status === 'blocked' ? ' (blocked)' : ''}`));
  $('question').addEventListener('change', () => showQuestion(questions.find((q) => q.id === $('question').value)));
  $('ask').addEventListener('click', askQuestion);
  function showQuestion(q) {
    $('question-params').replaceChildren();
    $('question-why').replaceChildren();
    $('question-status').textContent = '';
    $('ask').hidden = !q;
    if (!q) return;
    for (const [name, spec] of Object.entries(q.params ?? {}))
      $('question-params').append(el('label', {for: `qp-${name}`, title: spec.description ?? ''}, name),
        el('input', {id: `qp-${name}`, type: 'text', value: spec.default ?? '', title: spec.description ?? '', 'data-param': name,
          on: {keydown: (ev) => { if (ev.key === 'Enter') askQuestion(); }}}));
    $('question-why').append(q.why);
    if (q.status === 'blocked') $('question-why').append(el('span', {class: 'kgv-blocked'}, `Blocked: ${q.blocked_by}`));
  }
  async function askQuestion() {
    const q = questions.find((x) => x.id === $('question').value);
    if (!q) return;
    const params = Object.fromEntries([...$('question-params').querySelectorAll('input')].map((i) => [i.dataset.param, i.value.trim()]));
    $('question-status').textContent = 'asking…';
    applyView(q.view);
    $('drawer').hidden = false;
    $('cypher').value = `// ${q.id}${Object.keys(params).length ? ' ' + Object.entries(params).map(([k, v]) => `${k}=${v}`).join(' ') : ''}\n${q.cypher.trim()}`;
    try {
      const r = await api.ask(q.id, params);
      $('question-status').textContent = `${r.rows.length} row${r.rows.length === 1 ? '' : 's'} in ${r.ms} ms`;
      showAnswer(r, q.highlight);
    }
    catch (e) {
      $('question-status').textContent = 'error';
      $('query-status').textContent = 'error';
      $('query-result').replaceChildren(el('div', {class: 'kgv-error'}, e.message));
    }
  }
  /** A question's view: the preset or the types to show while its answer is lit, and whether only the answer stays. */
  function applyView(view) {
    if (!view) return;
    if (view.preset && PRESETS[view.preset]) {
      filter.nodeTypes.clear();
      for (const t of PRESETS[view.preset].types ?? data.nodeTypes) filter.nodeTypes.add(t);
      state.preset = view.preset;
    }
    if (view.types) {
      filter.nodeTypes.clear();
      for (const t of view.types) filter.nodeTypes.add(t);
      state.preset = 'custom';
    }
    if (view.edges) {
      filter.edgeKinds.clear();
      for (const k of view.edges) filter.edgeKinds.add(k);
    }
    leaveExplore();
    filter.pinned.clear();
    $('highlight').checked = true;
    $('restrict').checked = view.restrict === true;
    syncPresets();
  }
  function showAnswer(r, highlightColumns) {
    $('query-status').textContent = `${r.rows.length} row${r.rows.length === 1 ? '' : 's'}${r.truncated ? ' (truncated)' : ''} in ${r.ms} ms`;
    $('query-result').replaceChildren(
      r.notes?.length ? el('div', {class: 'kgv-notes'}, ...r.notes.map((n) => el('div', {}, n))) : null,
      r.rows.length ? panel.table(r.columns, r.rows) : el('div', {class: 'kgv-section-empty'}, 'no rows'));
    hits = new Set();
    for (const row of r.rows)
      for (const [column, v] of Object.entries(row))
        if (!highlightColumns || highlightColumns.includes(column)) collectIds(v, hits);
    applyQueryHits();
  }

  // facets
  for (const [facet, label] of FACETS) {
    const values = data.enums[facet].filter((v) => v !== '');
    if (values.length < 2) continue;
    const row = el('div', {class: 'kgv-facet', 'data-facet': facet}, el('span', {class: 'kgv-facet-name'}, label));
    for (const v of values) {
      const b = el('button', {class: 'kgv-on', 'data-value': v, on: {click: () => {
        if (filter.facets[facet].has(v)) filter.facets[facet].delete(v);
        else filter.facets[facet].add(v);
        b.classList.toggle('kgv-on', filter.facets[facet].has(v));
        refresh();
      }}}, v);
      row.append(b);
    }
    $('facets').append(row);
  }
  function syncFacets() {
    for (const row of $('facets').querySelectorAll('.kgv-facet'))
      for (const b of row.querySelectorAll('button')) b.classList.toggle('kgv-on', filter.facets[row.dataset.facet].has(b.dataset.value));
    $('conf').value = String(filter.minConfidence);
    $('conf-value').textContent = (filter.minConfidence / 100).toFixed(2);
    $('scope').value = filter.scope ? data.ids[filter.scope.index] : '';
    if (filter.scope) $('hops').value = String(filter.scope.hops);
  }
  $('conf').addEventListener('input', () => {
    filter.minConfidence = Number($('conf').value);
    $('conf-value').textContent = (filter.minConfidence / 100).toFixed(2);
    refresh();
  });

  // scope
  $('scope').addEventListener('change', () => setScope($('scope').value.trim()));
  $('hops').addEventListener('change', () => setScope($('scope').value.trim()));
  $('scope-clear').addEventListener('click', () => setScope(''));
  function setScope(id) {
    const i = id ? index(id) : -1;
    $('scope').value = id;
    if (id && i < 0) {
      $('scope').setCustomValidity(`no node ${id}`);
      $('scope').reportValidity();
      return;
    }
    $('scope').setCustomValidity('');
    filter.scope = i < 0 ? null : {index: i, hops: Number($('hops').value)};
    refresh({resimulate: true});
  }

  // explore, isolate and the view stack
  function pushView() {
    viewStack.push({filter: filter.snapshot(), preset: state.preset, isolated: state.isolated});
    syncMode();
  }
  function popView() {
    const s = viewStack.pop();
    if (!s) return;
    filter.restore(s.filter);
    state.preset = s.preset;
    state.isolated = s.isolated;
    syncPresets();
    syncFacets();
    typeTree.enabled = filter.nodeTypes;
    kindTree.enabled = filter.edgeKinds;
    typeTree.soloed = null;
    kindTree.soloed = null;
    refresh({resimulate: true});
  }
  /** Nodes about to appear next to [i] start at its position, so opening a node grows the picture from it. */
  function seedAround(i) {
    const x = positions[i * 2], y = positions[i * 2 + 1];
    for (const j of filter.neighbours(i)) {
      if (graph.isVisible(j)) continue;
      const a = hash(j * 7 + 1) * Math.PI * 2, r = 20 + hash(j * 7 + 2) * 40;
      positions[j * 2] = x + Math.cos(a) * r;
      positions[j * 2 + 1] = y + Math.sin(a) * r;
    }
  }
  function explore(i) {
    pushView();
    graph.savePositions(positions);
    seedAround(i);
    filter.explore = {open: new Set([i]), expanded: new Set([i])};
    filter.scope = null;
    filter.restrict = null;
    state.isolated = null;
    state.selection = data.ids[i];
    syncFacets();
    refresh({resimulate: true});
    panel.node(data.ids[i]);
  }
  function isolate(indices) {
    pushView();
    filter.explore = {open: new Set(indices), expanded: new Set()};
    filter.scope = null;
    filter.restrict = null;
    state.isolated = indices[0];
    state.selection = data.ids[indices[0]];
    syncFacets();
    refresh({resimulate: true});
    panel.node(data.ids[indices[0]]);
  }
  function syncMode() {
    $('back').hidden = !viewStack.length;
    const ex = filter.explore;
    $('mode').hidden = !ex;
    if (ex) $('mode').textContent = ex.expanded.size ? `exploring · ${ex.expanded.size} opened · click a node to open its neighbours` : `${ex.open.size === 1 ? data.names[[...ex.open][0]] : `${ex.open.size} nodes`} alone · click to open neighbours`;
  }
  $('back').addEventListener('click', popView);

  // toolbar
  $('sim').addEventListener('click', () => { graph.toggleSimulation(); syncSim(); });
  function syncSim() {
    $('sim').textContent = graph.simulating ? 'Freeze' : 'Simulate';
  }
  $('fit').addEventListener('click', () => graph.fit());
  $('zoom-sel').addEventListener('click', () => { if (graph.selected !== null) graph.zoomTo(graph.selected); });
  $('unpin').addEventListener('click', () => { filter.pinned.clear(); refresh({resimulate: true}); });
  $('labels-on').addEventListener('change', () => graph.setLabels($('labels-on').checked));
  $('drawer-toggle').addEventListener('click', () => { $('drawer').hidden = !$('drawer').hidden; if (!$('drawer').hidden) $('cypher').focus(); });
  document.addEventListener('keydown', (ev) => {
    if (ev.target.matches('input, textarea, select')) return;
    if (ev.key === ' ') { ev.preventDefault(); $('sim').click(); }
    else if (ev.key === 'f') graph.fit();
    else if (ev.key === 'z') $('zoom-sel').click();
    else if (ev.key === '/') { ev.preventDefault(); $('search').focus(); }
    else if (ev.key === 'Backspace') { ev.preventDefault(); popView(); }
    else if (ev.key === 'Escape') { graph.select(null); panel.clear(); }
  });

  // Cypher drawer
  const queryHistory = JSON.parse(localStorage.getItem('kgv-history') ?? '[]');
  const syncHistory = () => $('history').replaceChildren(el('option', {value: ''}, 'history…'), ...queryHistory.map((h) => el('option', {value: h}, h.slice(0, 80))));
  syncHistory();
  $('history').addEventListener('change', () => { if ($('history').value) $('cypher').value = $('history').value; $('history').value = ''; });
  $('run').addEventListener('click', runQuery);
  $('cypher').addEventListener('keydown', (ev) => { if (ev.key === 'Enter' && (ev.ctrlKey || ev.metaKey)) { ev.preventDefault(); runQuery(); } });
  $('highlight').addEventListener('change', () => applyQueryHits());
  $('restrict').addEventListener('change', () => applyQueryHits());
  let hits = null;
  async function runQuery() {
    const cypher = $('cypher').value.trim();
    if (!cypher) return;
    $('query-status').textContent = 'running…';
    try {
      const r = await api.query(cypher, 1000);
      showAnswer(r);
      if (!queryHistory.includes(cypher)) {
        queryHistory.unshift(cypher);
        queryHistory.splice(20);
        localStorage.setItem('kgv-history', JSON.stringify(queryHistory));
        syncHistory();
      }
    }
    catch (e) {
      $('query-status').textContent = 'error';
      $('query-result').replaceChildren(el('div', {class: 'kgv-error'}, e.message));
    }
  }
  function collectIds(v, into) {
    if (typeof v === 'string') {
      const i = index(v);
      if (i >= 0) into.add(i);
    }
    else if (Array.isArray(v)) for (const x of v) collectIds(x, into);
    else if (v && typeof v === 'object') for (const x of Object.values(v)) collectIds(x, into);
  }
  /** A hit the current filter hides is pinned into the view, so an answer is never lit on nothing; the next answer
   * takes those pins back. */
  let answerPins = new Set();
  function applyQueryHits() {
    const restrict = $('restrict').checked && hits ? hits : null;
    let changed = (restrict ?? null) !== (filter.restrict ?? null);
    filter.restrict = restrict;
    for (const i of answerPins)
      if (!hits?.has(i)) {
        filter.pinned.delete(i);
        changed = true;
      }
    answerPins = new Set();
    if ($('highlight').checked && hits)
      for (const i of hits)
        if (!masks.nodeOn[i] && !filter.pinned.has(i)) {
          filter.pinned.add(i);
          answerPins.add(i);
          changed = true;
        }
    if (changed) refresh({resimulate: true});
    highlighted = $('highlight').checked && hits?.size ? hits : null;
    graph.highlight(highlighted);
  }

  // search
  let results = [];
  let active = -1;
  $('search').addEventListener('input', () => search($('search').value.trim().toLowerCase()));
  $('search').addEventListener('focus', () => { if (results.length) $('results').hidden = false; });
  $('search').addEventListener('keydown', (ev) => {
    if (ev.key === 'ArrowDown') { active = Math.min(active + 1, results.length - 1); renderResults(); ev.preventDefault(); }
    else if (ev.key === 'ArrowUp') { active = Math.max(active - 1, 0); renderResults(); ev.preventDefault(); }
    else if (ev.key === 'Enter' && results[active >= 0 ? active : 0] !== undefined) pick(results[active >= 0 ? active : 0]);
    else if (ev.key === 'Escape') { $('results').hidden = true; $('search').blur(); }
  });
  document.addEventListener('click', (ev) => { if (!ev.target.closest('.kgv-search')) $('results').hidden = true; });
  function search(q) {
    results = [];
    active = -1;
    if (q.length < 2) {
      $('results').hidden = true;
      return;
    }
    const exact = [], starts = [], contains = [];
    for (let i = 0; i < data.nodes && contains.length < 200; i++) {
      const name = lower.names[i], id = lower.ids[i];
      if (name === q || id === q) exact.push(i);
      else if (name.startsWith(q) || id.startsWith(q)) starts.push(i);
      else if (name.includes(q) || id.includes(q)) contains.push(i);
    }
    results = [...exact, ...starts, ...contains].slice(0, 40);
    renderResults();
  }
  function renderResults() {
    $('results').hidden = !results.length;
    $('results').replaceChildren(...results.map((i, k) => el('div', {class: 'kgv-result' + (k === active ? ' kgv-active' : '') + (graph.isVisible(i) ? '' : ' kgv-hidden-node'),
      on: {click: () => pick(i)}}, badge(data.nodeTypes[data.type[i]], style), el('span', {}, data.names[i]), el('span', {class: 'kgv-mono'}, data.ids[i]))));
  }
  function pick(i) {
    $('results').hidden = true;
    if (filter.explore) {
      graph.savePositions(positions);
      seedAround(i);
      filter.explore.open.add(i);
      filter.explore.expanded.add(i);
      refresh({resimulate: true});
    }
    else if (!graph.isVisible(i)) {
      filter.pinned.add(i);
      refresh({resimulate: true});
    }
    graph.select(i);
    graph.zoomTo(i);
    panel.node(data.ids[i]);
    state.selection = data.ids[i];
    writeHash();
  }

  function expand(i) {
    graph.savePositions(positions);
    seedAround(i);
    if (filter.explore) {
      filter.explore.open.add(i);
      filter.explore.expanded.add(i);
    }
    else {
      for (const j of filter.neighbours(i)) filter.pinned.add(j);
      filter.pinned.add(i);
    }
    refresh({resimulate: true});
    graph.select(i);
  }

  function tooltip(i, event) {
    const t = $('tooltip');
    if (i === null || !event) {
      t.hidden = true;
      return;
    }
    t.replaceChildren(el('div', {}, badge(data.nodeTypes[data.type[i]], style), ' ', data.names[i]), el('span', {class: 'kgv-mono'}, data.ids[i]));
    const rect = $('graph').getBoundingClientRect();
    t.style.left = `${Math.min(event.clientX - rect.left + 12, rect.width - 300)}px`;
    t.style.top = `${event.clientY - rect.top + 12}px`;
    t.hidden = false;
  }

  function refresh(options = {}) {
    graph.savePositions(positions);
    masks = filter.apply();
    view = compact(data, masks, style, positions);
    graph.setData(view, options.resimulate === true);
    if (state.selection) graph.select(index(state.selection));
    if (highlighted) graph.highlight(highlighted);
    typeTree.refresh(totals, counts(data.nodeTypes, data.type, masks.nodeOn));
    kindTree.refresh(kindTotals, new Map(data.edgeKinds.map((k, i) => [k, masks.kindCounts[i]])));
    $('visible').textContent = `${masks.nodeCount.toLocaleString()} nodes · ${masks.edgeCount.toLocaleString()} edges`;
    $('stats').textContent = `${data.nodes.toLocaleString()} nodes · ${data.edges.toLocaleString()} edges`;
    $('pin-count').textContent = filter.pinned.size ? `(${filter.pinned.size})` : '';
    syncSim();
    syncMode();
    writeHash();
  }

  // URL state
  function writeHash() {
    const p = new URLSearchParams();
    p.set('p', state.preset);
    if (state.preset === 'custom') p.set('t', [...filter.nodeTypes].join(','));
    const kindsOff = data.edgeKinds.filter((k) => !filter.edgeKinds.has(k));
    if (kindsOff.length) p.set('eoff', kindsOff.join(','));
    for (const [facet] of FACETS) {
      const off = data.enums[facet].filter((v) => v && !filter.facets[facet].has(v));
      if (off.length) p.set(`${facet}off`, off.join(','));
    }
    if (filter.minConfidence) p.set('c', String(filter.minConfidence));
    if (filter.scope) { p.set('s', data.ids[filter.scope.index]); p.set('h', String(filter.scope.hops)); }
    if (filter.explore) {
      p.set('open', [...filter.explore.open].map((i) => data.ids[i]).join(','));
      p.set('x', [...filter.explore.expanded].map((i) => data.ids[i]).join(','));
    }
    if (state.selection) p.set('sel', state.selection);
    history.replaceState(null, '', `#${p.toString()}`);
  }
  function readHash() {
    const p = new URLSearchParams(location.hash.slice(1));
    const preset = p.get('p') ?? 'spine';
    if (preset === 'custom' && p.get('t')) {
      for (const t of p.get('t').split(',')) filter.nodeTypes.add(t);
      state.preset = 'custom';
    }
    else {
      const key = PRESETS[preset] ? preset : 'spine';
      for (const t of PRESETS[key].types ?? data.nodeTypes) filter.nodeTypes.add(t);
      state.preset = key;
    }
    for (const k of (p.get('eoff') ?? '').split(',').filter(Boolean)) filter.edgeKinds.delete(k);
    for (const [facet] of FACETS)
      for (const v of (p.get(`${facet}off`) ?? '').split(',').filter(Boolean)) filter.facets[facet].delete(v);
    if (p.get('c')) filter.minConfidence = Number(p.get('c'));
    if (p.get('s') && index(p.get('s')) >= 0) filter.scope = {index: index(p.get('s')), hops: Number(p.get('h') ?? 1)};
    const ids = (key) => (p.get(key) ?? '').split(',').map(index).filter((i) => i >= 0);
    if (p.get('open')) filter.explore = {open: new Set(ids('open')), expanded: new Set(ids('x'))};
    if (p.get('sel') && index(p.get('sel')) >= 0) state.selection = p.get('sel');
    syncFacets();
  }

  window.kgv = {data, filter, graph, panel, questions, typeTree, kindTree, refresh: (o) => refresh(o), explore, isolate, popView};
  collapsibleSections();
  resizablePanes();
  readHash();
  syncPresets();
  await graph.ready();
  refresh({resimulate: true});
  if (state.selection) panel.node(state.selection);
  $('loading').hidden = true;
}

function counts(names, codes, only) {
  const c = new Uint32Array(names.length);
  for (let i = 0; i < codes.length; i++)
    if (!only || only[i]) c[codes[i]]++;
  return new Map(names.map((n, i) => [n, c[i]]));
}

function makeStyle(data) {
  const byName = new Map(data.schema.nodeTypes.map((t) => [t.name, t]));
  const rootOf = (name) => byName.get(name)?.root ?? 'type';
  const roots = data.schema.roots;
  const perRoot = new Map();
  const nodeColor = data.nodeTypes.map((name) => {
    const root = rootOf(name);
    const k = perRoot.get(root) ?? 0;
    perRoot.set(root, k + 1);
    return shade(ROOT_COLORS[root] ?? '#888888', ((k % 5) - 2) * 0.08);
  });
  const nodeCss = data.nodeTypes.map((name, i) => rgbCss(nodeColor[i]));
  const cluster = data.nodeTypes.map((name) => roots.indexOf(rootOf(name)));
  const edgeColor = data.edgeKinds.map((name, i) => hslToRgb((i * 137.508) % 360, 0.45, 0.5));
  const edgeCss = data.edgeKinds.map((name, i) => rgbCss(edgeColor[i]));
  return {nodeColor, nodeCss, cluster, edgeColor, edgeCss, rootOf, size: (degree) => Math.min(10, 2.2 + Math.log2(1 + degree) * 0.8),
    typeCss: (type) => nodeCss[data.nodeTypes.indexOf(type)] ?? ROOT_COLORS[rootOf(type)] ?? '#888'};
}

/** Roots on a ring, each type's nodes scattered around its root's spot, so the simulation starts sorted. */
function seedPositions(data, style) {
  const pos = new Float32Array(data.nodes * 2);
  const roots = data.schema.roots.length;
  for (let i = 0; i < data.nodes; i++) {
    const r = style.cluster[data.type[i]];
    const angle = (r / roots) * Math.PI * 2;
    const h1 = hash(i * 2 + 1), h2 = hash(i * 2 + 2);
    const radius = 8192 * 0.16 * Math.sqrt(h1);
    const theta = h2 * Math.PI * 2;
    pos[i * 2] = 4096 + Math.cos(angle) * 8192 * 0.28 + Math.cos(theta) * radius;
    pos[i * 2 + 1] = 4096 + Math.sin(angle) * 8192 * 0.28 + Math.sin(theta) * radius;
  }
  return pos;
}

function hash(n) {
  const x = Math.sin(n * 12.9898 + 78.233) * 43758.5453;
  return x - Math.floor(x);
}

function nodeTypeItems(data, style) {
  const present = new Set(data.nodeTypes);
  return data.schema.nodeTypes.filter((t) => t.name !== 'node').map((t) => ({
    name: t.name, parent: t.extends === 'node' ? undefined : t.extends, abstract: t.abstract, description: t.description,
    color: t.abstract ? ROOT_COLORS[t.root ?? t.name] ?? '#888' : style.typeCss(t.name), line: false, empty: !t.abstract && !present.has(t.name)}));
}

const EDGE_GROUPS = {tree: 'the hierarchy and the links along it', ownership: 'what code a feature owns', code: 'what the AST sees between components',
  evidence: 'what vouches for a feature', work: 'tickets, releases and what they carry', people: 'who asked, who belongs, who is served',
  infra: 'pipelines, environments and hosts'};

/** Edge kinds under their group folder (conventions §4), the reference properties under `ref`. */
function edgeKindItems(data, style) {
  const items = [];
  const groups = data.schema.edgeGroups ?? [...new Set(data.schema.edgeTypes.map((e) => e.group).filter(Boolean))];
  for (const g of groups) items.push({name: g, parent: undefined, abstract: true, description: EDGE_GROUPS[g] ?? g, color: '#888', line: true});
  for (const e of data.schema.edgeTypes) {
    if (e.abstract) continue;
    items.push({name: e.name, parent: e.group || e.extends, abstract: false, description: e.description,
      color: style.edgeCss[data.edgeKinds.indexOf(e.name)] ?? '#bbb', line: true, empty: !data.edgeKinds.includes(e.name)});
  }
  items.push({name: 'ref', parent: undefined, abstract: true, description: 'reference properties, one rel table each', color: '#888', line: true});
  for (const r of data.schema.refs)
    items.push({name: r.name, parent: 'ref', abstract: false, description: `${r.declaredBy.join(', ')} → ${r.to.join(' | ')}`, color: style.edgeCss[data.edgeKinds.indexOf(r.name)] ?? '#bbb', line: true,
      empty: !data.edgeKinds.includes(r.name)});
  return items;
}

function badge(type, style) {
  const b = el('span', {class: 'kgv-badge'}, type ?? '?');
  b.style.background = style.typeCss(type);
  return b;
}

function shade(hex, amount) {
  const [r, g, b] = [1, 3, 5].map((i) => parseInt(hex.slice(i, i + 2), 16) / 255);
  return [r, g, b].map((c) => Math.max(0, Math.min(1, amount >= 0 ? c + (1 - c) * amount : c * (1 + amount))));
}

function hslToRgb(h, s, l) {
  const c = (1 - Math.abs(2 * l - 1)) * s, x = c * (1 - Math.abs((h / 60) % 2 - 1)), m = l - c / 2;
  const [r, g, b] = h < 60 ? [c, x, 0] : h < 120 ? [x, c, 0] : h < 180 ? [0, c, x] : h < 240 ? [0, x, c] : h < 300 ? [x, 0, c] : [c, 0, x];
  return [r + m, g + m, b + m];
}

function rgbCss([r, g, b]) {
  return `rgb(${Math.round(r * 255)}, ${Math.round(g * 255)}, ${Math.round(b * 255)})`;
}

/** Every left-pane section folds on its heading; the folded set is remembered per browser. */
function collapsibleSections() {
  let folded;
  try {
    folded = new Set(JSON.parse(localStorage.getItem('kgv-folded') ?? '[]'));
  }
  catch {
    folded = new Set();
  }
  for (const section of document.querySelectorAll('.kgv-left .kgv-section')) {
    const h3 = section.querySelector('h3');
    if (!h3) continue;
    const key = h3.firstChild.textContent.trim().toLowerCase();
    section.classList.toggle('kgv-collapsed', folded.has(key));
    h3.addEventListener('click', () => {
      const collapsed = section.classList.toggle('kgv-collapsed');
      if (collapsed) folded.add(key);
      else folded.delete(key);
      try {
        localStorage.setItem('kgv-folded', JSON.stringify([...folded]));
      }
      catch {
        // no storage: the fold still works for this page
      }
    });
  }
}

function resizablePanes() {
  const main = $('main');
  for (const gutter of document.querySelectorAll('.kgv-gutter')) {
    const side = gutter.dataset.side;
    gutter.addEventListener('mousedown', (ev) => {
      ev.preventDefault();
      const move = (e) => {
        const width = side === 'left' ? e.clientX : window.innerWidth - e.clientX;
        main.classList.toggle(`kgv-${side}-closed`, width < 60);
        if (width >= 60) main.style.setProperty(`--${side}`, `${Math.min(Math.max(width, 180), window.innerWidth / 2)}px`);
      };
      const up = () => { window.removeEventListener('mousemove', move); window.removeEventListener('mouseup', up); };
      window.addEventListener('mousemove', move);
      window.addEventListener('mouseup', up);
    });
    gutter.addEventListener('dblclick', () => main.classList.toggle(`kgv-${side}-closed`));
  }
}

main().catch((e) => {
  const box = $('loading');
  box.classList.add('kgv-error');
  box.textContent = `${e.message}\n\nIs grok kg serve running over a built generation?`;
  console.error(e);
});
