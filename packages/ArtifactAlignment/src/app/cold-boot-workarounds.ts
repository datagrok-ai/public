import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

/* Everything in this file works around ONE core defect and is deletable wholesale
 * once it is fixed: on a cold client boot (a deep link straight into the app) the
 * first Domain View can build its grid before the domain registry and per-table
 * entity metas have arrived, and dies with a NullError in
 * DataSourceCardView.refreshGrid / loadNextPage — a dead view that refresh()
 * cannot revive. Core ask: the cold-boot grid build awaits its own dependencies
 * (an awaitable registry-ready signal for JS apps).
 *
 * Strategy: PREVENT (await the metas the grid reads before creating the view),
 * then WATCH (grid content must materialize within a budget once docked), then
 * HEAL (placement-aware full rebuild — the preview is reassigned, an embedded
 * root is swapped in place). */

export interface GuardedDomainViewSpec {
  schema: string;
  table: string;
  permanentFilter?: string;
  name?: string;
  showFilters?: boolean;
}

/** Creates a Domain View that survives the cold-boot race — the only entry
 * point the app code needs; everything else here is plumbing. */
export async function createGuardedDomainView(spec: GuardedDomainViewSpec): Promise<DG.ViewBase> {
  await whenRegistryLoaded(`${spec.schema}.${spec.table}`);
  const view = build(spec);
  void healIfGridDied(view, spec);
  return view;
}

/** dartMeta is protected on the handler — this shim only exposes "has the
 * registry meta for the table arrived", the exact cache the grid is built from. */
class MetaProbe extends DG.DomainObjectHandler {
  get loaded(): boolean {
    return this.dartMeta != null;
  }
}

async function whenRegistryLoaded(table: string, timeoutMs: number = 10000): Promise<void> {
  const probe = new MetaProbe(table);
  const started = Date.now();
  while (!probe.loaded && Date.now() - started < timeoutMs)
    await delay(100);
  // The registry meta alone is not enough: the grid build also needs the
  // per-table entity metas, which the platform otherwise registers lazily
  // during the first Domain View creation — racing refreshGrid. Forcing the
  // registration up front is awaitable and idempotent.
  await (window as any).grok_DomainRowMeta_RegisterPerTableMetas?.();
}

function build(spec: GuardedDomainViewSpec): DG.DomainView {
  const view = DG.DomainView.create(
    {schema: spec.schema, table: spec.table, permanentFilter: spec.permanentFilter});
  if (spec.name != null)
    view.name = spec.name;
  if (spec.showFilters)
    whenDocked(view, () => view.showFilters());
  return view;
}

function whenDocked(view: DG.ViewBase, action: () => void, tries: number = 100): void {
  view.root.isConnected ? action() :
    tries > 0 ? void setTimeout(() => whenDocked(view, action, tries - 1), 100) : null;
}

/** Proof of life: the grid canvas (rendered even for an empty result) or a card. */
function materialized(view: DG.ViewBase): boolean {
  return view.root.querySelector('canvas, .d4-item-card') != null;
}

const delay = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

async function healIfGridDied(view: DG.DomainView, spec: GuardedDomainViewSpec,
  attempt: number = 0): Promise<void> {
  // the budget starts once the view is actually docked — cold boots dock late
  for (let waited = 0; !view.root.isConnected; waited += 200) {
    if (waited >= 10000)
      return;
    await delay(200);
  }
  for (let waited = 0; waited < 8000; waited += 300) {
    await delay(300);
    if (materialized(view))
      return;
  }
  if (attempt >= 2)
    return;
  const replacement = build(spec);
  void healIfGridDied(replacement, spec, attempt + 1);
  if (grok.shell.preview?.root === view.root) {
    // the dead view is the current Browse preview — a full, properly wired swap
    grok.shell.preview = replacement as DG.View;
    view.close();
  } else if (view.root.isConnected) {
    // docked somewhere the package does not control — swap the DOM in place
    replacement.root.style.width = '100%';
    replacement.root.style.height = '100%';
    ui.empty(view.root);
    view.root.appendChild(replacement.root);
  }
}
