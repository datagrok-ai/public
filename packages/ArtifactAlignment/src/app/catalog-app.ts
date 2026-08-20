import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SCHEMA, T_ALIGNMENT, T_COMPOUND, T_PROGRAM, T_PROGRAM_COMPOUND, T_STUDY} from '../domain/constants';
import {approvePublication, getPublicationHistory} from '../service/publication-service';
import {showCurateDialog, showProgramDialog, showRejectDialog} from './admin-dialogs';

/** Opens the frozen workflow or function run in Compute2. The tree opens runnable —
 * a new run can be started from the artifact's state. */
export async function openArtifact(frozenMetaCallId: string): Promise<void> {
  try {
    const func = DG.Func.byName('Compute2:OpenWorkflowRun');
    await func.prepare({id: frozenMetaCallId}).call();
  } catch (_) {
    grok.shell.warning('Opening workflow runs requires the Compute2 package');
  }
}

/** Internal bookkeeping columns kept out of the catalog grid; all of them stay
 * reachable through the row's context panel (History pane, Open artifact). */
export const HIDDEN_GRID_COLUMNS = ['publication_id', 'artifact_id', 'source_artifact_id',
  'approved_by', 'approved_on', 'reject_reason'];

export class AlignmentHandler extends DG.DomainObjectHandler {
  constructor() {
    super(T_ALIGNMENT);
  }

  override renderGrid(grid: DG.Grid, options?: {items?: DG.DataFrame}): void {
    super.renderGrid(grid, options);
    for (const name of HIDDEN_GRID_COLUMNS) {
      const column = grid.columns.byName(name);
      if (column != null)
        column.visible = false;
    }
  }

  /** Double-click / Open on a catalog row opens the published artifact itself;
   * the row's entity view stays reachable via its deep link (Copy link). */
  override openRow(x: DG.DomainRow): void {
    const artifactId = x?.values?.['artifact_id'];
    if (artifactId != null)
      void openArtifact(artifactId);
    else
      super.openRow(x);
  }

  /** Whether the client's domain registry meta for this table has arrived —
   * dartMeta resolves from exactly the cache the Domain View's grid is built
   * from, and js-api re-resolves misses on every read. */
  get registryLoaded(): boolean {
    return this.dartMeta != null;
  }

  override renderProperties(x: DG.DomainRow, context: any = null): HTMLElement {
    const root = super.renderProperties(x, context);
    const values = x.values ?? {};
    const pane = ui.divV([], {style: {marginTop: '8px'}});
    pane.appendChild(ui.button('Open artifact', () => openArtifact(values['artifact_id'])));
    // The service re-checks approvers-group membership and author != approver, so the
    // buttons need no client-side gating for correctness.
    if (values['status'] === 'pending') {
      pane.appendChild(ui.divH([
        ui.bigButton('Approve', async () => {
          try {
            await approvePublication(values['id']);
            grok.shell.info('Publication approved — reopen the view to refresh');
          } catch (e: any) {
            grok.shell.error(`Approve failed: ${e?.message ?? e}`);
          }
        }),
        ui.button('Reject…', () => showRejectDialog(values['id'])),
      ], {style: {gap: '8px'}}));
    } else if (values['status'] === 'approved')
      pane.appendChild(ui.button('Curate…', () => void showCurateDialog(values['id'])));
    const historyHost = ui.divV([ui.divText('Loading history…', 'aa-history-loading')]);
    pane.appendChild(ui.h3('History'));
    pane.appendChild(historyHost);
    getPublicationHistory(values['publication_id']).then((entries) => {
      ui.empty(historyHost);
      if (entries.length === 0) {
        historyHost.appendChild(ui.divText('No archived versions'));
        return;
      }
      for (const entry of entries) {
        const reason = entry.final_status === 'rejected' && entry.reject_reason ?
          ` — "${entry.reject_reason}"` : '';
        historyHost.appendChild(ui.divH([
          ui.divText(`v${entry.revision} · ${entry.final_status}${reason}`),
          ui.link('open', () => openArtifact(entry.artifact_id), 'Open this archived version'),
        ], {style: {gap: '8px'}}));
      }
    }).catch(() => {
      ui.empty(historyHost);
      historyHost.appendChild(ui.divText('History unavailable'));
    });
    root.appendChild(pane);
    return root;
  }
}

function whenDocked(view: DG.ViewBase, action: () => void, tries: number = 100): void {
  view.root.isConnected ? action() :
    tries > 0 ? void setTimeout(() => whenDocked(view, action, tries - 1), 100) : null;
}

/** On a cold start (deep link into the app) the app function can run before the
 * client boot delivers the domain registry meta the Domain View's grid is built
 * from (NullError in refreshGrid otherwise). Wait on that meta directly before
 * creating the view. Platform ask stays: an awaitable registry-ready signal. */
async function whenRegistryLoaded(timeoutMs: number = 10000): Promise<void> {
  const probe = new AlignmentHandler();
  const started = Date.now();
  while (!probe.registryLoaded && Date.now() - started < timeoutMs)
    await new Promise((resolve) => setTimeout(resolve, 100));
  // The registry meta alone is not enough: the view's grid build also needs the
  // per-table entity metas, which the platform otherwise registers lazily during
  // the first Domain View creation — racing refreshGrid at cold boot (NullError
  // in refreshGrid/loadNextPage, dead view). Force the registration up front;
  // awaitable and idempotent.
  await (window as any).grok_DomainRowMeta_RegisterPerTableMetas?.();
}

/** Cold-boot core race: the first Domain View of a deep-linked session can die
 * while building its grid (NullError in DataSourceCardView.refreshGrid /
 * loadNextPage) — not preventable from a package (neither the registry meta nor
 * the per-table entity metas close the window) and not recoverable via
 * refresh(). Once the client is warm, creation is reliable — so when the grid
 * never materializes, transplant a freshly created view into the docked root.
 * Platform ask stays: make the cold-boot grid build await its dependencies. */
async function healIfGridDied(view: DG.DomainView, permanentFilter: string, attempt: number = 0): Promise<void> {
  for (let waited = 0; waited < 8000; waited += 500) {
    await new Promise((resolve) => setTimeout(resolve, 500));
    if (view.root.querySelector('canvas') != null)
      return;
  }
  if (!view.root.isConnected || attempt >= 2)
    return;
  const replacement = DG.DomainView.create({schema: SCHEMA, table: 'alignment', permanentFilter});
  replacement.root.style.width = '100%';
  replacement.root.style.height = '100%';
  ui.empty(view.root);
  view.root.appendChild(replacement.root);
  void healIfGridDied(replacement, permanentFilter, attempt + 1);
}

async function alignmentView(permanentFilter: string, name: string): Promise<DG.ViewBase> {
  await whenRegistryLoaded();
  const view = DG.DomainView.create({schema: SCHEMA, table: 'alignment', permanentFilter});
  view.name = name;
  whenDocked(view, () => view.showFilters());
  void healIfGridDied(view, permanentFilter);
  return view;
}

export function catalogView(): Promise<DG.ViewBase> {
  return alignmentView('status = "approved"', 'Artifact Catalog');
}

export function programView(programCode: string): Promise<DG.ViewBase> {
  return alignmentView(`status = "approved" and program_id.code = "${programCode}"`, programCode);
}

export function reviewQueueView(): Promise<DG.ViewBase> {
  return alignmentView('status != "approved"', 'Review queue');
}

export function myPublicationsView(): Promise<DG.ViewBase> {
  return alignmentView('artifact_author = @current', 'My publications');
}

async function registryView(table: string, name: string): Promise<DG.ViewBase> {
  await whenRegistryLoaded();
  const view = DG.DomainView.create({schema: SCHEMA, table});
  view.name = name;
  return view;
}

export async function buildTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  treeNode.item('My publications').onSelected.subscribe(async () =>
    grok.shell.preview = await myPublicationsView() as DG.View);
  treeNode.item('Review queue').onSelected.subscribe(async () =>
    grok.shell.preview = await reviewQueueView() as DG.View);
  const programsNode = treeNode.group('Programs');
  programsNode.item('New program…').onSelected.subscribe(() => showProgramDialog());
  programsNode.item('Edit program…').onSelected.subscribe(async () => {
    const picked = await DG.DomainObjectHandler.pickRow(T_PROGRAM);
    if (picked?.values != null)
      showProgramDialog(picked.values);
  });
  try {
    const programs = await grok.dapi.domains.table(T_PROGRAM).query({sort: 'code', limit: 100});
    const programHandler = new DG.DomainObjectHandler(T_PROGRAM);
    for (const program of programs) {
      programsNode.item(program.code).onSelected.subscribe(async () => {
        grok.shell.preview = await programView(program.code) as DG.View;
        // the program's context panel renders the audience group columns as
        // clickable links — the member-management entry point for admins
        grok.shell.o = programHandler.rowFrom(program);
        grok.shell.windows.showContextPanel = true;
      });
    }
  } catch (_) {/* no visible programs */}
  const registryNode = treeNode.group('Registry');
  registryNode.item('New study…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_STUDY));
  registryNode.item('New compound…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_COMPOUND));
  registryNode.item('Link compound to program…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_PROGRAM_COMPOUND));
  registryNode.item('Studies').onSelected.subscribe(async () =>
    grok.shell.preview = await registryView('study', 'Studies') as DG.View);
  registryNode.item('Compounds').onSelected.subscribe(async () =>
    grok.shell.preview = await registryView('compound', 'Compounds') as DG.View);
}
