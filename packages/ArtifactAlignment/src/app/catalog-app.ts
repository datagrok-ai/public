import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SCHEMA, T_ALIGNMENT, T_COMPOUND, T_PROGRAM, T_PROGRAM_COMPOUND, T_STUDY} from '../domain/constants';
import {approvePublication, getPublicationHistory} from '../service/publication-service';
import {showCurateDialog, showProgramDialog, showRejectDialog} from './admin-dialogs';
import {createGuardedDomainView} from './cold-boot-workarounds';

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

function alignmentView(permanentFilter: string, name: string): Promise<DG.ViewBase> {
  return createGuardedDomainView({schema: SCHEMA, table: 'alignment', permanentFilter, name, showFilters: true});
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

function registryView(table: string, name: string): Promise<DG.ViewBase> {
  return createGuardedDomainView({schema: SCHEMA, table, name});
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
