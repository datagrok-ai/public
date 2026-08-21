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

/** study id → protocol code. Render members are synchronous, so the map is warmed
 * lazily by the first card and cards re-render their facts once it arrives. */
let studyCodes: Map<string, string> | null = null;
let studiesLoading: Promise<Map<string, string>> | null = null;

function warmStudyCodes(): Promise<Map<string, string>> {
  return studiesLoading ??= grok.dapi.domains.table(T_STUDY).query({limit: 1000})
    .then((rows: any[]) => studyCodes = new Map(rows.map((r) => [r.id, r.protocol_code])))
    .catch((e) => {
      studiesLoading = null; // a failed warm must not disable study labels forever
      throw e;
    });
}

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

  override renderCard(x: DG.DomainRow, context: any = null): HTMLElement {
    const row = this.rowOf(x);
    if (row == null)
      return super.renderCard(x, context);
    const values = row.values ?? {};
    const detail: {[label: string]: string | undefined} = {};
    const detailHost = ui.div([], 'aa-card-details');
    const renderDetail = () => {
      ui.empty(detailHost);
      const filled = ['Status', 'Study', 'Workstream', 'Tags']
        .filter((label) => detail[label])
        .map((label) => [label, detail[label]!] as const);
      if (filled.length > 0)
        detailHost.appendChild(ui.tableFromMap(Object.fromEntries(filled)));
    };
    if (values['status'] != null && values['status'] !== 'approved')
      detail['Status'] = values['status'];
    detail['Workstream'] = values['workstream'];
    renderDetail();
    void warmStudyCodes().then(() => {
      detail['Study'] = studyCodes?.get(values['study_id']);
      renderDetail();
    }).catch(() => {});
    if (row.id != null) {
      grok.dapi.domains.table(T_ALIGNMENT)
        .query({filter: DG.cond('id', '=', row.id), expand: ['tags']})
        .then(([full]: any[]) => {
          detail['Tags'] = (full?.tags ?? []).map((t: any) => t.name).join(', ');
          renderDetail();
        }).catch(() => {});
    }
    // deliberately NOT ui.bind: the gallery wrapper already selects the row on
    // click, and a second binding turns one left click into select + open
    const name = values['name'] ?? row.displayName ?? '';
    return ui.divV([
      ui.divText(values['revision'] != null ? `${name} · v${values['revision']}` : name,
        {style: {fontWeight: 'bold'}}),
      detailHost,
    ], 'd4-gallery-item');
  }

  override renderProperties(x: DG.DomainRow, context: any = null): HTMLElement {
    const root = super.renderProperties(x, context);
    const values = this.rowOf(x)?.values ?? {};
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

/** URL base of the app's views when the package sets an absolute address itself
 * (tree navigation); views returned from the app function get this prefix from
 * the platform and carry a RELATIVE path instead — never both, that duplicates
 * the suffix. */
const APP_BASE = 'apps/ArtifactAlignment';

function alignmentView(permanentFilter: string, name: string): DG.ViewBase {
  const view = DG.DomainView.create({schema: SCHEMA, table: 'alignment', permanentFilter});
  view.name = name;
  return view;
}

export function catalogView(): DG.ViewBase {
  return alignmentView('status = "approved"', 'Artifact Catalog');
}

export function programView(programCode: string): DG.ViewBase {
  return alignmentView(`status = "approved" and program_id.code = "${programCode}"`, programCode);
}

export function reviewQueueView(): DG.ViewBase {
  return alignmentView('status != "approved"', 'Review queue');
}

export function myPublicationsView(): DG.ViewBase {
  return alignmentView('artifact_author = @current', 'My publications');
}

function registryView(table: string, name: string): DG.ViewBase {
  const view = DG.DomainView.create({schema: SCHEMA, table});
  view.name = name;
  return view;
}

/** Resolves the app's `path` URL input (everything after the app base) to the view
 * it addresses. The paths set here are relative — the platform prefixes the app
 * base for views returned from the app function. */
export function routeView(path?: string): DG.ViewBase {
  const seg = (path ?? '').split('?')[0].split('/').filter((s) => s !== '').map(decodeURIComponent);
  const routed = (view: DG.ViewBase, route: string): DG.ViewBase => {
    view.path = `/${route}`;
    return view;
  };
  if (seg[0] === 'programs' && seg[1] != null)
    return routed(programView(seg[1]), `programs/${encodeURIComponent(seg[1])}`);
  if (seg[0] === 'review-queue')
    return routed(reviewQueueView(), 'review-queue');
  if (seg[0] === 'my-publications')
    return routed(myPublicationsView(), 'my-publications');
  if (seg[0] === 'registry' && seg[1] === 'studies')
    return routed(registryView('study', 'Studies'), 'registry/studies');
  if (seg[0] === 'registry' && seg[1] === 'compounds')
    return routed(registryView('compound', 'Compounds'), 'registry/compounds');
  return catalogView();
}

/** Shows [view] as the Browse preview and publishes its address. The path is
 * absolute here: a tree-created view has no app call behind it, so nothing else
 * supplies the prefix. */
function preview(view: DG.ViewBase, route: string): void {
  grok.shell.preview = view as DG.View;
  view.path = `${APP_BASE}/${route}`;
}

export async function buildTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  treeNode.item('My publications').onSelected.subscribe(() =>
    preview(myPublicationsView(), 'my-publications'));
  treeNode.item('Review queue').onSelected.subscribe(() =>
    preview(reviewQueueView(), 'review-queue'));
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
    // keyed by the node's DOM root — the wrapper identity is not stable across events
    const nodePrograms = new Map<HTMLElement, any>();
    const editableIds = new Set<string>();
    for (const program of programs) {
      const item = programsNode.item(program.code);
      nodePrograms.set(item.root, program);
      // session-truth row permission gates the edit affordance; the save path
      // re-checks server-side, so this is presentation only
      void programHandler.rowFrom(program).permissions()
        .then((p) => p.edit && editableIds.add(program.id))
        .catch(() => {/* no edit affordance */});
      item.onSelected.subscribe(() => {
        preview(programView(program.code), `programs/${encodeURIComponent(program.code)}`);
        // the program's context panel renders the audience group columns as
        // clickable links — the member-management entry point for admins
        grok.shell.o = programHandler.rowFrom(program);
        grok.shell.windows.showContextPanel = true;
      });
    }
    // subscribed on the browse-tree root: context-menu events surface there;
    // the node map scopes the handler to this build's program items
    treeNode.rootNode.onNodeContextMenu.subscribe((event: any) => {
      const menu: DG.Menu | null = event?.args?.menu;
      const program = nodePrograms.get((event?.args?.item as DG.TreeViewNode)?.root);
      if (menu != null && program != null && editableIds.has(program.id))
        menu.item('Edit program…', () => showProgramDialog(program));
    });
  } catch (_) {/* no visible programs */}
  const registryNode = treeNode.group('Registry');
  registryNode.item('New study…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_STUDY));
  registryNode.item('New compound…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_COMPOUND));
  registryNode.item('Link compound to program…').onSelected.subscribe(() =>
    void DG.DomainObjectHandler.createRow(T_PROGRAM_COMPOUND));
  registryNode.item('Unlink compound from program…').onSelected.subscribe(async () => {
    const picked = await DG.DomainObjectHandler.pickRow(T_PROGRAM_COMPOUND);
    const values = picked?.values;
    if (values == null)
      return;
    const [program, compound] = await Promise.all([
      grok.dapi.domains.table(T_PROGRAM).get(values['program_id']),
      grok.dapi.domains.table(T_COMPOUND).get(values['compound_id']),
    ]);
    ui.dialog('Unlink compound')
      .add(ui.divText(`Remove ${compound?.registration_code ?? 'the compound'} from ` +
        `${program?.code ?? 'the program'}? Published artifacts keep their molecule links.`))
      .onOK(async () => {
        await grok.dapi.domains.table(T_PROGRAM_COMPOUND).delete(values['id']);
        grok.shell.info('Compound unlinked');
      })
      .show();
  });
  registryNode.item('Studies').onSelected.subscribe(() =>
    preview(registryView('study', 'Studies'), 'registry/studies'));
  registryNode.item('Compounds').onSelected.subscribe(() =>
    preview(registryView('compound', 'Compounds'), 'registry/compounds'));
}
