import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SCHEMA, T_ALIGNMENT, T_PROGRAM} from '../domain/constants';
import {getPublicationHistory} from '../service/publication-service';

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

export class AlignmentHandler extends DG.DomainObjectHandler {
  constructor() {
    super(T_ALIGNMENT);
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
}

async function alignmentView(permanentFilter: string, name: string): Promise<DG.ViewBase> {
  await whenRegistryLoaded();
  const view = DG.DomainView.create({schema: SCHEMA, table: 'alignment', permanentFilter});
  view.name = name;
  whenDocked(view, () => view.showFilters());
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

export async function buildTreeBrowser(treeNode: DG.TreeViewGroup): Promise<void> {
  treeNode.item('My publications').onSelected.subscribe(async () =>
    grok.shell.preview = await myPublicationsView() as DG.View);
  treeNode.item('Review queue').onSelected.subscribe(async () =>
    grok.shell.preview = await reviewQueueView() as DG.View);
  const programsNode = treeNode.group('Programs');
  try {
    const programs = await grok.dapi.domains.table(T_PROGRAM).query({sort: 'code', limit: 100});
    for (const program of programs) {
      programsNode.item(program.code).onSelected.subscribe(async () =>
        grok.shell.preview = await programView(program.code) as DG.View);
    }
  } catch (_) {/* no visible programs */}
}
