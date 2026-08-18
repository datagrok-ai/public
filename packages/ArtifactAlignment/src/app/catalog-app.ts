import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {SCHEMA, T_ALIGNMENT, T_PROGRAM} from '../domain/constants';
import {getPublicationHistory} from '../service/publication-service';

/** Opens the frozen (readonly) workflow run in Compute2; the serialized config of a
 * frozen copy carries isReadonly on every node, so the whole tree opens readonly. */
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

async function gridMaterialized(view: DG.ViewBase, timeoutMs: number = 3000): Promise<boolean> {
  const started = Date.now();
  while (Date.now() - started < timeoutMs) {
    if (view.root.querySelector('canvas, .d4-item-card') != null)
      return true;
    await new Promise((resolve) => setTimeout(resolve, 150));
  }
  return false;
}

async function alignmentView(permanentFilter: string, name: string,
  attempt: number = 0): Promise<DG.ViewBase> {
  // On a cold start (deep link into the app) creating a DomainView while the client
  // is still booting intermittently dies with a NullError in refreshGrid: the view's
  // column set is built from entity-type/property-schema meta that arrives with the
  // boot sequence, and a not-yet-loaded column comes through as null. The /domains
  // routes are sequenced after that load; a JS app function is not, and no JS call
  // can await the Dart-side cache. Warm what is awaitable, then verify the grid
  // actually materialized and rebuild the view once the boot has settled.
  // Platform ask: an awaitable registry-ready signal for JS apps.
  await grok.dapi.domains.table(T_ALIGNMENT).capabilities();
  const view = DG.DomainView.create({schema: SCHEMA, table: 'alignment', permanentFilter});
  view.name = name;
  whenDocked(view, async () => {
    view.showFilters();
    if (!await gridMaterialized(view) && attempt < 4) {
      view.close();
      grok.shell.addView(await alignmentView(permanentFilter, name, attempt + 1) as DG.View);
    }
  });
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
