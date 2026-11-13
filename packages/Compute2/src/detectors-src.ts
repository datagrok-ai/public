/* eslint-disable */
// TypeScript source of ../detectors.js — compiled by `tsc -p tsconfig.detectors.json` (plain
// script output, no bundler: the platform evaluates detectors.js as a global script and scans
// its annotation comments, so the emitted file must stay module-free with comments intact).
declare const DG: typeof import('datagrok-api/dg');
declare const grok: typeof import('datagrok-api/grok');
declare const ui: typeof import('datagrok-api/ui');
declare const dayjs: typeof import('dayjs');

/**
 * Loaded by the platform at startup without the main package bundle. Besides semantic type
 * detectors (none yet), it hosts the autostart that registers the ComputeRun object handler,
 * so saved Compute2 runs (e.g. linked into a space) render and open through their own editor
 * even before any Compute2 function has run. The heavy work (preview, open) goes through
 * DG.Func.byName('Compute2:...'), which loads the main bundle on demand.
 */
class Compute2PackageDetectors extends DG.Package {
  //name: autostart
  //meta.role: autostart
  //meta.autostartImmediate: true
  //description: Registers the Compute2 run handler at platform startup
  autostart() {
    if (DG.ObjectHandler.list().some((h) => h.type === 'ComputeRun'))
      return;
    DG.ObjectHandler.register(new ComputeRunHandler());
  }
}

function openComputeRun(id: string): Promise<any> {
  return DG.Func.byName('Compute2:OpenWorkflowRun').prepare({id: id}).call();
}

function computeRunEditorOf(call: any): string | null {
  const editor = call.func && call.func.options ? call.func.options['editor'] : null;
  return typeof editor === 'string' && editor.startsWith('Compute2:') ? editor : null;
}

/** Same contract as the ModelHandler one entity level up: makes saved Compute2 runs render
 * and open through their editor instead of the platform's generic funccall view. */
class ComputeRunHandler extends DG.ObjectHandler {
  get type(): string {
    return 'ComputeRun';
  }

  isApplicable(x: any): boolean {
    try {
      const call = DG.toJs(x);
      return call instanceof DG.FuncCall && computeRunEditorOf(call) != null;
    } catch (_) {
      return false;
    }
  }

  async getById(id: string): Promise<any> {
    return grok.dapi.functions.calls.allPackageVersions().include('session.user,func.package').find(id);
  }

  async refresh(x: any): Promise<any> {
    return this.getById(DG.toJs(x).id);
  }

  // The run name: the legacy title when set, otherwise the friendly name, else the model name.
  getCaption(x: any): string {
    const call = DG.toJs(x);
    return call.options['title'] || call.friendlyName || this.modelName(call);
  }

  modelName(call: any): string {
    return (call.func && (call.func.friendlyName || call.func.name)) || 'Run';
  }

  // Copies of the run browser rules (compute-utils getRunCreated, Compute2 useRunDetails and runVersion):
  // runs saved before the createdById option fall back to their start time and platform author.
  savedText(call: any): string {
    try {
      const saved = call.options['createdById'] ? dayjs(call.options['createdAt']) : call.started;
      return saved && saved.isValid() ? saved.format('YYYY-MM-DD HH:mm') : '';
    } catch (_) {
      return '';
    }
  }

  // after someone else edited the details, the platform author is the last editor
  async authors(run: any): Promise<{by: string, modifiedBy?: string}> {
    const authorId = run.options['createdById'];
    const lastEditor = run.author;
    if (!authorId || authorId === (lastEditor && lastEditor.id))
      return {by: (lastEditor && lastEditor.friendlyName) || ''};
    const author = await grok.dapi.users.find(authorId).then((user: any) => (user && user.friendlyName) || '', () => '');
    return {by: author, modifiedBy: lastEditor && lastEditor.friendlyName};
  }

  renderIcon(_x: any): HTMLElement {
    return ui.iconFA('file-chart-line');
  }

  renderTooltip(x: any): HTMLElement {
    const call = DG.toJs(x);
    return ui.divV([
      ui.divText(this.getCaption(call)),
      ui.divText([this.modelName(call), this.savedText(call)].filter((s) => s).join(' · ')),
    ]);
  }

  renderProperties(x: any): HTMLElement {
    const call = DG.toJs(x);
    return ui.wait(async () => {
      const run = await this.getById(call.id);
      const {by, modifiedBy} = await this.authors(run);
      const map: {[key: string]: string} = {'Model': this.modelName(run), 'Saved': this.savedText(run), 'By': by};
      if (modifiedBy)
        map['Modified by'] = modifiedBy;
      const tags = ((run.options['tags'] || []) as string[]).join(', ');
      if (tags)
        map['Tags'] = tags;
      const version = run.options['version'] || (run.func && run.func.package && run.func.package.version);
      if (version)
        map['Version'] = version;
      const description = run.options['description'];
      return ui.divV([
        ui.h2(this.getCaption(run)),
        ...description ? [ui.divText(description)] : [],
        ui.tableFromMap(map),
        ui.div([ui.bigButton('Open run', () => void openComputeRun(run.id))], {style: {marginTop: '8px'}}),
      ]);
    });
  }

  async renderPreview(x: any): Promise<any> {
    const call = DG.toJs(x);
    const editorCall = call.func.prepare({});
    editorCall.aux.initialRunId = call.id;
    // Return the editor view directly so the docked preview is the view RFV/TreeWizard controls
    // (wrapping it makes a separate view the tracked preview — see compute-utils ModelHandler).
    const res = await DG.Func.byName(computeRunEditorOf(call)!).prepare({call: editorCall}).call();
    return res.getOutputParamValue();
  }

  init(): void {
    this.registerParamFunc('Open run', (x: any) => void openComputeRun(DG.toJs(x).id));
  }
}
