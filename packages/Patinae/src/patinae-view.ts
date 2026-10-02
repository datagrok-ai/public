import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {getPatinae} from './patinae-loader';
import {CommandMessage, PatinaeViewer} from './patinae-types';
import {baseName, runPml} from './pml-script';
import '../css/patinae-panels.css';
import '../css/patinae.css';

export interface PymolPreview {
  view: DG.View;
  /** Resolves once the file is loaded; rejects if WebGPU is unavailable or loading failed. */
  loaded: Promise<{viewer: PatinaeViewer, messages: CommandMessage[]}>;
}

const fileReader = {
  readText: (path: string) => grok.dapi.files.readAsText(path),
  readBytes: (path: string) => grok.dapi.files.readAsBytes(path),
};

// Builds a preview view for a PyMOL session (.pse), Patinae session (.prs), or script (.pml).
export function openPymol(file: DG.FileInfo): PymolPreview {
  const view = DG.View.create();
  view.name = file.fileName;

  if (!('gpu' in navigator)) {
    view.append(ui.divText('Patinae requires WebGPU. Use Chrome or Edge 120 or newer.', 'patinae-message'));
    const loaded = Promise.reject(new Error('WebGPU is not available'));
    loaded.catch(() => {});
    return {view, loaded};
  }

  const host = ui.div([], 'patinae-host');
  const top = ui.div([], 'patinae-slot-top');
  const right = ui.div([], 'patinae-slot-right');
  const bottom = ui.div([], 'patinae-slot-bottom');
  view.append(ui.div([top, ui.div([host, right], 'patinae-middle'), bottom], 'patinae-view'));

  let closed = false;
  let viewer: PatinaeViewer | null = null;
  const sub = grok.events.onViewRemoved.subscribe((v: DG.View) => {
    if (v.id !== view.id) return;
    closed = true;
    viewer?.destroy();
    sub.unsubscribe();
  });

  const loaded = (async () => {
    const {PatinaeViewer} = await getPatinae();
    viewer = new PatinaeViewer(host, {
      layout: [
        {name: 'repl', slot: 'top'}, {name: 'objects', slot: 'right'},
        {name: 'movie', slot: 'bottom'}, {name: 'sequence', slot: 'bottom'},
      ],
      slots: {top, right, bottom},
      picking: true,
    });
    if (closed) {
      viewer.destroy();
      throw new Error('Preview closed before the viewer initialized');
    }
    await viewer.init();

    const messages: CommandMessage[] = [];
    if (file.extension === 'pml') {
      const dir = file.fullPath.slice(0, file.fullPath.lastIndexOf('/'));
      const output = top.querySelector('.repl-output');
      const echo = (command: string, lines: CommandMessage[]) => {
        for (const [cls, text] of [['cmd', `Patinae> ${command}`], ...lines.map((m) => [m.level, m.text])])
          output?.appendChild(ui.divText(text, `repl-line,repl-${cls}`));
        output?.scrollTo(0, output.scrollHeight);
      };
      messages.push(...await runPml(await file.readAsString(), dir, fileReader, viewer, echo));
    } else
      viewer.loadData(await file.readAsBytes(), baseName(file.fileName), file.extension);
    const errors = messages.filter((m) => m.level === 'error');
    if (errors.length > 0)
      grok.shell.error(errors.map((m) => m.text).join('\n'));
    return {viewer, messages};
  })();
  loaded.catch((e) => {
    if (!closed) grok.shell.error(e instanceof Error ? e.message : String(e));
  });

  return {view, loaded};
}
