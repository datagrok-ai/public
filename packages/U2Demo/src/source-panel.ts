/* The current sub-demo's source in the context panel — a `DemoPage` holder pushed through the
   current-object channel (the PropRowHandler shape), the text bundled into this build by webpack
   (`?raw` → `asset/source`) and sliced to the page factory. So the panel shows the source that is
   actually running, and GitHub stays a link rather than the transport. Plain text for now; a code
   editor comes later. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {divV, link, span} from '@datagrok-libraries/u2';
import type {DemoLeaf} from './nav';
import {SRC_ROOT, setPageShownHook} from './nav';
import convergenceSource from './convergence?raw';
import funcConvergenceSource from './func-convergence?raw';
import funcsSource from './funcs?raw';
import funcHistorySource from './func-history?raw';
import collectionsSource from './pages/collections?raw';
import containersSource from './pages/containers?raw';
import displaySource from './pages/display?raw';
import formsSource from './pages/forms?raw';
import inputsSource from './pages/inputs?raw';
import overviewSource from './pages/overview?raw';
import platformSource from './pages/platform?raw';
import msaWorkbenchSource from './pages/msa-workbench?raw';
import '../css/u2demo-source.css';

const BLOB = 'https://github.com/datagrok-ai/public/blob/master/';

/** Every file a `SourceRef` can name, keyed the way `nav.ts` names it. A leaf pointing anywhere
 * else renders the "not bundled" state — the one place this map can drift. */
const SOURCES: Record<string, string> = {
  [`${SRC_ROOT}/convergence.ts`]: convergenceSource,
  [`${SRC_ROOT}/func-convergence.ts`]: funcConvergenceSource,
  [`${SRC_ROOT}/funcs.ts`]: funcsSource,
  [`${SRC_ROOT}/func-history.ts`]: funcHistorySource,
  [`${SRC_ROOT}/pages/collections.ts`]: collectionsSource,
  [`${SRC_ROOT}/pages/containers.ts`]: containersSource,
  [`${SRC_ROOT}/pages/display.ts`]: displaySource,
  [`${SRC_ROOT}/pages/forms.ts`]: formsSource,
  [`${SRC_ROOT}/pages/inputs.ts`]: inputsSource,
  [`${SRC_ROOT}/pages/overview.ts`]: overviewSource,
  [`${SRC_ROOT}/pages/platform.ts`]: platformSource,
  [`${SRC_ROOT}/pages/msa-workbench.ts`]: msaWorkbenchSource,
};

/** From `function <symbol>(` to the first line that is exactly `}` — every page factory in this
 * package is a top-level function closed at column 0. Both line endings: the checkout is CRLF,
 * GitHub raw serves LF. */
export function sliceSymbol(text: string, symbol?: string): string {
  if (!symbol)
    return text;
  const start = text.search(new RegExp(`^(export )?(async )?function ${symbol}\\b`, 'm'));
  if (start < 0)
    return text;
  const closing = /\r?\n\}\r?\n/g;
  closing.lastIndex = start;
  const end = closing.exec(text);
  return end == null ? text.slice(start) : text.slice(start, end.index + end[0].length);
}

/** What navigation makes current: the sub-demo whose source the panel shows. */
export class DemoPage {
  constructor(readonly leaf: DemoLeaf) {}
}

function sourceText(text: string, symbol?: string): HTMLElement {
  const pre = document.createElement('pre');
  pre.className = 'u2demo-source-text';
  pre.textContent = sliceSymbol(text, symbol);
  return pre;
}

class DemoSourceHandler extends DG.ObjectHandler<DemoPage> {
  get type(): string {
    return 'u2demo-source';
  }

  isApplicable(x: any): boolean {
    return x instanceof DemoPage;
  }

  getCaption(x: DemoPage): string {
    return `${x.leaf.group.label} / ${x.leaf.label}`;
  }

  renderProperties(x: DemoPage): HTMLElement {
    const ref = x.leaf.source;
    const text = SOURCES[ref.file];
    return divV([
      divV([
        span(this.getCaption(x), 'u2demo-source-caption'),
        span(x.leaf.description, 'u2demo-source-about'),
        link(ref.file, BLOB + ref.file),
        span('The source of the build you are running; the link opens master, which may differ.',
          'u2demo-source-note'),
      ], 'u2demo-source-header'),
      text == null ? span(`Not bundled: ${ref.file}`, 'u2demo-source-missing') :
        sourceText(text, ref.symbol),
    ], 'u2demo-source');
  }
}

let registered = false;

/** Registered once per session (the platform keeps every handler it is given), and wires
 * navigation: every page shown becomes the panel's current object. */
export function registerDemoSourceHandler(): void {
  if (registered)
    return;
  registered = true;
  DG.ObjectHandler.register(new DemoSourceHandler());
  setPageShownHook((leaf: DemoLeaf, open = false) => {
    if (open)
      grok.shell.windows.showContextPanel = true;
    // forced: `grok.shell.o` mutes the channel for a second and drops muted writes, so two nav
    // clicks under 1s apart would leave the panel on the previous sub-demo's source for good
    grok.shell.setCurrentObject(new DemoPage(leaf), false, true);
  });
}
