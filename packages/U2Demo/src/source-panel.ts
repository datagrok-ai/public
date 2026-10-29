/* The current sub-demo's source in the context panel — a `DemoPage` holder pushed through the
   current-object channel (the PropRowHandler shape), the text bundled into this build by rspack
   (`?raw` → `asset/source`) and sliced to the page factory. So the panel shows the source that is
   actually running, and GitHub stays a link rather than the transport. Plain text for now; a code
   editor comes later. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {divV, link, span} from '@datagrok-libraries/u2';
import type {DemoLeaf} from './nav';
import {setPageShownHook} from './nav';
import {SRC_ROOT} from './const';
import '../css/u2demo-source.css';

const BLOB = 'https://github.com/datagrok-ai/public/blob/master/';

/** Every module of this folder and of `pages/` as text (the `?raw` rule applies to each), keyed the way `nav.ts`
 * names a `SourceRef`: a page added there is bundled with nothing to list. */
const bundled = import.meta.webpackContext('./?raw', {recursive: true, regExp: /^\.\/(pages\/)?[\w-]+\.ts$/});
const SOURCES: Record<string, string> = Object.fromEntries(bundled.keys().map((key) =>
  [`${SRC_ROOT}/${key.slice(2)}`, bundled<string>(key)]));

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

/** A page's own code shown on the page: the bundled text of `file` (keyed as `nav.ts` names a `SourceRef`),
 * sliced to `symbol`. */
export function sourceOf(file: string, symbol?: string): HTMLElement {
  return SOURCES[file] == null ? span(`Not bundled: ${file}`, 'u2demo-source-missing') :
    sourceText(SOURCES[file], symbol);
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
