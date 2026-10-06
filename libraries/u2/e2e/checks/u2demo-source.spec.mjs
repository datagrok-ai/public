/* The current sub-demo's source in the context panel (packages/U2Demo/src/source-panel.ts):
   navigation pushes a DemoPage through grok.shell.o, and the handler renders the file text that
   webpack inlined into this build (`?raw` → asset/source), sliced to the page factory. Nothing is
   fetched — GitHub is a link in the header, not the transport — so both lanes assert the same
   ready state, and a proxy request during navigation is itself a failure. */
import {ok, shot} from '../local.mjs';
import {openDemoPage} from '../lib.mjs';

const PANEL = '.u2demo-source';
const CAPTION = '.u2demo-source-caption';

let proxyRequests = 0;
let counting = false;

const state = (page) => page.evaluate((sel) => {
  const panel = document.querySelector(sel);
  if (!panel)
    return {present: false, caption: '', about: '', github: '', text: 0, head: '', missing: false};
  const pre = panel.querySelector('.u2demo-source-text');
  return {
    present: true,
    caption: panel.querySelector('.u2demo-source-caption')?.textContent ?? '',
    about: panel.querySelector('.u2demo-source-about')?.textContent ?? '',
    github: panel.querySelector('.u2demo-source-header a[href*="github.com"]')?.getAttribute('href') ?? '',
    missing: panel.querySelector('.u2demo-source-missing') != null,
    text: pre?.textContent.length ?? 0,
    head: (pre?.textContent ?? '').slice(0, 60),
  };
}, PANEL);

const clickLeaf = (page, label) => page.evaluate((label) => {
  const row = [...document.querySelectorAll('.u2demo-nav .u2-tree-row')]
    .find((r) => r.querySelector('.u2-tree-label')?.textContent.trim() === label && r.offsetParent != null);
  if (row == null)
    return false;
  row.click();
  return true;
}, label);

const waitCaption = (page, text) => page.waitForFunction(({sel, text}) =>
  document.querySelector(sel)?.textContent === text,
{sel: CAPTION, text}, {timeout: 10000}).then(() => true).catch(() => false);

export async function fixture(page) {
  if (!counting) {
    counting = true;
    page.on('request', (r) => {
      if (r.url().includes('/connectors/proxy'))
        proxyRequests++;
    });
  }
  await openDemoPage(page, 'lists');
  await page.evaluate(() => {
    grok.shell.windows.showContextPanel = true;
  });
}

async function checkCaption(page) {
  const clicked = await clickLeaf(page, 'Trees');
  const landed = await waitCaption(page, 'Collections / Trees');
  const s = await state(page);
  await shot(page, 'u2demo-source-1-caption');
  ok('u2demo-source/1a/navigating-to-trees-captions-the-panel-and-describes-what-the-page-covers',
    clicked && landed && s.about.length > 0 &&
    s.github === 'https://github.com/datagrok-ai/public/blob/master/packages/U2Demo/src/pages/collections.ts',
    `clicked=${clicked} caption="${s.caption}" about="${s.about}" github="${s.github}"`);
}

/** The whole point of bundling: the text is the running build's, available with no round trip, so
 * the panel is never a 404 and never a spinner — in either lane. */
async function checkBundled(page) {
  const before = proxyRequests;
  await clickLeaf(page, 'Lists');
  await waitCaption(page, 'Collections / Lists');
  const lists = await state(page);
  await clickLeaf(page, 'Trees');
  await waitCaption(page, 'Collections / Trees');
  const trees = await state(page);
  await shot(page, 'u2demo-source-2-bundled');
  // the head proves the slice landed on the factory: an LF-only close-brace search returns the
  // whole file from a CRLF checkout, which `text > 0` happily accepts
  ok('u2demo-source/2a/the-panel-shows-the-bundled-factory-source-without-fetching-anything',
    /^(export )?(async )?function listsPage\b/.test(lists.head) &&
    /^(export )?(async )?function treesPage\b/.test(trees.head) &&
    !lists.missing && !trees.missing && trees.text > 200 && proxyRequests === before,
    `lists=${JSON.stringify(lists.head.slice(0, 32))} trees=${JSON.stringify(trees.head.slice(0, 32))} ` +
    `missing=${lists.missing || trees.missing} chars=${trees.text} requests ${before}->${proxyRequests}`);
}

/** M7: `white-space: pre` put the horizontal scrollbar at the foot of an 1885 px `<pre>`, so a
 * third of every long line was unreachable. The text wraps now — nothing scrolls sideways. */
async function checkWrapping(page) {
  const box = await page.evaluate((sel) => {
    const pre = document.querySelector(`${sel} .u2demo-source-text`);
    if (pre == null)
      return null;
    return {white: getComputedStyle(pre).whiteSpace, scroll: pre.scrollWidth, client: pre.clientWidth};
  }, PANEL);
  await shot(page, 'u2demo-source-3-wrapping');
  ok('u2demo-source/3a/the-source-text-wraps-instead-of-scrolling-past-the-panel-width',
    box != null && box.white === 'pre-wrap' && box.scroll <= box.client + 1,
    box == null ? 'no source text' : `white-space=${box.white} scrollWidth=${box.scroll} clientWidth=${box.client}`);
}

async function checkSlice(page) {
  const result = await page.evaluate(async () => {
    const sibling = DG.Func.find({package: 'U2demo'})[0];
    if (sibling == null)
      return {error: 'U2demo is not in the fixture'};
    await sibling.package.load();
    const mod = sibling.package.getModule('package.js');
    if (typeof mod?.sliceSymbol !== 'function')
      return {error: 'sliceSymbol is not exported from the package module'};
    const canned = 'import x;\n\nexport function other(): void {\n  return;\n}\n\n' +
      'function asyncPage(): HTMLElement {\n  const el = div();\n  return el;\n}\n\nconst tail = 1;\n';
    return {
      sliced: mod.sliceSymbol(canned, 'asyncPage'),
      // the bundled text comes from the checkout, which is CRLF; only GitHub raw is LF
      crlf: mod.sliceSymbol(canned.replace(/\n/g, '\r\n'), 'asyncPage'),
      whole: mod.sliceSymbol(canned) === canned,
      missing: mod.sliceSymbol(canned, 'nowhere') === canned,
    };
  });
  ok('u2demo-source/4a/sliceSymbol-cuts-the-factory-out-of-a-canned-file-and-falls-back-whole',
    result.error == null && result.sliced.startsWith('function asyncPage(') &&
    result.sliced.endsWith('}\n') && result.sliced.includes('return el;') &&
    result.crlf.startsWith('function asyncPage(') && result.crlf.endsWith('}\r\n') &&
    !result.crlf.includes('const tail') &&
    result.whole && result.missing,
    result.error ?? `sliced=${JSON.stringify(result.sliced?.slice(0, 40))} ` +
      `crlf=${JSON.stringify(result.crlf?.slice(-20))} whole=${result.whole} missing=${result.missing}`);
}

export const checks = [
  {id: 'u2demo-source/1 navigation captions and describes the sub-demo in the context panel', run: checkCaption},
  {id: 'u2demo-source/2 the source is bundled into the build, not fetched', run: checkBundled},
  {id: 'u2demo-source/3 the source text wraps to the panel width', run: checkWrapping},
  {id: 'u2demo-source/4 sliceSymbol slices a top-level function, whole file otherwise', run: checkSlice},
];
