/* The boot half of the js-api bundle gate (the eval half is js-api/scripts/smoke-bundle.cjs):
   boots the real client in local mode and requires that the platform globals came up and the
   typed js-api surface is live. Catches what a stubbed eval cannot — a bundle that evaluates
   but breaks the client. Run after regenerating js-api.js; needs `grok-core local up`. */
import {consoleErrors, ok, openLocal, report} from './local.mjs';

const {browser, page} = await openLocal();

const probe = await page.evaluate(() => ({
  grok: typeof window.grok,
  classes: Object.fromEntries(['Viewer', 'Grid', 'FormViewer', 'Point', 'DataFrame', 'Column']
    .map((n) => [n, typeof window.DG?.[n]])),
  registerCleanup: typeof window.DG?.Widget?.registerCleanup,
  u2core: typeof window.DG?.U2?.Control,
  widgetIsControl: typeof window.DG?.U2?.Control === 'function' &&
    window.DG.Widget.prototype instanceof window.DG.U2.Control,
  user: window.grok?.shell?.user?.login ?? null,
}));

ok('boot/globals', probe.grok === 'object' && probe.user !== null,
  `grok=${probe.grok} user=${probe.user}`);
ok('boot/core-classes', Object.values(probe.classes).every((t) => t === 'function'),
  JSON.stringify(probe.classes));
ok('boot/widget-registerCleanup', probe.registerCleanup === 'function',
  `typeof=${probe.registerCleanup}`);
ok('boot/u2core', probe.u2core === 'function', `typeof=${probe.u2core}`);
ok('boot/widget-is-control', probe.widgetIsControl === true,
  `widgetIsControl=${probe.widgetIsControl}`);
ok('boot/console-clean', consoleErrors.length === 0, consoleErrors.join(' | '));

await browser.close();
process.exit(report('boot-smoke'));
