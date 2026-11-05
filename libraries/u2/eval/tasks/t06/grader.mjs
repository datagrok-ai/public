import {errors, expect, flush, input, liveScopes, mount, need, type, until, wait} from '../../lib/harness.mjs';

const CAMPAIGNS = [{code: 'KHL-1', name: 'Kinase hits'}, {code: 'KHL-2', name: 'PI3K quick wins'}];

function deferred() {
  let resolve;
  const promise = new Promise((r) => resolve = r);
  return {promise, resolve};
}

async function build(mod) {
  const pending = deferred();
  const base = liveScopes();
  const panel = need(mod, 'campaignsPanel')(() => pending.promise);
  expect(typeof panel?.dispose === 'function', 'campaignsPanel() did not return a control');
  const root = mount(panel);
  await flush();
  return {panel, root, base, pending};
}

const rowsShown = (root) => root.querySelectorAll('.u2-input-root').length >= CAMPAIGNS.length;
const saysLoading = (root) => /loading/i.test(root.textContent) || [root, ...root.querySelectorAll('*')]
  .some((el) => el.getAttribute('aria-busy') === 'true' || el.getAttribute('role') === 'status');

async function loaded(mod) {
  const b = await build(mod);
  b.pending.resolve(CAMPAIGNS.map((c) => ({...c})));
  expect(await until(() => rowsShown(b.root)), 'the rows never show after the load resolves');
  return b;
}

export const checks = {
  'shows-rows': async (mod) => {
    const {root, pending} = await build(mod);
    expect(await until(() => saysLoading(root), 500), 'nothing shows that it is loading');
    pending.resolve(CAMPAIGNS.map((c) => ({...c})));
    expect(await until(() => rowsShown(root)), 'the rows never show after the load resolves');
    for (const c of CAMPAIGNS)
      expect(input(c.code, root).getState().text === c.name, `the ${c.code} row does not show "${c.name}"`);
  },
  'live-count': async (mod) => {
    const {root} = await loaded(mod);
    expect(/0 renamed/.test(root.textContent), 'the panel does not say "0 renamed"');
    await type(input('KHL-1', root), 'Kinase hits v2');
    expect(await until(() => /1 renamed/.test(root.textContent), 500),
      'after a rename the panel does not say "1 renamed"');
  },
  'dispose-releases': async (mod) => {
    const {panel, base} = await loaded(mod);
    panel.dispose();
    await flush();
    expect(liveScopes() === base, `${liveScopes() - base} scopes outlive the disposed panel`);
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
  'dispose-while-loading': async (mod) => {
    const {panel, root, base, pending} = await build(mod);
    panel.dispose();
    pending.resolve(CAMPAIGNS.map((c) => ({...c})));
    await wait(300);
    expect(liveScopes() === base, `${liveScopes() - base} scopes were built for a panel disposed while it loaded`);
    expect(root.querySelectorAll('.u2-input-root').length === 0,
      'rows were built into the panel after it was disposed');
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
};
