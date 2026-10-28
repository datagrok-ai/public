import * as grok from 'datagrok-api/grok';
import {View} from '../../../tests/platform-doubles.mjs';
import {errors, expect, flush, liveScopes, mount, need} from '../../lib/harness.mjs';

function current(name) {
  const view = name === null ? null : View.fromRoot(document.createElement('div'));
  if (view !== null)
    view.name = name;
  grok.shell.v = view;
  grok.events.onCurrentViewChanged.fire(view);
}

async function build(mod) {
  grok.shell.v = null;
  const base = {scopes: liveScopes(), subs: grok.events.onCurrentViewChanged.count};
  const badge = need(mod, 'currentViewBadge')();
  expect(typeof badge?.dispose === 'function', 'currentViewBadge() did not return a control');
  const root = mount(badge);
  await flush();
  return {badge, root, base};
}

export const checks = {
  'follows-current-view': async (mod) => {
    const {root} = await build(mod);
    expect(/no view/i.test(root.textContent), `with no view it shows "${root.textContent}"`);
    current('Campaigns');
    await flush();
    expect(root.textContent.includes('Campaigns'), `after switching to Campaigns it shows "${root.textContent}"`);
    current('KHL-3');
    await flush();
    expect(root.textContent.includes('KHL-3'), `after switching to KHL-3 it shows "${root.textContent}"`);
  },
  'stops-listening': async (mod) => {
    const {badge, base} = await build(mod);
    badge.dispose();
    await flush();
    const left = grok.events.onCurrentViewChanged.count - base.subs;
    expect(left === 0, `${left} subscriptions to onCurrentViewChanged outlive the disposed badge`);
    expect(liveScopes() === base.scopes, `${liveScopes() - base.scopes} scopes outlive the disposed badge`);
  },
  'quiet-after-dispose': async (mod) => {
    const {badge, root} = await build(mod);
    badge.dispose();
    const before = root.textContent;
    let reads = 0;
    const shell = grok.shell;
    Object.defineProperty(shell, 'v', {configurable: true, get() {
      reads++;
      return shell.dart.v;
    }, set(x) { shell.dart.v = x; }});
    try {
      current('Campaigns');
      await flush();
    }
    finally {
      delete shell.v;
    }
    expect(root.textContent === before && reads === 0, 'the disposed badge still reacts to view changes');
    expect(errors().length === 0, `logged: ${errors()[0]}`);
  },
};
