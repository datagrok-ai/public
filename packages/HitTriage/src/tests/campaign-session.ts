import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import * as api from '../package-api';
import {HitDesignApp} from '../app/hit-design-app';
import {designCampaign, designCsv, noComputeTemplate, Sandbox} from './fixtures';

type Query = (name: string, ...args: unknown[]) => Promise<DG.DataFrame>;
type SessionInternals = {_session?: {host: {scope: {isDisposed: boolean}}, needsRefresh: {peek(): boolean}}};

const session = (app: HitDesignApp) => (app as unknown as SessionInternals)._session!;

/** The next GetLastModified call of the sandbox — a campaign view's first read of the last save — runs [first]
 * instead, given the real call; the calls after it run as they are. Answers the undo. */
function replaceFirstLastModified(sb: Sandbox,
  first: (query: () => Promise<DG.DataFrame>) => Promise<DG.DataFrame>): () => void {
  const data = Object.getPrototypeOf(grok.data) as {query: Query};
  const query = data.query;
  let replaced = false;
  data.query = function(this: unknown, name: string, ...args: unknown[]) {
    if (replaced || name !== 'HitTriage:GetLastModified' || (args[0] as {appName?: string})?.appName !== sb.appName)
      return query.call(this, name, ...args);
    replaced = true;
    return first(() => query.call(this, name, ...args));
  };
  return () => data.query = query;
}

/** A campaign of the sandbox on {@link designCsv}, open in a new Hit Design app. */
async function openCampaign(sb: Sandbox): Promise<HitDesignApp> {
  const t = noComputeTemplate('ZZB3 Shared', 'ZZB3');
  await sb.writeTemplate(t);
  await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1}), designCsv);
  const app = await sb.hitDesign();
  await sb.openDesign(app, 'ZZB3-1');
  return app;
}

/** Another session's save as the campaign view sees it: the lock taken and let go, which logs the save time. */
async function anotherSessionSaves(sb: Sandbox): Promise<void> {
  await api.queries.acquireCampaignLock(sb.appName, 'ZZB3-1', 'another session');
  await api.queries.releaseCampaignLock(sb.appName, 'ZZB3-1');
}

const ownerDisposedReports = () => DG.U2.Usage.reports().filter((r) => r.code === 'owner-disposed')
  .reduce((n, r) => n + r.count, 0);

category('Campaign view: other sessions', () => {
  test('a failed first read of the last save still starts the check for other sessions\' saves', async () => {
    const sb = await Sandbox.create('poll');
    let failed = false;
    const restore = replaceFirstLastModified(sb, async () => {
      failed = true;
      throw new Error('ZZB3 refused the first read');
    });
    try {
      const app = await openCampaign(sb);
      await waitUntil(() => failed, 'the campaign view never read the last save', 10000);
      restore();
      await anotherSessionSaves(sb);
      await waitUntil(() => session(app).needsRefresh.peek(),
        'the view never noticed another session\'s save', 10000);
    } finally {
      restore();
      await sb.dispose();
    }
  }, {timeout: 120000});

  test('a view closed while it reads the last save starts no check; the next view notices other saves', async () => {
    const sb = await Sandbox.create('poll');
    let letGo!: () => void;
    const gate = new Promise<void>((resolve) => letGo = resolve);
    let reading = false;
    let read = false;
    const restore = replaceFirstLastModified(sb, async (query) => {
      reading = true;
      await gate;
      const df = await query();
      read = true;
      return df;
    });
    try {
      const app = await openCampaign(sb);
      await waitUntil(() => reading, 'the campaign view never read the last save', 10000);
      const closed = session(app);
      await sb.openDesign(app, 'ZZB3-1');
      await waitUntil(() => closed.host.scope.isDisposed, 'the first view\'s session did not end with its view', 10000);
      const disposedBefore = ownerDisposedReports();
      letGo();
      await waitUntil(() => read, 'the closed view\'s read never ended', 10000);
      await anotherSessionSaves(sb);
      await waitUntil(() => session(app).needsRefresh.peek(),
        'the open view never noticed another session\'s save', 10000);
      expect(ownerDisposedReports(), disposedBefore,
        'u2 owner-disposed reports: the closed view\'s session started no check');
    } finally {
      letGo();
      restore();
      await sb.dispose();
    }
  }, {timeout: 120000});
});
