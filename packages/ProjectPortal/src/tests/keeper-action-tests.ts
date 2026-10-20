import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, syncPortfolio} from '../portal';
import {recordedSync} from '../sync-runs';
import {buildPortfolioWorkspace} from '../workspace';
import {exportPortfolioConfig} from '../profile';
import {autoAdvance} from '../auto-advance';
import {moveToSecuredTables} from '../secured-rows';
import {KEEPER_ACTIONS, keeperRefusal, showsBoardButton} from '../model/keeper-actions';
import {ADVANCED_BOARD_BUTTONS} from '../model/personas';
import {withStandUser} from './stand-users';
import {_package} from '../package-test';

category('Portal: who keeps a portfolio (rules)', () => {
  test('sync, workspace, park, export and advance by data are the approvers\'', async () => {
    expect(keeperRefusal('sync', 'EX-SMD', false), 'Only EX-SMD\'s approvers sync it');
    expect(keeperRefusal('export', 'EX-SMD', true), null);
    expect(Object.keys(KEEPER_ACTIONS).join(), 'sync,workspace,park,export,auto_advance');
  });

  test('a viewer sees no Sync, Workspace or Export; Simple mode hides the advanced ones from everyone', async () => {
    const shown = (simple: boolean, isApprover: boolean) => ['sync', 'sync_history', 'workspace', 'export', 'all_portfolios']
      .filter((k) => showsBoardButton(k, {simple, advanced: ADVANCED_BOARD_BUTTONS, isApprover})).join();
    expect(shown(false, false), 'sync_history,all_portfolios');
    expect(shown(false, true), 'sync,sync_history,workspace,export,all_portfolios');
    expect(shown(true, true), 'sync', 'Simple mode keeps Sync for an approver');
  });
});

/** The error a call raised, or 'allowed'. */
async function refusal(f: () => Promise<unknown>): Promise<string> {
  try {
    await f();
    return 'allowed';
  } catch (e: any) {
    return String(e?.message ?? e);
  }
}

category('Portal: a viewer cannot keep a portfolio', () => {
  before(async () => {
    await loadAllPortfolios(_package);
  });

  test('the functions refuse a viewer, before writing anything', async () => {
    const smd = await portfolioByCode('EX-SMD');
    const runs = async () => (await portalDb.syncLogs.query({filter: `portfolio_id = "${smd.id}"`, limit: 10000})).length;
    const before = await runs();
    const tried = await withStandUser('viewer', (u) => u.asUser(async () => [
      await refusal(() => recordedSync('EX-SMD', 'manual')),
      await refusal(() => syncPortfolio('EX-SMD')),
      await refusal(() => buildPortfolioWorkspace('EX-SMD')),
      await refusal(() => exportPortfolioConfig('EX-SMD')),
      await refusal(() => autoAdvance('EX-SMD', false)),
      await refusal(() => moveToSecuredTables('EX-SMD', true)),
    ]));
    expect(tried !== null, true, 'no plain user could be had on this stand');
    expect(tried!.slice(0, 5).join('\n'), [
      'Only EX-SMD\'s approvers sync it', 'Only EX-SMD\'s approvers sync it', 'Only EX-SMD\'s approvers build its workspace',
      'Only EX-SMD\'s approvers export it', 'Only EX-SMD\'s approvers advance its programs by data',
    ].join('\n'));
    expect(tried![5].startsWith('Only admins move records'), true, tried![5]);
    expect(await runs(), before, 'a refused sync wrote a history row');
    // dry runs read only: still a viewer's
    const dry = await withStandUser('viewer', (u) => u.asUser(() => refusal(() => buildPortfolioWorkspace('EX-SMD', true))));
    expect(dry, 'allowed');
  });
});
