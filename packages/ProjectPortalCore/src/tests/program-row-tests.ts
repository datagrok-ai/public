import * as grok from 'datagrok-api/grok';
import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {portalDb} from '../generated/db';
import {loadAllPortfolios, portfolioByCode, syncPortfolio} from '../portal';
import {describeRowLinks, planRowLinks} from '../model/program-rows';
import {folderLinks, linkProgramRows} from '../program-rows';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

category('Portal: program rows through their folders', () => {
  test('plan: link what has a folder, skip the rest and say why', async () => {
    const links = new Map([['s1', new Set(['p2'])], ['s2', new Set<string>()]]);
    const plan = planRowLinks([
      {id: 'p1', code: 'A', spaceId: 's2'},
      {id: 'p2', code: 'B', spaceId: 's1'},
      {id: 'p3', code: 'C', spaceId: null},
      {id: 'p4', code: 'D', spaceId: 'gone'},
    ], links);
    expect(plan.link.map((p) => p.code).join(), 'A');
    expect(plan.linked.map((p) => p.code).join(), 'B');
    expect(describeRowLinks(plan, true).join('\n'), [
      'would link 1 program row(s) into their folders: A',
      '1 already linked',
      'C: no workspace folder: private to its author',
      'D: its workspace folder is gone',
    ].join('\n'));
  });
});

category('Portal: a plain viewer sees their portfolio\'s programs', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    await syncPortfolio('EX-SMD');
    await syncPortfolio('EX-BIOD');
  });

  test('a viewer sees the linked programs of their portfolio, read-only, and no other portfolio\'s', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const bd = await portfolioByCode('EX-BIOD');
    await linkProgramRows('EX-SMD');
    const mine = (await portalDb.programs.query({filter: `portfolio_id = "${sm.id}"`, limit: 1000}))
      .filter((p) => p.space_id);
    const links = await folderLinks(mine.map((p) => p.space_id!));
    const linked = mine.filter((p) => links.get(p.space_id!)?.has(p.id)).map((p) => p.code).sort();
    expect(linked.length > 0, true, 'no EX-SMD program has a folder: build the EX-SMD workspace');
    const ran = await withRestrictedUser('pprow', async (u) => {
      expect(await u.asUser(() => portalDb.programs.count({})), 0, 'a user in no portfolio group saw programs');
      const viewers = await grok.dapi.groups.find(sm.viewers!);
      const personal = await grok.dapi.groups.find(u.group);
      await grok.dapi.groups.includeTo(personal, viewers);
      try {
        const seen = await u.asUser(() => portalDb.programs.query({limit: 1000}));
        expect(seen.every((p) => p.portfolio_id === sm.id), true, `saw another portfolio's: ${seen.map((p) => p.code).join()}`);
        expect(seen.map((p) => p.code).sort().join(), linked.join());
        expect(await u.asUser(() => portalDb.programs.count({filter: `portfolio_id = "${bd.id}"`})), 0);
        // a link reaches the row for View only: the viewer cannot rename it
        let refused = false;
        await u.asUser(() => portalDb.programs.update(seen[0].id, {name: 'renamed by a viewer'})).catch(() => refused = true);
        expect(refused, true, 'a viewer wrote a program row');
      } finally {
        await grok.dapi.groups.excludeFrom(personal, viewers);
      }
    });
    expect(ran !== null, true, 'no restricted session could be created on this stand');
  });

  test('linking is idempotent', async () => {
    await linkProgramRows('EX-SMD');
    const again = await linkProgramRows('EX-SMD', true);
    expect(again.some((l) => l.startsWith('would link')), false, again.join('; '));
  });
});
