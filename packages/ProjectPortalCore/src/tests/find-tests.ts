import {before, category, expect, test} from '@datagrok-libraries/test/src/test';
import {loadAllPortfolios, syncPortfolio} from '../portal';
import {addRecent, FindProgram, findPrograms, findTerms, matchedBy, MAX_RECENT, parseRecents, resolveRecents} from '../model/find';
import {findablePrograms, recentPrograms, rememberProgram} from '../find';
import {readPreference, writePreference} from '../personas';
import {withRestrictedUser} from './restricted-user';
import {_package} from '../package-test';

const prog = (code: string, name: string, extra: Partial<FindProgram> = {}): FindProgram =>
  ({id: `id-${code}`, code, name, portfolio: 'EX', portfolioName: 'Example', ...extra});

category('Portal: find a program', () => {
  const programs = [
    prog('EXS-KIN1', 'EXS-KIN1', {target: 'Kinase K1'}),
    prog('EXS-KIN12', 'Kinase follow-up'),
    prog('EXB-201', 'Antibody Alpha', {portfolio: 'EXB', target: 'ALP-1', details: {molecule_code: 'MC-778', count: 3,
      note: 'a long free-text remark that is prose, not an alias, and runs past the limit'}}),
  ];

  test('terms: code, name, target and short details values', async () => {
    const terms = findTerms(programs[2]).map((t) => `${t.field}=${t.value}`);
    expect(terms.join('|'), 'code=EXB-201|name=Antibody Alpha|target=ALP-1|molecule code=MC-778|count=3');
  });

  test('best match first: the exact code, then the code\'s start, then a word', async () => {
    expect(findPrograms('exs-kin1', programs).map((h) => h.program.code).join(), 'EXS-KIN1,EXS-KIN12');
    expect(findPrograms('exs', programs).map((h) => h.program.code).join(), 'EXS-KIN1,EXS-KIN12');
    expect(findPrograms('follow', programs).map((h) => h.program.code).join(), 'EXS-KIN12', 'a word of the name');
    expect(findPrograms('alpha', programs).map((h) => h.program.code).join(), 'EXB-201');
    expect(findPrograms('k', programs).length, 0, 'one letter finds nothing');
  });

  test('an alias finds its program and says which', async () => {
    const [hit] = findPrograms('mc 778', programs);
    expect(hit.program.code, 'EXB-201');
    expect(matchedBy(hit), 'molecule code: MC-778');
    expect(findPrograms('prose', programs).length, 0, 'long text is not an alias');
    expect(matchedBy(findPrograms('antibody', programs)[0]), null, 'the name says it');
  });

  test('recents: newest first, once each, capped; unreadable or unseen ones dropped', async () => {
    let list = parseRecents('not json');
    expect(list.length, 0);
    for (let i = 0; i < MAX_RECENT + 2; i++)
      list = addRecent(list, {id: `p${i}`, portfolio: 'EX', at: String(i)});
    list = addRecent(list, {id: 'p5', portfolio: 'EX', at: 'again'});
    expect(list.length, MAX_RECENT);
    expect(list.slice(0, 2).map((r) => r.id).join(), 'p5,p9');
    expect(parseRecents(JSON.stringify([{id: 'a', portfolio: 'X'}, {bad: 1}])).length, 1);
    const seen = resolveRecents([{id: 'id-EXB-201', portfolio: 'EXB', at: ''}, {id: 'gone', portfolio: 'EX', at: ''}], programs);
    expect(seen.map((p) => p.code).join(), 'EXB-201');
  });
});

category('Portal: find a program on the stand', () => {
  before(async () => {
    await loadAllPortfolios(_package);
    await syncPortfolio('EX-SMD');
    await syncPortfolio('EX-BIOC');
  });

  test('finds across portfolios by code, target and name', async () => {
    const all = await findablePrograms();
    expect(findPrograms('EXS-KIN1', all)[0]?.program.portfolio, 'EX-SMD');
    const byTarget = findPrograms('ALP-1', all).find((h) => h.program.code === 'EXB-201');
    expect(byTarget?.program.portfolio, 'EX-BIOC');
    expect(matchedBy(byTarget!), 'target: ALP-1');
    expect(findPrograms('Antibody Alpha', all).some((h) => h.program.code === 'EXB-201'), true);
  });

  test('a user finds only what they can see', async () => {
    const ran = await withRestrictedUser('ppfind', async (u) => {
      expect((await u.asUser(() => findablePrograms())).length, 0, 'a user in no portfolio group found programs');
    });
    expect(ran !== null, true, 'no restricted session could be created on this stand');
  });

  test('opening a program puts it first in the user\'s recents', async () => {
    const before = readPreference('recent');
    try {
      const all = await findablePrograms();
      const kin = all.find((p) => p.code === 'EXS-KIN1')!;
      const alpha = all.find((p) => p.code === 'EXB-201')!;
      rememberProgram(kin.portfolio, kin.id);
      rememberProgram(alpha.portfolio, alpha.id);
      rememberProgram(kin.portfolio, kin.id);
      expect(resolveRecents(recentPrograms(), all).slice(0, 2).map((p) => p.code).join(), 'EXS-KIN1,EXB-201');
    } finally {
      writePreference('recent', before);
    }
  });
});
