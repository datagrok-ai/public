import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {FollowRow, portalDb} from '../generated/db';
import {portfolioByCode} from '../portal';
import {feedOfFollows} from '../feed';
import {provenanceRecords} from '../provenance';
import {inFilter, queryIn} from '../rows';

category('Portal: lists in filters', () => {
  test('an empty list makes no filter and asks nothing', async () => {
    expect(inFilter('id', []), null);
    expect(inFilter('id', [null, undefined, '']), null);
    expect(inFilter('id', ['a', 'b', 'a']), 'id in ("a", "b")');
    expect(inFilter('code', ['say "hi"']), 'code in ("say \\"hi\\"")');
    let asked = 0;
    const table = {query: async () => { asked++; return [1]; }};
    expect((await queryIn(table, 'id', [])).length, 0);
    expect(asked, 0, 'an empty list went to the server');
    expect((await queryIn(table, 'id', ['x'], 'status = "ok"')).length, 1);
  });

  test('my feed over a followed program that is gone is empty, not an error', async () => {
    const sm = await portfolioByCode('EX-SMD');
    const code = `GONE-${Date.now() % 1e9}`;
    const [ins] = await portalDb.programs.insert({portfolio_id: sm.id, code, name: code,
      stage_model_id: sm.stage_model_id, source: 'manual'});
    await portalDb.programs.delete(ins.id);
    // the follow as myFollowsAll reads it; not written, so nobody's feed changes
    const follow = {id: 'f-gone', user_id: 'u', program_id: ins.id, active: true} as unknown as FollowRow;
    const feed = await feedOfFollows([follow]);
    expect(feed.programs, 0);
    expect(feed.entries.length, 0);
    const records = await provenanceRecords([ins.id]);
    expect(records.programs.length, 0, 'provenance of a gone program');
  });
});
