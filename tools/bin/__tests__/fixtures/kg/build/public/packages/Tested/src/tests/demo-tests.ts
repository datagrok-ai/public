import {category, test, expect} from '@datagrok-libraries/utils/src/test';

// test('commented out', async () => {});
/* test('also out', async () => {}); */

const re = /a/;
const name = 'dynamic';

category('~domains/bio', () => {
  test('detects sequences', async () => {
    expect(re.test('a'), true);
  });

  test('renders slowly', async () => {
    const options = {retries: 2};
  }, {timeout: 60000, benchmark: true, tags: ['slow', 'render']});

  test('skipped one', async () => {
  }, {skipReason: 'GROK-100: flaky on CI'});
});

category('Tested: Utils', () => {
  test('tagged', async () => {
    await helper(1, {a: 1});
  }, {tags: ['~domains/bio', 'x']});

  test('conditional skip', async () => {}, {skipReason: typeof process !== 'undefined' ? 'NodeJS environment' : undefined});

  test(`template ${name}`, async () => {});
});

async function helper(n: number, o: object): Promise<void> {}
