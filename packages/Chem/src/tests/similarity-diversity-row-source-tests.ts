import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {waitUntil, before, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {ChemSearchBaseViewer, sameRows} from '../analysis/chem-search-base-viewer';
import {ChemSimilarityViewer} from '../analysis/chem-similarity-viewer';

/** The rows the viewer's cards show, as its status reports them. */
function cardRows(v: ChemSearchBaseViewer): number[] {
  const set = String(v.getWidgetStatus().values?.['card row set'] ?? '');
  return set === '' ? [] : set.split(', ').map(Number);
}

/** Waits for the viewer to search the rows its row source has now and settle on cards that satisfy [holds]. */
async function cardsWhere(v: ChemSearchBaseViewer, claim: string, holds: (rows: number[]) => boolean):
  Promise<number[]> {
  let rows: number[] = [];
  await waitUntil(() => {
    rows = cardRows(v);
    return !v.isRenderPending && sameRows(v.getRowSourceIndexes(), v['searchedRows']) && rows.length > 0 &&
      holds(rows);
  }, () => `${claim}, but the cards show rows ${rows.join(', ')}`, 15000);
  return rows;
}

async function withViewer(type: string, body: (df: DG.DataFrame, v: ChemSearchBaseViewer) => Promise<void>) {
  const df = grok.data.demo.molecules(100);
  await grok.data.detectSemanticTypes(df);
  const tv = grok.shell.addTableView(df);
  try {
    tv.addViewer(type);
    // the menu's addViewer returns before the package has created the viewer, so it hands back the platform's
    // wrapper; the viewer itself comes from the package's bundle, whose classes are not this bundle's
    let v: ChemSearchBaseViewer | undefined;
    await waitUntil(() => (v = Array.from(tv.viewers)
      .find((x) => x.type === type && 'getRowSourceIndexes' in x) as unknown as ChemSearchBaseViewer) !== undefined,
    `${type} was not added`, 10000);
    await cardsWhere(v!, `${type} opens with its cards`, () => true);
    await body(df, v!);
  } finally {
    tv.close();
  }
}

category('similarity/diversity: row source', () => {
  before(async () => {
    grok.shell.closeAll();
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  // the similarity viewer always shows its reference molecule, the current row: row 0 here, which every
  // row set below keeps
  for (const [name, type] of [['similarity', 'Chem Similarity Search'], ['diversity', 'Chem Diversity Search']]) {
    test(`${name}: Filtered is the default row source`, async () => {
      await withViewer(type, async (_, v) => {
        expect(v.rowSource, 'Filtered');
        expect(v.getWidgetStatus().values?.['row source'], 'Filtered');
      });
    });

    test(`${name}: Filtered searches the rows that pass the filter, and follows the filter`, async () => {
      await withViewer(type, async (df, v) => {
        df.filter.init((i) => i % 2 === 0);
        await cardsWhere(v, 'with the odd rows filtered out, every card should show an even row',
          (rows) => rows.every((r) => r % 2 === 0));
        df.filter.init((i) => i % 3 === 0);
        await cardsWhere(v, 'with every third row passing the filter, every card should show one of them',
          (rows) => rows.every((r) => r % 3 === 0));
      });
    });

    test(`${name}: Selected searches the selected rows, and follows the selection`, async () => {
      await withViewer(type, async (df, v) => {
        v.setOptions({rowSource: 'Selected'});
        df.selection.init((i) => i % 2 === 0);
        await cardsWhere(v, 'with the even rows selected, every card should show an even row',
          (rows) => rows.every((r) => r % 2 === 0));
        df.selection.init((i) => i % 3 === 0);
        await cardsWhere(v, 'with every third row selected, every card should show one of them',
          (rows) => rows.every((r) => r % 3 === 0));
      });
    });

    test(`${name}: FilteredSelected searches the selected rows that pass the filter`, async () => {
      await withViewer(type, async (df, v) => {
        v.setOptions({rowSource: 'FilteredSelected'});
        df.filter.init((i) => i % 2 === 0);
        df.selection.init((i) => i % 3 === 0);
        await cardsWhere(v, 'with the even rows passing the filter and every third row selected, ' +
          'every card should show a row that is both', (rows) => rows.every((r) => r % 6 === 0));
      });
    });

    test(`${name}: All ignores the filter`, async () => {
      await withViewer(type, async (df, v) => {
        v.setOptions({rowSource: 'All'});
        const before = await cardsWhere(v, 'with the row source All, the viewer should show its cards', () => true);
        const shown = new Set(before);
        df.filter.init((i) => i === 0 || !shown.has(i));
        await delay(1000);
        await cardsWhere(v,
          `with the rows of the cards filtered out, the cards should stay on rows ${before.join(', ')}`,
          (rows) => rows.join(', ') === before.join(', '));
      });
    });
  }

  test('similarity: Select all similar selects the similar rows of the row source', async () => {
    await withViewer('Chem Similarity Search', async (df, v) => {
      df.filter.init((i) => i % 2 === 0);
      await cardsWhere(v, 'with the odd rows filtered out, every card should show an even row',
        (rows) => rows.every((r) => r % 2 === 0));
      (v as ChemSimilarityViewer).selectSimilar();
      const selected = Array.from(df.selection.getSelectedIndexes());
      expect(selected.length > 0, true, 'Select all similar selected nothing');
      const odd = selected.filter((r) => r % 2 === 1);
      expect(odd.length, 0, `Select all similar selected the filtered-out rows ${odd.join(', ')}`);
    });
  });
});
