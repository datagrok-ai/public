/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {layoutObjects} from '../app/shared/campaign-session';
import {designCampaign, designCsv, expectSame, hasBalloon, noComputeTemplate, Sandbox, startDesignCampaign} from './fixtures';

/** The columns the Filters viewers of a layout filter on. */
function filterColumns(viewState: string): string[] {
  return layoutObjects(viewState).filter((o) => o['#type'] === 'FiltersLook')
    .flatMap((o) => ((o.filters ?? []) as {column?: string}[]).map((f) => f.column as string));
}

/** A layout whose Filters viewer filters on Stage and on two descriptors a new campaign calculates later. */
async function filtersLayout(): Promise<string> {
  const df = DG.DataFrame.fromCsv('Molecule,Stage,V-iD,HeavyAtomCount,qed\nCCO,Idea,,3,0.4\n');
  const tv = grok.shell.addTableView(df);
  try {
    const filters = tv.getFiltersGroup({createDefaultFilters: false});
    filters.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Stage'});
    filters.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'HeavyAtomCount'});
    filters.updateOrAdd({type: DG.FILTER_TYPE.HISTOGRAM, column: 'qed'});
    let state = '';
    await waitUntil(() => filterColumns(state = tv.saveLayout().viewState).length === 3,
      'the filters never reached the layout', 10000);
    return state;
  } finally {
    tv.close();
    grok.shell.closeTable(df);
  }
}

category('Campaign view: layout', () => {
  test('a new campaign whose template layout filters on columns it calculates later opens without errors and keeps those filters', async () => {
    const sb = await Sandbox.create('layout');
    try {
      const t = {...noComputeTemplate('ZZB3 Layout', 'ZZB3'), layoutViewState: await filtersLayout(),
        compute: {descriptors: {enabled: true, args: ['HeavyAtomCount', 'qed']}, functions: [], scripts: [], queries: []}};
      await sb.writeTemplate(t);
      const app = await sb.hitDesign();
      const id = await startDesignCampaign(app, t, {Chemist: 'Ann'});
      expect(hasBalloon('Error adding filter'), false, 'an "Error adding filter" balloon');
      expect(document.querySelector('[data-u2-item^="Apply template layout"] i')?.classList.contains('svg-view-layout'), true,
        'Apply template layout shows the platform\'s view-layout picture');
      expectSame(filterColumns((await sb.readCampaign(id)).layout).sort(), ['HeavyAtomCount', 'Stage', 'qed'],
        'the filters the created campaign keeps');
      await app.performSingleCellCalculations(0, 'CCO');
      await waitUntil(() => app.campaign?.version === 2 && !app.isSaving, 'the calculated row was not saved', 60000);
      expect(app.dataFrame!.get('HeavyAtomCount', 0), 3, 'the calculated HeavyAtomCount');
      expectSame(filterColumns((await sb.readCampaign(id)).layout).sort(), ['HeavyAtomCount', 'Stage', 'qed'],
        'the filters kept through the first calculation');
    } finally {
      await sb.dispose();
    }
  }, {timeout: 180000});

  test('Add new row in the ribbon adds a row and logs no error', async () => {
    const sb = await Sandbox.create('layout');
    const logged: string[] = [];
    const error = console.error;
    try {
      const t = noComputeTemplate('ZZB3 Layout', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t, {version: 1}), designCsv);
      const app = await sb.hitDesign();
      const view = await sb.openDesign(app, 'ZZB3-1');
      const add = view.root.querySelector<HTMLElement>('[data-u2-item="Add new row"]') ??
        document.querySelector<HTMLElement>('[data-u2-item="Add new row"]');
      expect(add != null, true, 'the ribbon\'s Add new row');
      console.error = (...args: unknown[]) => {
        logged.push(args.map(String).join(' '));
        error(...args);
      };
      add!.click();
      await waitUntil(() => app.dataFrame!.rowCount === 3, 'no row was added', 10000);
    } finally {
      console.error = error;
      await sb.dispose();
    }
    expectSame(logged.filter((l) => l.includes('No such column')), [], 'errors logged by the click');
  }, {timeout: 120000});
});
