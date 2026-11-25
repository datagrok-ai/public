import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {waitUntil, category, delay, expect, test} from '@datagrok-libraries/test/src/test';
import {initTemplates} from '../search/templates-search';
import {appSearch, exactAppFuncSearch, functionSearch, helpSearch} from '../search/entity-search';
import {createFuncTableViewWidget, powerSearch, searchSuggestions} from '../search/power-search';

/** What the menu under the Home search box lists for the text, top to bottom. A test session loads the package
 * again for its test module (js_func.dart, _loadPackage), so the handlers its init registers come twice. */
async function suggestions(s: string): Promise<string> {
  return [...new Set((await searchSuggestions(s)).sort((a, b) => (a.priority ?? 0) - (b.priority ?? 0))
    .map((sug) => `${sug.text.trim()} -> ${sug.value ?? ''}`))].join(' | ');
}

/** The Home page's search for the text, run into a host of its own until every path has answered: the category
 * headers and the text it shows, and the paths that failed. Its suggestion menu stays up until the next search. */
async function search(s: string): Promise<{categories: string[], text: string, failures: string[]}> {
  const failures: string[] = [];
  const onRejection = (e: PromiseRejectionEvent) => failures.push(`${e.reason}`);
  window.addEventListener('unhandledrejection', onRejection);
  try {
    const input = ui.element('input') as HTMLInputElement;
    input.value = s;
    ui.div([input]);
    const host = ui.div();
    powerSearch(s, host, input);
    await waitUntil(() => host.getAttribute('aria-busy') === 'false', `the search for "${s}" has not finished`, 30000);
    await delay(10);
    return {categories: [...host.querySelectorAll('.power-pack-search-list-header')].map((h) => h.textContent ?? ''),
      text: (host.textContent ?? '').trim(), failures};
  }
  finally {
    window.removeEventListener('unhandledrejection', onRejection);
  }
}

category('Search', () => {
  test('InitTemplates', async () => await initTemplates());

  test('AppSearchIgnoresSpaces', async () => {
    const app = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}}).find((f) => f.friendlyName?.includes(' '));
    if (app == null)
      return;
    const found = await appSearch(app.friendlyName.replaceAll(' ', ''));
    expect(found.some((f) => f.id === app.id), true);
  });

  test('ExactAppSearchIgnoresSpaces', async () => {
    const app = DG.Func.find({meta: {role: DG.FUNC_TYPES.APP}, returnType: 'view'})
      .find((f) => f.friendlyName?.includes(' '));
    if (app == null)
      return;
    expect(exactAppFuncSearch(app.friendlyName.replaceAll(' ', ''))?.id, app.id);
  });

  test('NewSuggestsNewUsersAndFindsFunctionsAndHelp', async () => {
    expect(await suggestions('new'), 'New Users Today ->  | New users This Month ->  | New users This Year ->  | ' +
      'New users last 3 months ->  | New user last 7 days ->  | New users yesterday -> ');
    expect((await functionSearch('new')).some((f) => f.friendlyName === 'Add New Column'), true,
      'Add New Column among the functions found');
    expect((await helpSearch('new')).length > 0, true, 'help pages found');
  });

  test('UseSuggestsUserLoginAndFindsFunctionsAndHelp', async () => {
    expect(await suggestions('use'), 'DGUSER-{User Login} -> DGUSER-');
    expect((await functionSearch('use')).length > 0, true, 'functions found');
    expect((await helpSearch('use')).length > 0, true, 'help pages found');
  });

  // A formula is evaluated only with parentheses, and a pattern is searched for as text.
  test('FormulaAndPatternAreSearchedAsText', async () => {
    expect(await suggestions('1+1'), 'PDB ID, e.g. 4AKZ -> 4AKZ');
    const sum = await search('1+1');
    expect(sum.failures.join('; '), '', 'the paths of the search for "1+1" that failed');
    expect(sum.categories.join(', '), '', 'the categories the search for "1+1" shows');
    expect(sum.text, '', 'what the search for "1+1" shows');

    expect(await suggestions('Project[0-9]+'), '');
    expect((await search('Project[0-9]+')).failures.join('; '), '',
      'the paths of the search for "Project[0-9]+" that failed');
  });

  test('TableResultFillsItsWidget', async () => {
    const script = DG.Script.create('//name: Table\n//language: javascript\n//output: dataframe df\n' +
      'df = grok.data.demo.demog(10);');
    const widget = createFuncTableViewWidget(script, {});
    const view = grok.shell.newView('table search result', [ui.block([widget], 'power-pack-search-host')]);
    const height = (selector: string) => widget.querySelector(selector)?.clientHeight;
    const gridFills = (h: number) => waitUntil(() => widget.clientHeight === h &&
      height('.d4-grid canvas') !== undefined && height('.d4-grid') === height('.power-pack-widget-content'),
    () => `widget ${widget.clientHeight} of ${h} px: the grid takes ${height('.d4-grid')} of its ` +
      `${height('.power-pack-widget-content')} px`);
    try {
      await gridFills(window.innerHeight - 200);
      for (const h of [900, 350]) {
        widget.style.setProperty('height', `${h}px`, 'important');
        await gridFills(h);
      }
    }
    finally {
      view.close();
    }
  });
});
