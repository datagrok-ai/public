/** Toolbox top tabs (Files / Queries / Workflows / Favorites) and the
 *  star-driven favorites store: tab structure, workflow listing, star toggling
 *  from a catalog row, localStorage persistence, and the Favorites tab
 *  round-trip (star → listed → create node → unstar → gone). */
import * as DG from 'datagrok-api/dg';
import {category, test, expect, before, after} from '@datagrok-libraries/utils/src/test';

import {registerBuiltinNodes, registerAllFunctions, getRegisteredFuncs, isWorkflowFunc} from '../rete/node-factory';
import {FunctionBrowser, TOOLBOX_TABS} from '../panel/function-browser';
import {getFavorites, isFavorite, toggleFavorite, clearFavorites, onFavoritesChanged} from '../panel/favorites';
import {supportedUploadExtensions} from '../utils/uploaded-files';

function makeBrowser(onAdd?: (type: string) => void, onFiles?: (files: File[]) => void): FunctionBrowser {
  return new FunctionBrowser({
    onFunctionDoubleClick: (info) => onAdd?.(info.nodeTypeName),
    onBuiltinNodeDoubleClick: (type) => onAdd?.(type),
    onFileDoubleClick: () => {},
    onLocalFilesPicked: (files) => onFiles?.(files),
  });
}

category('Flow: toolbox tabs', () => {
  before(async () => {
    registerBuiltinNodes();
    registerAllFunctions();
    clearFavorites();
  });

  after(async () => {
    clearFavorites();
  });

  test('the top strip has the four collection tabs', async () => {
    const browser = makeBrowser();
    document.body.appendChild(browser.root);
    try {
      browser.render();
      expect(!!browser.topTabs, true, 'tab control created');
      const names = browser.topTabs!.panes.map((p) => p.name);
      expect(names.join(','), TOOLBOX_TABS.join(','), 'Files, Queries, Workflows, Favorites in order');
      for (const name of TOOLBOX_TABS) {
        expect(!!browser.root.querySelector(`[data-testid="ff-browser-tab-${name.toLowerCase()}"]`),
          true, `${name} tab header carries its test id`);
      }
    } finally {
      browser.root.remove();
      browser.destroy();
    }
  });

  test('the Workflows tab lists saved flows (and only flows)', async () => {
    const browser = makeBrowser();
    document.body.appendChild(browser.root);
    try {
      browser.render();
      browser.showTab('Workflows');
      const tab = browser.root.querySelector('[data-testid="ff-browser-workflows"]') as HTMLElement;
      expect(!!tab, true, 'Workflows tab content present');
      const flows = getRegisteredFuncs().filter((f) => isWorkflowFunc(f.func));
      const items = Array.from(tab.querySelectorAll('.funcflow-func-item')) as HTMLElement[];
      expect(items.length, flows.length, 'one row per saved flow');
      // Workflows are OUT of the accordion below.
      expect(browser.accordion!.panes.some((p) => p.name === 'Workflows'), false,
        'no Workflows accordion section');
    } finally {
      browser.root.remove();
      browser.destroy();
    }
  });

  test('favorites store: toggle, persistence, change notifications', async () => {
    clearFavorites();
    let notified = 0;
    const unsub = onFavoritesChanged(() => notified++);
    try {
      expect(toggleFavorite({type: 'Inputs/Table Input', label: 'Table Input'}), true, 'starred');
      expect(isFavorite('Inputs/Table Input'), true);
      expect(notified > 0, true, 'listeners notified');
      // Round-trips through localStorage.
      const raw = localStorage.getItem('funcflow.favorites.v1');
      expect(!!raw && raw.includes('Table Input'), true, 'persisted to localStorage');
      expect(toggleFavorite({type: 'Inputs/Table Input', label: 'Table Input'}), false, 'unstarred');
      expect(isFavorite('Inputs/Table Input'), false);
      expect(getFavorites().length, 0);
    } finally {
      unsub();
      clearFavorites();
    }
  });

  test('clicking a row star pins the node into the Favorites tab', async () => {
    clearFavorites();
    const added: string[] = [];
    const browser = makeBrowser((t) => added.push(t));
    document.body.appendChild(browser.root);
    try {
      browser.render();
      // Materialize the Inputs section, then star Table Input via its row star.
      browser.accordion!.getPane('Inputs').expanded = true;
      const star = browser.root.querySelector(
        '[data-testid="ff-browser-item-star-inputs-table-input"]') as HTMLElement;
      expect(!!star, true, 'row star rendered');
      star.click();
      expect(isFavorite('Inputs/Table Input'), true, 'star click favorites the type');
      expect(star.classList.contains('funcflow-item-star-active'), true, 'star repainted as active');

      // The Favorites tab lists it; double-click creates the node type.
      browser.showTab('Favorites');
      const fav = browser.root.querySelector(
        '[data-testid="ff-browser-fav-item-inputs-table-input"]') as HTMLElement;
      expect(!!fav, true, 'favorite listed in the Favorites tab');
      expect(fav.textContent!.includes('Table Input'), true, 'favorite keeps its label');
      fav.dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      expect(added.join(','), 'Inputs/Table Input', 'double-click adds the starred node type');

      // Unstar from the Favorites tab row → the entry disappears in place.
      const favStar = fav.querySelector('.funcflow-item-star') as HTMLElement;
      favStar.click();
      expect(isFavorite('Inputs/Table Input'), false, 'unstarred from the Favorites tab');
      expect(browser.root.querySelector('[data-testid="ff-browser-fav-item-inputs-table-input"]') == null,
        true, 'entry removed from the Favorites tab');
    } finally {
      browser.root.remove();
      browser.destroy();
      clearFavorites();
    }
  });

  test('searching badges the tab headers with match counts and dims 0-match tabs', async () => {
    clearFavorites();
    const browser = makeBrowser();
    document.body.appendChild(browser.root);
    try {
      browser.render();
      const input = browser.root.querySelector('[data-testid="ff-browser-search"]') as HTMLInputElement;
      const qHeader = browser.root.querySelector('[data-testid="ff-browser-tab-queries"]') as HTMLElement;

      input.value = 'zzzz-no-such-thing';
      input.dispatchEvent(new Event('input', {bubbles: true}));
      expect(qHeader.querySelector('.funcflow-tab-badge') == null, true, '0 matches → no badge, just dim');
      expect(qHeader.classList.contains('funcflow-tab-dim'), true, '0-match tab dims');

      // A query with real query matches shows a >0 count badge.
      const queries = getRegisteredFuncs().filter((f) => f.func instanceof DG.DataQuery);
      if (queries.length > 0) {
        input.value = queries[0].name.slice(0, 6);
        input.dispatchEvent(new Event('input', {bubbles: true}));
        const badge = qHeader.querySelector('.funcflow-tab-badge');
        expect(badge != null && Number(badge.textContent) > 0, true, 'matching tab shows its count');
        expect(qHeader.classList.contains('funcflow-tab-dim'), false, 'matching tab not dimmed');
      }

      input.value = '';
      input.dispatchEvent(new Event('input', {bubbles: true}));
      expect(qHeader.querySelector('.funcflow-tab-badge') == null, true, 'badges clear with the query');
      expect(qHeader.classList.contains('funcflow-tab-dim'), false, 'dim clears with the query');
    } finally {
      browser.root.remove();
      browser.destroy();
    }
  });

  test('group-by is a compact catalog-header button; setGroupBy regroups and persists', async () => {
    const saved = localStorage.getItem('funcflow.browser.v1');
    const browser = makeBrowser();
    document.body.appendChild(browser.root);
    try {
      browser.render();
      const header = browser.root.querySelector('[data-testid="ff-browser-catalog-header"]');
      const btn = browser.root.querySelector('[data-testid="ff-browser-groupby"]') as HTMLElement;
      expect(!!header && !!btn && header!.contains(btn), true, 'button lives in the catalog header');
      expect(btn.textContent!.includes('what it does'), true, 'shows the current mode');

      browser.setGroupBy('package');
      expect(btn.textContent!.includes('package'), true, 'label follows the mode');
      const stored = JSON.parse(localStorage.getItem('funcflow.browser.v1') ?? '{}');
      expect(stored.groupBy, 'package', 'mode persisted');
    } finally {
      browser.root.remove();
      browser.destroy();
      if (saved == null) localStorage.removeItem('funcflow.browser.v1');
      else localStorage.setItem('funcflow.browser.v1', saved);
    }
  });

  test('the Files tab Upload chip hands picked files to the view', async () => {
    const picked: File[][] = [];
    const browser = makeBrowser(undefined, (files) => picked.push(files));
    document.body.appendChild(browser.root);
    try {
      browser.render();
      browser.showTab('Files');
      const btn = browser.root.querySelector('[data-testid="ff-browser-upload"]') as HTMLElement;
      expect(!!btn, true, 'Upload chip rendered in the Files tab');
      const hint = btn.closest('.funcflow-tab-hint');
      expect(!!hint, true, 'chip shares the hint row — no extra vertical space');

      // Simulate the OS picker returning a file (input.files is read-only —
      // a DataTransfer is the standard way to set it programmatically).
      const input = btn.querySelector('input[type="file"]') as HTMLInputElement;
      const dt = new DataTransfer();
      dt.items.add(new File(['a,b\n1,2\n'], 'picked.csv', {type: 'text/csv'}));
      input.files = dt.files;
      input.dispatchEvent(new Event('change', {bubbles: true}));
      expect(picked.length, 1, 'callback fired once');
      expect(picked[0].map((f) => f.name).join(','), 'picked.csv', 'the picked file came through');
      expect(input.value, '', 'input cleared so re-picking the same file fires again');
    } finally {
      browser.root.remove();
      browser.destroy();
    }
  });

  test('supported upload extensions cover the natively-parsed formats', async () => {
    const exts = supportedUploadExtensions();
    for (const e of ['csv', 'tsv', 'txt', 'd42'])
      expect(exts.includes(e), true, `.${e} is uploadable`);
  });

  test('a DG-function favorite resolves back to its FuncInfo', async () => {
    clearFavorites();
    const info = getRegisteredFuncs().find((f) => f.func.name === 'JoinTables') ?? getRegisteredFuncs()[0];
    const added: string[] = [];
    const browser = makeBrowser((t) => added.push(t));
    document.body.appendChild(browser.root);
    try {
      browser.render();
      toggleFavorite({type: info.nodeTypeName, label: info.name});
      browser.showTab('Favorites');
      const rows = Array.from(browser.root.querySelectorAll(
        '[data-testid="ff-browser-favorites"] .funcflow-func-item')) as HTMLElement[];
      expect(rows.length, 1, 'one favorite listed');
      expect(rows[0].dataset.nodeTypeName, info.nodeTypeName, 'row targets the function node type');
      expect(rows[0].dataset.func, info.func.name, 'row resolved the live FuncInfo');
      rows[0].dispatchEvent(new MouseEvent('dblclick', {bubbles: true, cancelable: true}));
      expect(added.join(','), info.nodeTypeName, 'double-click routes through onFunctionDoubleClick');
    } finally {
      browser.root.remove();
      browser.destroy();
      clearFavorites();
    }
  });
});
