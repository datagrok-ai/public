import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';

const nodeSkip = typeof process !== 'undefined' ? 'NodeJS environment' : undefined;

/** The crumbs and the name the ribbon shows for the current view. */
function viewRibbon(): HTMLElement | null {
  return document.querySelector('.d4-view-ribbon');
}

function parentCrumb(): HTMLElement | null {
  return viewRibbon()?.querySelector(':scope > span.d4-ribbon-name') ?? null;
}

function viewName(): HTMLElement | null {
  return viewRibbon()?.querySelector('.grok-ns-entity-name') ?? null;
}

function saveButton(tv: DG.TableView): HTMLElement | null {
  return tv.root.parentElement!.querySelector('.d4-ribbon-panel[name="Save"] button');
}

const isShown = (e: HTMLElement | null): boolean => e != null && e.style.display !== 'none';

/** The same platform view, whichever wrapper each side holds. */
const sameView = (a: DG.ViewBase | null | undefined, b: DG.ViewBase): boolean => a?.dart === b.dart;

const isOpen = (v: DG.ViewBase): boolean => Array.from(grok.shell.views).some((x) => sameView(x, v));

const textLook = (e: HTMLElement): string => `${getComputedStyle(e).color} ${getComputedStyle(e).fontWeight}`;

/** How a ui.breadcrumbs link looks where apps put one: in the ribbon's name (the Demo app, Chem MPO). */
function breadcrumbsLinkLook(): string {
  const host = ui.div([ui.breadcrumbs(['Ribbon test app', 'Ribbon test campaign']).root], 'd4-ribbon-name');
  document.body.append(host);
  try {
    return textLook(host.querySelector<HTMLElement>('.ui-breadcrumbs-element a')!);
  } finally {
    host.remove();
  }
}

async function withAppTableView(body: (tv: DG.TableView, home: DG.View) => Promise<void>): Promise<void> {
  const home = grok.shell.newView('Ribbon test app', [ui.divText('home')]);
  const df = grok.data.demo.demog(20);
  const tv = grok.shell.addTableView(df);
  try {
    await body(tv, home);
  } finally {
    for (const v of [tv, home]) {
      if (isOpen(v))
        v.close();
    }
    grok.shell.closeTable(df);
  }
}

category('TableView: Ribbon', () => {
  test('showParentCrumb is off by default: a parent alone adds no crumb', async () => {
    await withAppTableView(async (tv, home) => {
      expect(tv.showParentCrumb, false);
      const name = viewName();
      tv.parentView = home;
      tv.name = 'Ribbon test campaign';
      await waitUntil(() => viewName() !== name && viewName()?.textContent === tv.name, 'no re-render');
      expect(parentCrumb(), null);
    });
  }, {skipReason: nodeSkip});

  test('showParentCrumb shows the parent before the name on the current view', async () => {
    await withAppTableView(async (tv, home) => {
      tv.parentView = home;
      tv.showParentCrumb = true;
      expect(tv.showParentCrumb, true);
      await waitUntil(() => parentCrumb()?.textContent === home.name, 'no parent crumb');
      expect(parentCrumb()!.nextElementSibling, viewName());
      expect(viewName()!.textContent, tv.name);

      const other = DG.View.create();
      other.name = 'Ribbon test other app';
      tv.parentView = other;
      await waitUntil(() => parentCrumb()?.textContent === other.name, 'the crumb did not follow the new parent');
      tv.parentView = home;
      await waitUntil(() => parentCrumb()?.textContent === home.name, 'the crumb did not follow the parent back');

      tv.showParentCrumb = false;
      await waitUntil(() => parentCrumb() == null, 'the crumb stayed after showParentCrumb = false');
    });
  }, {skipReason: nodeSkip});

  test('the parent crumb is a link as in ui.breadcrumbs; the name keeps the ribbon\'s look', async () => {
    await withAppTableView(async (tv, home) => {
      await waitUntil(() => viewName()?.textContent === tv.name, 'the table view is not shown');
      const nameLook = textLook(viewName()!);
      tv.parentView = home;
      tv.showParentCrumb = true;
      await waitUntil(() => parentCrumb()?.textContent === home.name, 'no parent crumb');
      const crumbText = parentCrumb()!.querySelector<HTMLElement>('a') ?? parentCrumb()!;
      expect(textLook(crumbText), breadcrumbsLinkLook(), 'the parent is not drawn as a breadcrumbs link');
      expect(getComputedStyle(crumbText).cursor, 'pointer');
      expect(textLook(viewName()!), nameLook, 'the table\'s name changed its look');
    });
  }, {skipReason: nodeSkip});

  test('a click on the parent crumb makes the parent current', async () => {
    await withAppTableView(async (tv, home) => {
      tv.parentView = home;
      tv.showParentCrumb = true;
      await waitUntil(() => parentCrumb() != null, 'no parent crumb');
      parentCrumb()!.click();
      await waitUntil(() => sameView(grok.shell.v, home), 'the parent did not become current');
      grok.shell.v = tv;
      await waitUntil(() => parentCrumb()?.textContent === home.name, 'no crumb after switching back');
    });
  }, {skipReason: nodeSkip});

  test('closing a table view makes its parent current', async () => {
    await withAppTableView(async (tv, home) => {
      tv.parentView = home;
      grok.shell.v = tv;
      tv.close();
      await waitUntil(() => sameView(grok.shell.v, home), 'the parent did not become current');
    });
  }, {skipReason: nodeSkip});

  test('a parent closed before is not brought back when its table view closes', async () => {
    await withAppTableView(async (tv, home) => {
      tv.parentView = home;
      home.close();
      grok.shell.v = tv;
      tv.close();
      expect(isOpen(home), false);
    });
  }, {skipReason: nodeSkip});

  test('other table views keep no crumb and their Save', async () => {
    await withAppTableView(async (tv, home) => {
      tv.parentView = home;
      tv.showParentCrumb = true;
      tv.showSaveButton = false;
      const plainDf = grok.data.demo.demog(10);
      const plain = grok.shell.addTableView(plainDf);
      try {
        await waitUntil(() => viewName()?.textContent === plain.name, 'the plain view is not shown');
        expect(parentCrumb(), null);
        expect(plain.showParentCrumb, false);
        expect(plain.showSaveButton, true);
        expect(isShown(saveButton(plain)), true);
      } finally {
        plain.close();
        grok.shell.closeTable(plainDf);
      }
    });
  }, {skipReason: nodeSkip});

  test('showSaveButton = false hides only the Save, through a ribbon round-trip', async () => {
    await withAppTableView(async (tv) => {
      expect(tv.showSaveButton, true);
      expect(isShown(saveButton(tv)), true);
      const panels = () => Array.from(tv.root.parentElement!.querySelectorAll('.d4-ribbon-panel[name]'))
        .map((p) => p.getAttribute('name')).join(',');
      const before = panels();

      tv.showSaveButton = false;
      expect(tv.showSaveButton, false);
      expect(isShown(saveButton(tv)), false);

      const submit = ui.bigButton('Submit', () => {});
      tv.setRibbonPanels([...tv.getRibbonPanels(), [submit]]);
      expect(isShown(saveButton(tv)), false);
      expect(panels(), before);
      expect(tv.root.parentElement!.contains(submit), true);

      tv.showSaveButton = true;
      expect(isShown(saveButton(tv)), true);
    });
  }, {skipReason: nodeSkip});
}, {owner: 'askalkin@datagrok.ai'});
