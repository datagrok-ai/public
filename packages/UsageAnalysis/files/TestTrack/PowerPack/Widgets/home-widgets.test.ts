import { test, expect } from './fixtures';
import {
  WELCOME_VIEW,
  SEARCH_INPUT,
  WIDGETS_PANEL,
  SEARCH_HOST,
  widgetByTitle,
  widgetContent,
  customizeToggle,
  WIDGET_TITLE,
  WIDGET_HEADER_HIDDEN,
  SPOTLIGHT_TAB,
  DEMO_OF_THE_DAY,
  TAB_SELECTED,
  MARK_ALL_AS_READ,
  spotlightTab,
  widgetLink,
  communityLinks,
} from './selectors';
import { BROWSE_HEADER_HOME } from '../browse/selectors';
import {
  resetHome,
  reloadHome,
  widgetTitlesInOrder,
  readWidgetIgnored,
  withSettingsSync,
  closeWidgetViaIcon,
  openCustomizeForm,
  setWidgetVisible,
  restoreWidgetVisible,
  cleanWidgetSettings,
  watchErrors,
  expectNoErrors,
} from './helpers';
import { ensureBrowsePanelOpen } from '../browse/helpers';

test.describe.configure({ mode: 'serial' });

const EXPECTED_WIDGETS = ['Spotlight', 'Community', 'Usage', 'Reports'];
const SPOTLIGHT_TABS = ['Workspace', 'Spotlight', 'Favorites', 'Notifications', 'My Activity', 'Learn'];

test.describe('Home page Widgets (Widgets-*)', () => {
  test.beforeEach(async ({ homePage: page }) => {

    test.setTimeout(180_000);
    await resetHome(page);

    await restoreWidgetVisible(page, 'Community');
  });

  test.afterAll(async ({ _app }) => {

    await cleanWidgetSettings(_app.page);
  });

  test('Widgets-Display-01 — Home page renders all widgets @smoke', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    await expect(page.locator(WELCOME_VIEW)).toBeVisible();

    const search = page.locator(SEARCH_INPUT);
    await expect(search).toBeVisible();
    expect(await search.getAttribute('placeholder')).toMatch(/^Search everywhere/);

    for (const title of EXPECTED_WIDGETS)
      await expect(widgetByTitle(page, title), `widget "${title}" should render`).toBeVisible({ timeout: 15_000 });

    const order = await widgetTitlesInOrder(page);
    expect(order.indexOf('Spotlight')).toBeLessThan(order.indexOf('Community'));

    await expect(widgetByTitle(page, 'Community').locator(WIDGET_TITLE)).toHaveText('Community');
    await expect(widgetByTitle(page, 'Spotlight').locator(WIDGET_HEADER_HIDDEN)).toHaveCount(1);

    for (const title of EXPECTED_WIDGETS) {
      await expect
        .poll(async () => widgetContent(page, title).evaluate((el) => el.textContent?.trim().length ?? 0),
          { timeout: 15_000, message: `widget "${title}" content should become non-empty` })
        .toBeGreaterThan(0);
    }

    await expectNoErrors(page, sink);
  });

  test('Widgets-Search-01 — typing hides widgets, clearing restores them', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const panel = page.locator(WIDGETS_PANEL);
    const results = page.locator(SEARCH_HOST);
    await expect(panel).toBeVisible();

    await page.locator(SEARCH_INPUT).fill('aspirin');
    await expect(panel, 'widgets panel hidden while searching').toBeHidden({ timeout: 10_000 });
    await expect(results, 'search results host shown').toBeVisible();
    expect(page.url()).toContain('search');

    await page.locator(SEARCH_INPUT).fill('');
    await expect(panel, 'widgets panel restored after clearing').toBeVisible({ timeout: 10_000 });

    await expectNoErrors(page, sink);
  });

  test('Widgets-Nav-01 — Home icon returns to Home, widget set stable across reload @smoke', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const before = await widgetTitlesInOrder(page);

    await page.evaluate(async () => {
      const grok = (window as any).grok;
      const df = await grok.dapi.files.readCsv('System:DemoFiles/demog.csv');
      grok.shell.addTableView(df);
    });
    await expect(page.locator(WELCOME_VIEW)).toBeHidden({ timeout: 10_000 });

    await ensureBrowsePanelOpen(page);
    await page.locator(BROWSE_HEADER_HOME).click();
    await expect(page.locator(WELCOME_VIEW)).toBeVisible({ timeout: 10_000 });

    await reloadHome(page);
    const after = await widgetTitlesInOrder(page);
    expect(after, 'widget set should be identical before/after reload').toEqual(before);

    await expectNoErrors(page, sink);
  });

  test('Widgets-Spotlight-01 — Spotlight renders with tabs and Demo of the day @smoke', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const spotlight = widgetByTitle(page, 'Spotlight');
    await expect(spotlight).toBeVisible();

    const tabTexts = (await page.locator(SPOTLIGHT_TAB).allInnerTexts())
      .map((t) => t.replace(/\d+$/, '').trim());
    for (const tab of SPOTLIGHT_TABS)
      expect(tabTexts, `Spotlight should have a "${tab}" tab`).toContain(tab);

    await expect(spotlight.getByText(DEMO_OF_THE_DAY, { exact: false })).toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Widgets-Spotlight-02 — switching tabs (incl. Learn sub-tabs) works', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const content = widgetContent(page, 'Spotlight');

    for (const name of ['Workspace', 'Spotlight', 'Favorites', 'My Activity']) {
      const tab = spotlightTab(page, name);
      await tab.click();
      await expect(tab, `tab "${name}" should be selected`).toHaveClass(TAB_SELECTED);
      expect((await content.evaluate((el) => el.textContent?.trim().length ?? 0))).toBeGreaterThan(0);
    }

    await spotlightTab(page, 'Learn').click();
    for (const sub of ['VIDEO', 'WIKI', 'DEMO', 'TUTORIALS']) {
      const subTab = spotlightTab(page, sub);
      await expect(subTab, `Learn sub-tab "${sub}" should be present`).toBeVisible();
      await subTab.click();
    }

    await expectNoErrors(page, sink);
  });

  test('Widgets-Spotlight-03 — Notifications tab and Mark all as read', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    await spotlightTab(page, 'Notifications').click();
    const content = widgetContent(page, 'Spotlight');
    await expect(content).toBeVisible();

    const markAllAsRead = content.getByText(MARK_ALL_AS_READ, { exact: false });
    if (await markAllAsRead.isVisible().catch(() => false)) {
      await markAllAsRead.click();
      await expect(markAllAsRead, 'unread controls should clear after Mark all as read')
        .toBeHidden({ timeout: 10_000 });
    }

    await expectNoErrors(page, sink);
  });

  test('Widgets-Community-01 — links render and point to the community forum @smoke', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const links = communityLinks(page);

    await expect(links.first()).toBeVisible({ timeout: 15_000 });
    const count = await links.count();
    expect(count, 'Community should list several links').toBeGreaterThan(2);

    const hrefs = await links.evaluateAll((els) => els.map((a) => (a as HTMLAnchorElement).href));
    for (const href of hrefs)
      expect(href, 'every link should point to community.datagrok.ai').toContain('community.datagrok.ai');

    await expectNoErrors(page, sink);
  });

  test('Widgets-Usage-01 — Usage renders and Open Usage Analysis navigates', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const usage = widgetByTitle(page, 'Usage');
    await expect(usage).toBeVisible();

    await expect(usage.locator('canvas').first()).toBeVisible({ timeout: 15_000 });
    const openLink = widgetLink(page, 'Usage', 'Open Usage Analysis');
    await expect(openLink).toBeVisible();

    await expectNoErrors(page, sink);

    await openLink.dispatchEvent('click');
    await expect
      .poll(() => page.evaluate(() => (window as any).grok?.shell?.v?.name),
        { timeout: 25_000, message: 'clicking Open Usage Analysis should open the Usage Analysis app' })
      .not.toBe('Home');
  });

  test('Widgets-Reports-01 — Reports renders and Open Reports navigates', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const reports = widgetByTitle(page, 'Reports');
    await expect(reports).toBeVisible();
    const openLink = widgetLink(page, 'Reports', 'Open Reports');
    await expect(openLink).toBeVisible();

    await expectNoErrors(page, sink);

    await openLink.dispatchEvent('click');
    await expect
      .poll(() => page.evaluate(() => (window as any).grok?.shell?.v?.name),
        { timeout: 25_000, message: 'clicking Open Reports should open the reports list' })
      .not.toBe('Home');
  });

  test('Widgets-Perm-02 — admin sees Usage and Reports', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    for (const title of EXPECTED_WIDGETS)
      await expect(widgetByTitle(page, title), `admin should see "${title}"`).toBeVisible({ timeout: 15_000 });

    await expectNoErrors(page, sink);
  });

  test('Widgets-Controls-01 — hover reveals close icon, click removes widget', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    const community = widgetByTitle(page, 'Community');
    await expect(community).toBeVisible();

    const closeIcon = community.locator('i.grok-font-icon-close');
    await expect(closeIcon).toBeHidden();
    await community.hover();
    await expect(closeIcon).toBeVisible();

    await closeWidgetViaIcon(page, 'Community');
    await expect(community).toHaveCount(0);
    await expect(widgetByTitle(page, 'Spotlight')).toBeVisible();

    await expectNoErrors(page, sink);
  });

  test('Widgets-Controls-02 — closed widget stays hidden after reload @regression', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    await withSettingsSync(page, () => closeWidgetViaIcon(page, 'Community'));

    await reloadHome(page);

    await expect(widgetByTitle(page, 'Community'), 'Community must stay hidden after reload').toHaveCount(0);
    expect(await readWidgetIgnored(page, 'Community')).toBe(true);

    await expectNoErrors(page, sink);
  });

  test('Widgets-Customize-01 — Customize form hides and shows a widget', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    await openCustomizeForm(page);

    for (const title of EXPECTED_WIDGETS)
      await expect(customizeToggle(page, title), `toggle "${title}" should exist`).toBeVisible();

    await setWidgetVisible(page, 'Community', false);
    await expect(widgetByTitle(page, 'Community')).toHaveCount(0);
    expect(await readWidgetIgnored(page, 'Community')).toBe(true);

    await setWidgetVisible(page, 'Community', true);
    await expect(widgetByTitle(page, 'Community')).toBeVisible();
    expect(await readWidgetIgnored(page, 'Community')).toBe(false);

    await expectNoErrors(page, sink);
  });

  test('Widgets-Customize-02 — visibility setting persists after reload @regression', async ({ homePage: page }) => {
    const sink = watchErrors(page);

    await openCustomizeForm(page);
    await withSettingsSync(page, () => setWidgetVisible(page, 'Community', false));

    await reloadHome(page);
    await expect(widgetByTitle(page, 'Community'), 'Community stays hidden after reload').toHaveCount(0);

    await openCustomizeForm(page);
    await expect(customizeToggle(page, 'Community')).not.toBeChecked();

    await expectNoErrors(page, sink);
  });
});
