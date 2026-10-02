/* The PowerPack Home page: its widgets (each host is named `widget-<friendly name>`), the search box
   above them, and the widget settings the page keeps on the server. Signing in as another account, a
   reload, tabs and the placeholder of a field are library vocabulary. */
import {Page} from '@playwright/test';
import {Given, Then, element, kind} from '@datagrok-libraries/bdd';
import {reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {type ElementRef, atFeatureEnd, expect, locate, pollMs} from '@datagrok-libraries/bdd/runtime';

declare const grok: any;
declare const DG: any;

kind('home widget', {
  selector: '.power-pack-widget-host',
  match: ['dart'],
  dartNames: ['widget-{q}'],
  parts: {'close icon': '.grok-font-icon-close', title: '.d4-dialog-title', content: '.power-pack-widget-content',
    badge: '.pp-notification-badge', 'tip of the day': '.power-pack-activity-widget-spotlight-tip',
    'spotlight page': '[data-source="tab-content-Spotlight"]', 'notifications page': '[data-source="tab-content-Notifications"]'},
  description: 'a widget of the PowerPack Home page by its friendly name ("Spotlight", "Community", "Usage", "Reports")',
});

element('home search', {selector: '.power-search-search-everywhere-input',
  description: 'the "Search everywhere" box at the top of the Home page'});
element('home widgets panel', {selector: '.power-pack-widgets-panel',
  description: 'the widgets and the "Customize widgets..." link under them; hidden while a search is shown'});
element('home search results', {selector: '.power-pack-search-host',
  description: 'where a search shows what it found; aria-busy until every search path has settled'});

/* An account of its own for a feature that changes what an account keeps on the server — the widget
   settings of its Home page — so that no page of another feature, signed in as the running account,
   reads or writes them meanwhile. A user cannot be deleted: the account is made once per stand. It is
   a member of Administrators, whose Home page has every widget, only while the feature runs. */
const administrators = (page: Page, login: string, member: boolean): Promise<boolean> => page.evaluate(async ([l, want]) => {
  const user = await grok.dapi.users.filter(`login = ${JSON.stringify(l)}`).first();
  const admins = (await grok.dapi.groups.list({pageSize: 1000})).find((g: any) => g.friendlyName === 'Administrators' || g.name === 'Administrators');
  const isIn = async () => (await grok.dapi.groups.include('memberships').find(user.group.id)).memberships.some((g: any) => g.id === admins.id);
  if (await isIn() !== want) {
    const own = await grok.dapi.groups.find(user.group.id);
    await (want ? grok.dapi.groups.addMember(admins, own) : grok.dapi.groups.removeMember(admins, own));
  }
  return isIn();
}, [login, member] as [string, boolean]);

export const adminAccountOnServer = Given('an administrator account {string} is on the server', async (page: Page, login: string) => {
  await page.evaluate(async (l) => {
    if (await grok.dapi.users.filter(`login = ${JSON.stringify(l)}`).first())
      return;
    const made = DG.User.create();
    made.login = l;
    made.email = `${l}@datagrok.ai`;
    made.firstName = l;
    made.lastName = '';
    made.status = 'active';
    await grok.dapi.users.save(made);
  }, login);
  atFeatureEnd(page, async () => {
    await expect.poll(() => administrators(page, login, false), {message: `${login} still a member of Administrators`,
      timeout: pollMs(30000)}).toBe(false);
  });
  await expect.poll(() => administrators(page, login, true), {message: `${login} a member of Administrators`,
    timeout: pollMs(30000)}).toBe(true);
}, {tier: 'api', description: 'found by login or made (once per stand: users cannot be deleted), and a member of Administrators until the feature ends, read back both ways'});

/** The widgets as the page lays them out: the host orders them with CSS `order`, not in DOM order. */
function shownWidgets(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const hosts = [...document.querySelectorAll('.power-pack-widgets-host > .power-pack-widget-host')] as HTMLElement[];
    return hosts.filter((h) => h.offsetParent !== null)
      .map((h) => ({name: (h.getAttribute('name') ?? '').replace(/^widget-/, ''), box: h.getBoundingClientRect()}))
      .sort((a, b) => Math.round(a.box.top) - Math.round(b.box.top) || a.box.left - b.box.left)
      .map((w) => w.name);
  });
}

export const homeWidgetsAre = Then('the Home page should show the widgets {string}', async (page: Page, list: string) => {
  await expect.poll(() => shownWidgets(page).then((w) => w.join(', ')),
    {message: 'the widgets of the Home page in reading order (top to bottom, left to right)'}).toBe(list);
}, {description: 'every visible widget, in the order the page lays them out; a widget left out or added fails'});

export const everyWidgetHasContent = Then('every widget of the Home page should show content', async (page: Page) => {
  await expect.poll(() => page.evaluate(() => {
    const hosts = [...document.querySelectorAll('.power-pack-widgets-host > .power-pack-widget-host')] as HTMLElement[];
    // Community shows what community.datagrok.ai answers, an outside service (the scope rule)
    return hosts.filter((h) => h.offsetParent !== null && h.getAttribute('name') !== 'widget-Community').filter((h) => {
      const content = h.querySelector('.power-pack-widget-content') as HTMLElement | null;
      return !content || content.querySelector('.grok-loader') != null ||
        ((content.innerText ?? '').trim() === '' && content.querySelector('canvas') == null);
    }).map((h) => h.getAttribute('name'));
  }), {message: 'the widgets whose content is empty or still loading'}).toEqual([]);
}, {description: 'each visible widget but Community (the community site\'s answer) has finished loading and shows text or a chart'});

const TITLES = {Spotlight: 'activityDashboardWidget'} as Record<string, string>;

async function storedWidgetSetting(page: Page, widget: string): Promise<string> {
  return page.evaluate(async ([title, known]) => {
    const f = DG.Func.find({meta: {role: 'dashboard'}}).find((x: any) => x.friendlyName === title) ?? DG.Func.byName(known ?? title);
    if (!f)
      return `no Home widget "${title}"`;
    const stored = await grok.dapi.userDataStorage.get('widgets', true);
    const setting = stored?.[f.name];
    return setting == null ? 'not stored' : JSON.parse(setting).ignored ? 'hidden' : 'shown';
  }, [widget, TITLES[widget]] as [string, string | undefined]);
}

export const widgetStoredAs = Then('the {word} widget should be stored as {word}', async (page: Page, widget: string, state: string) => {
  await expect.poll(() => storedWidgetSetting(page, widget), {message: `the "${widget}" entry of the widget settings on the server`,
    timeout: pollMs(15000)}).toBe(state);
}, {tier: 'api', description: '"hidden" or "shown": the entry the Home page saves for the widget, read from the server — what a reload starts from'});

/** The widget settings of the account signed in now, all shown, when the feature starts and again when
 * it ends — the account is the feature's own, so a widget a killed run left hidden does not reach the
 * next run. Read and written through the server's user settings storage with the session of that
 * account, so they are reset even when the page has signed in as another one by the time the feature
 * ends; a page of the account is reloaded after a change, since a page writes back what it read. */
export const widgetsAllShown = Given('every widget of the Home page is stored as shown, now and when the feature ends', async (page: Page) => {
  const {root, token, login} = await page.evaluate(() => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: String(grok.dapi.token), login: String(grok.shell.user.login)}));
  const url = `${root}/user_settings_storage/widgets?currentUser=true`;
  const headers = {Authorization: token};
  const read = async (): Promise<Record<string, string>> => (await (await page.request.get(url, {headers})).json()) ?? {};
  const parsed = (v: string | undefined): Record<string, unknown> => {
    try {
      return JSON.parse(v ?? '{}') ?? {};
    }
    catch {
      return {};
    }
  };
  const hidden = async (): Promise<string[]> => Object.entries(await read()).filter(([, v]) => parsed(v).ignored === true).map(([k]) => k);
  const showAll = async (): Promise<boolean> => {
    const now = await read();
    if (!Object.values(now).some((v) => parsed(v).ignored === true))
      return false;
    const shown = Object.fromEntries(Object.entries(now).map(([k, v]) => [k, JSON.stringify({...parsed(v), ignored: false})]));
    const put = await page.request.put(url, {headers: {...headers, 'Content-Type': 'application/json'}, data: shown});
    if (!put.ok())
      throw new Error(`storing the widgets of ${login} as shown: HTTP ${put.status()}`);
    await expect.poll(hidden, {message: `the widgets of ${login} stored as hidden`, timeout: pollMs(15000)}).toEqual([]);
    return true;
  };
  const reloadIfSignedIn = async (): Promise<void> => {
    if (await page.evaluate(() => String(grok.shell.user.login)).catch(() => '') === login)
      await reloadPage(page);
  };
  if (await showAll())
    await reloadIfSignedIn();
  atFeatureEnd(page, async () => {
    if (await showAll())
      await reloadIfSignedIn();
  });
}, {tier: 'api', description: 'the widget settings of the signed-in account on the server: every widget shown now and when the feature ends, read back; the page reloaded when that changed them'});

export const tipOfTheDay = Then('the tip of the day of {element} should open what it names', async (page: Page, target: ElementRef) => {
  const tip = (await locate(page, target)).first().locator('.power-pack-activity-widget-spotlight-tip');
  await expect(tip, `the tip of the day of ${target.phrase}`).toContainText(/(Demo|Tutorial|Tip) of the day/);
  const kind = /(Demo|Tutorial|Tip) of the day/.exec(await tip.innerText())![1];
  const link = tip.locator('a.ui-link');
  if (kind === 'Tip') {
    await expect(link, 'a link in a plain tip of the day').toHaveCount(0);
    await expect(tip).toContainText(/^.*Tip of the day: \S.+/s);
    return;
  }
  const name = (await link.innerText()).trim();
  expect(name, `the ${kind.toLowerCase()} the tip names`).not.toBe('');
  if (kind === 'Tutorial') {
    const panel = page.locator('.panel-base').filter({has: page.locator('.panel-titlebar-text', {hasText: /^Tutorials$/})});
    await expect(panel, 'the Tutorials panel before the click').toHaveCount(0);
    await link.click();
    await expect(panel.first(), `the Tutorials panel the tutorial of the day "${name}" opens`).toBeVisible({timeout: pollMs(60000)});
    await panel.first().locator('.panel-titlebar-button-close, [name="icon-times"], .grok-font-icon-close').first().click();
    await expect(panel, 'the Tutorials panel, closed again').toHaveCount(0);
    await expect.poll(() => page.evaluate(() => String(grok.shell.v?.type ?? '')), {message: 'the current view'}).toBe('datagrok');
    return;
  }
  // a demo runs when opened, and the day's may start a container (scope rule): the claim is the demo it names
  const demos: string[] = await page.evaluate(() => DG.Func.find({meta: {demoPath: null}})
    .map((f: any) => String(f.options?.demoPath ?? '').split('|').pop()!.trim()));
  expect(demos, `the platform's demos, among them the demo of the day "${name}"`).toContain(name);
}, {tier: 'ui', description: 'the tip changes with the weekday (a demo Monday, Friday and the weekend, a tutorial Wednesday, a plain tip Tuesday and Thursday): a tutorial has to open the Tutorials panel, which is closed again with the Home page current; a demo is not opened (it may need a container), its link has to name one of the platform\'s demos; a plain tip must carry text and no link'});
