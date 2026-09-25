/* The PowerPack Home page: its widgets (each host is named `widget-<friendly name>`), the search box
   above them, the widget settings the page keeps on the server, and a reload of the same page — the
   one browser page a worker has. The second account the library's sharing features share with signs
   in on that same page (its token swapped in, the page reloaded) and the running account comes back
   the same way; the account that started the feature is back when the feature ends, whatever failed. */
import {Page} from '@playwright/test';
import {Given, Then, When, element, kind} from '@datagrok-libraries/bdd';
import {type ElementRef, atFeatureEnd, el, expect, gestures, locate, pollMs, viewers} from '@datagrok-libraries/bdd/runtime';

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
element('unread notifications counter', {selector: '.pp-notifications-unread-count',
  description: '"N unread" above the list of the Notifications tab of Spotlight, gone when nothing is unread'});

const SWITCH_BACK = new WeakMap<Page, {token: string; login: string}>();

async function shellSettled(page: Page): Promise<void> {
  await page.locator('[name="Browse"]').first().waitFor({timeout: 180000});
  await page.waitForFunction(() => grok.shell.v?.type === 'datagrok', null, {timeout: 60000});
  await expect.poll(() => page.evaluate(() => {
    const contents = [...grok.shell.v.root.querySelectorAll('.power-pack-widgets-host .power-pack-widget-content')] as HTMLElement[];
    if (contents.length > 0 && contents.every((c) => c.children.length > 0 && c.querySelector('.grok-loader') == null))
      return 'loaded';
    return `the "${grok.shell.v?.name}" view of ${grok.shell.user?.login} at ${location.pathname}, ${contents.length} widget(s), ` +
      `${DG.Func.find({package: 'PowerPack'}).length} PowerPack function(s) of ${new Set(DG.Func.find().map((f: any) => f.package?.name)).size} packages`;
  }), {message: 'the Home widgets, loaded after the page loaded', timeout: pollMs(60000)}).toBe('loaded');
}

/** The same page loaded again: the shell boots anew, the Home widgets load, and the in-page runtime
 * the checks read through is installed again. What the page logs while it boots stays on the error
 * floor — a widget that fails to load is the claim of a scenario that reloads. */
async function reload(page: Page, home = false, simpleMode?: boolean): Promise<void> {
  const simple = simpleMode ?? await page.evaluate(() => grok.shell.windows.simpleMode).catch(() => true);
  if (home)
    await page.goto(new URL('/', page.url()).href, {waitUntil: 'domcontentloaded', timeout: 180000});
  else
    await page.reload({waitUntil: 'domcontentloaded', timeout: 180000});
  await shellSettled(page);
  await page.evaluate((s) => {
    document.body.classList.add('selenium');
    grok.shell.windows.simpleMode = s;
  }, simple);
  await viewers.installViewerRuntime(page);
}

export const reloadPage = When('user reloads the page', (page: Page) => reload(page),
  {tier: 'ui', description: 'the browser\'s reload of the one page the feature has; done when the Home widgets have loaded again'});

/* The session is swapped the way a sign-out and a sign-in swap it: the auth cookie replaced, and what
   a sign-out clears (Auth.cleanup in user.dart) cleared — the function list the client cached for the
   account before, which the next account would otherwise start from. The cache is deleted from a
   document of the same origin the platform does not run in, so no open connection holds the deletion
   up and it has finished before the platform loads again. */
async function signInWith(page: Page, token: string): Promise<void> {
  const origin = new URL(page.url()).origin;
  const simple = await page.evaluate(() => grok.shell.windows.simpleMode).catch(() => true);
  await page.goto(`${origin}/favicon.ico`, {waitUntil: 'load', timeout: 60000});
  await page.context().clearCookies({name: 'auth'});
  await page.context().addCookies([{name: 'auth', value: token, domain: new URL(origin).hostname, path: '/'}]);
  await page.evaluate((t) => new Promise<void>((resolve, reject) => {
    localStorage.setItem('auth', t);
    const request = indexedDB.deleteDatabase('CachedFuncs');
    request.onsuccess = () => resolve();
    request.onerror = () => reject(new Error(`the client's function cache was not deleted: ${request.error}`));
  }), token);
  // the Home page is what a sign-in lands on
  await reload(page, true, simple);
}

/** The running account's session is remembered once per feature, and comes back when the feature ends. */
async function rememberOwnSession(page: Page): Promise<void> {
  if (!SWITCH_BACK.has(page)) {
    SWITCH_BACK.set(page, await page.evaluate(() => ({token: localStorage.getItem('auth') ?? String(grok.dapi.token),
      login: String(grok.shell.user.login)})));
    atFeatureEnd(page, async () => {
      const own = SWITCH_BACK.get(page);
      SWITCH_BACK.delete(page);
      if (own && await page.evaluate(() => String(grok.shell.user.login)).catch(() => '') !== own.login)
        await signInWith(page, own.token);
    });
  }
}

export const sharingUserAtHand = Given('the sharing user can sign in on this page', async (page: Page) => {
  if (!process.env.DATAGROK_SHARING_LOGIN)
    throw new Error('no second account: set DATAGROK_SHARING_LOGIN, or run with a dev key so the setup can create one');
  await rememberOwnSession(page);
}, {tier: 'api', description: 'remembers the running account\'s session, which comes back when the feature ends'});

/** A session of an account, made with the running account's rights through the account's developer key. */
async function sessionOf(page: Page, login: string): Promise<{root: string; token: string}> {
  const {root, token} = await page.evaluate(() => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: String(grok.dapi.token)}));
  const user = await (await page.request.get(`${root}/public/v1/users/${encodeURIComponent(login)}`, {headers: {Authorization: token}})).json();
  const key = await (await page.request.get(`${root}/users/${user.id}/dev_key`, {headers: {Authorization: token}})).json();
  const session = await (await page.request.post(`${root}/users/login/dev`, {headers: {Authorization: `Dev ${key}`}, data: ''})).json();
  if (!session.token)
    throw new Error(`could not sign in as "${login}": ${session.message ?? 'no token'}`);
  return {root, token: session.token};
}

export const sharingUserAllRead = Given('the sharing user has no unread notifications', async (page: Page) => {
  const {root, token} = await sessionOf(page, process.env.DATAGROK_SHARING_LOGIN!);
  const done = await page.request.post(`${root}/users/notifications/current/read`, {headers: {Authorization: token}, data: ''});
  if (!done.ok())
    throw new Error(`marking the notifications of the sharing user read: HTTP ${done.status()}`);
  await expect.poll(async () => (await page.request.get(`${root}/users/notifications/current/count_unread`, {headers: {Authorization: token}})).text(),
    {message: 'the unread notifications of the sharing user'}).toBe('0');
}, {tier: 'api', description: 'everything the second account was notified of before is marked read on the server, so what is unread afterwards is new'});

export const signInAsSharingUser = When('user signs in as the sharing user', async (page: Page) => {
  const login = process.env.DATAGROK_SHARING_LOGIN!;
  const session = await sessionOf(page, login);
  await signInWith(page, session.token);
  await expect.poll(() => page.evaluate(() => grok.shell.user.login), {message: 'the signed-in account'}).toBe(login);
}, {tier: 'api', description: 'a session of the second account through its developer key, put where the page keeps its own, and the page reloaded'});

/* An account of its own for a feature that changes what an account keeps on the server — the widget
   settings of its Home page — so that no page of another feature, signed in as the running account,
   reads or writes them meanwhile. A user cannot be deleted: the account is made once per stand, and
   kept a member of Administrators, whose Home page has every widget. */
export const adminAccountOnServer = Given('an administrator account {string} is on the server', async (page: Page, login: string) => {
  const member = await page.evaluate(async (l) => {
    let user = await grok.dapi.users.filter(`login = "${l}"`).first();
    if (!user) {
      const made = DG.User.create();
      made.login = l;
      made.email = `${l}@datagrok.ai`;
      made.firstName = l;
      made.lastName = '';
      made.status = 'active';
      await grok.dapi.users.save(made);
      user = await grok.dapi.users.filter(`login = "${l}"`).first();
    }
    const admins = (await grok.dapi.groups.list({pageSize: 1000})).find((g: any) => g.friendlyName === 'Administrators' || g.name === 'Administrators');
    const own = await grok.dapi.groups.include('memberships,adminMemberships').find(user.group.id);
    if (![...own.memberships, ...own.adminMemberships].some((g: any) => g.id === admins.id))
      await grok.dapi.groups.addMember(admins, own);
    const again = await grok.dapi.groups.include('memberships').find(user.group.id);
    return again.memberships.some((g: any) => g.id === admins.id);
  }, login);
  expect(member, `${login} a member of Administrators`).toBe(true);
}, {tier: 'api', description: 'found by login or made (once per stand: users cannot be deleted), and a member of Administrators'});

export const signedInAs = Given('user is signed in as {string} on this page', async (page: Page, login: string) => {
  await rememberOwnSession(page);
  await signInWith(page, (await sessionOf(page, login)).token);
  await expect.poll(() => page.evaluate(() => grok.shell.user.login), {message: 'the signed-in account'}).toBe(login);
}, {tier: 'api', description: 'a session of that account through its developer key, and the Home page loaded with it; the running account comes back when the feature ends'});

export const signBackIn = When('user signs back in', async (page: Page) => {
  const own = SWITCH_BACK.get(page);
  if (!own)
    throw new Error('"the sharing user can sign in on this page" did not run: no session to come back to');
  await signInWith(page, own.token);
  await expect.poll(() => page.evaluate(() => grok.shell.user.login), {message: 'the signed-in account'}).toBe(own.login);
}, {tier: 'api', description: 'the session the feature started with, and the page reloaded; the account is the one the feature started with'});

export const sharingUserNotIn = Then('the signed-in user should not be a member of {string}', async (page: Page, group: string) => {
  await expect.poll(() => page.evaluate(async (g) => {
    const own = await grok.dapi.groups.include('memberships,adminMemberships').find(grok.shell.user.group.id);
    return [...own.memberships, ...own.adminMemberships].map((m: any) => m.friendlyName ?? m.name).includes(g);
  }, group), {message: `whether ${group} is among the groups of the signed-in user`}).toBe(false);
}, {tier: 'api', description: 'the groups the server lists for the account, directly or as an admin member'});

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
    return hosts.filter((h) => h.offsetParent !== null).filter((h) => {
      const content = h.querySelector('.power-pack-widget-content') as HTMLElement | null;
      return !content || content.querySelector('.grok-loader') != null ||
        ((content.innerText ?? '').trim() === '' && content.querySelector('canvas') == null);
    }).map((h) => h.getAttribute('name'));
  }), {message: 'the widgets whose content is empty or still loading'}).toEqual([]);
}, {description: 'each visible widget has finished loading and shows text or a chart'});

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

/** The account's widget settings when the feature started, written back and the page reloaded when
 * the feature ends with any of them changed — a widget a failed scenario left hidden would otherwise
 * be hidden for every later feature of the account. */
/* Read and written through the server's user settings storage with the session of the account that
   is signed in now, so the settings come back to that account even when the page has signed in as
   another one by the time the feature ends. */
export const widgetSettingsRestored = Given('the widget settings of the Home page come back when the feature ends', async (page: Page) => {
  const {root, token, login} = await page.evaluate(() => ({root: new URL(grok.dapi.root, location.href).href.replace(/\/$/, ''),
    token: String(grok.dapi.token), login: String(grok.shell.user.login)}));
  const url = `${root}/user_settings_storage/widgets?currentUser=true`;
  const read = async (): Promise<Record<string, string>> => (await (await page.request.get(url, {headers: {Authorization: token}})).json()) ?? {};
  const before = await read();
  atFeatureEnd(page, async () => {
    const now = await read();
    const changed = [...new Set([...Object.keys(before), ...Object.keys(now)])].filter((k) => now[k] !== before[k]);
    if (changed.length === 0)
      return;
    // a widget the feature saw for the first time has no setting to go back to: it goes back to shown
    const restored = {...Object.fromEntries(Object.keys(now).map((k) => [k, JSON.stringify({ignored: false})])), ...before};
    const put = await page.request.put(url, {headers: {Authorization: token, 'Content-Type': 'application/json'}, data: restored});
    if (!put.ok())
      throw new Error(`putting back the widget settings of ${login}: HTTP ${put.status()}`);
    if (await page.evaluate(() => String(grok.shell.user.login)).catch(() => '') === login)
      await reload(page);
  });
}, {tier: 'api', description: 'the widget settings of the signed-in account read from the server now, and put back at feature end if they changed'});

export const linksPointTo = Then('every link of {element} should point to {string}', async (page: Page, target: ElementRef, prefix: string) => {
  const root = (await locate(page, target)).first();
  await expect.poll(() => root.locator('a').evaluateAll((all) => all.length), {message: `the links of ${target.phrase}`}).toBeGreaterThanOrEqual(3);
  const off = await root.locator('a').evaluateAll((all, p) => all.map((a) => (a as HTMLAnchorElement).href).filter((h) => !h.startsWith(p)), prefix);
  expect(off, `the links of ${target.phrase} that do not point to ${prefix}`).toEqual([]);
}, {description: 'at least three links, and the address of every one starts with the prefix; nothing is opened'});

export const tabShowing = Then('the {string} tab of {element} should be showing', async (page: Page, tab: string, target: ElementRef) => {
  const root = (await locate(page, target)).first();
  await expect(root.locator(`.d4-tab-header[name="${tab}"]`), `the ${tab} tab header of ${target.phrase}`).toHaveClass(/(^|\s)selected(\s|$)/);
  const content = root.locator(`[data-source="tab-content-${tab}"]`).first();
  await expect(content, `the ${tab} page of ${target.phrase}`).toBeVisible();
}, {description: 'the tab header marked selected and the page it names shown, not the one before; what the page lists is a claim of its own'});

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
  await link.click();
  // the demo app names its view after the demo, the last part of the path the link shows
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.name ?? '')),
    {message: `the view the demo of the day "${name}" opened`, timeout: pollMs(60000)}).toBe(name);
  await page.evaluate(() => grok.shell.v.close());
  await expect.poll(() => page.evaluate(() => String(grok.shell.v?.type ?? '')),
    {message: 'the current view once the demo is closed'}).toBe('datagrok');
}, {tier: 'ui', description: 'the tip changes with the weekday (a demo Monday, Friday and the weekend, a tutorial Wednesday, a plain tip Tuesday and Thursday): a tutorial has to open the Tutorials panel, a demo the view named as the demo; either is closed again and the Home page is current; a plain tip must carry text and no link'});

export const unreadNotificationsOnServer = Then('the signed-in user should have {int} unread notification(s) on the server', async (page: Page, count: number) => {
  await expect.poll(() => page.evaluate(async () => grok.dapi.users.notifications.countUnread()),
    {message: 'the unread notifications of the signed-in account on the server'}).toBe(count);
}, {tier: 'api'});

export const someUnreadNotifications = Then('the signed-in user should have unread notifications on the server', async (page: Page) => {
  await expect.poll(() => page.evaluate(async () => grok.dapi.users.notifications.countUnread()),
    {message: 'the unread notifications of the signed-in account on the server'}).toBeGreaterThan(0);
}, {tier: 'api'});

export const placeholderStarts = Then('the placeholder of {element} should start with {string}', async (page: Page, target: ElementRef, text: string) => {
  await expect((await locate(page, target)).first(), `the placeholder of ${target.phrase}`)
    .toHaveAttribute('placeholder', new RegExp(`^${text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}`));
});

/* What another account shared is listed under Shared with me by the name of the account that
   shared it: the name of the account that started the feature, taken before anyone signs in. */
const RUNNING_ACCOUNT = new WeakMap<Page, string>();

export const runningAccountKnown = Given('the name of the running account is remembered', async (page: Page) => {
  RUNNING_ACCOUNT.set(page, await page.evaluate(() => String(grok.shell.user.friendlyName)));
}, {tier: 'api', description: 'the display name the Browse tree lists the account\'s shares under'});

function sharedPhrase(page: Page, path: string): ElementRef {
  const account = RUNNING_ACCOUNT.get(page);
  if (!account)
    throw new Error('"the name of the running account is remembered" did not run');
  return el(`"My stuff > Shared with me > ${account}${path ? ` > ${path}` : ''}" tree node inside browse tree`);
}

export const expandShared = When('user expands {string} shared by the running account', async (page: Page, path: string) => {
  await gestures.setExpanded(page, sharedPhrase(page, path === '.' ? '' : path), true);
}, {tier: 'ui', description: 'a node under My stuff > Shared with me > <the running account> ("." for that account\'s own node)'});

export const openShared = When('user double-clicks on {string} shared by the running account', async (page: Page, path: string) => {
  await gestures.dblclick(page, sharedPhrase(page, path));
}, {tier: 'ui', description: 'a node under My stuff > Shared with me > <the running account>'});

export const viewOpen = Then('the {string} view should be open', async (page: Page, name: string) => {
  await expect.poll(() => page.evaluate(() => [...grok.shell.views].map((v: any) => String(v.name))),
    {message: 'the views open in the shell'}).toContain(name);
}, {description: 'among the open views, current or not'});

export const containsOneOf = Then('{element} should contain one of the texts {string}', async (page: Page, target: ElementRef, list: string) => {
  const texts = list.split('|').map((t) => t.trim()).filter(Boolean);
  const loc = (await locate(page, target)).first();
  await expect.poll(async () => {
    const shown = await loc.innerText().catch(() => '');
    return texts.some((t) => shown.includes(t)) ? 'one of them' : `none of them in: ${shown.slice(0, 120)}`;
  }, {message: `one of ${texts.join(' | ')} in ${target.phrase}`}).toBe('one of them');
}, {description: 'the texts separated by |: what an element shows in one of its known states'});
