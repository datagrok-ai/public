/* P3.5 WO-1/3/4 — spec metadata as the user meets it (defaults, humanized labels, 'Add …' actions,
   choices), canvas navigation (hover, one context menu, Esc, live tab headers, control bounds) and
   the icon ribbon with New form, samples and the gallery round-trip. Each check reopens EDIT_SPEC. */
import {ok, shot} from '../local.mjs';
import {EDIT_SPEC, confirmDiscard, drag, dumpViaDialog, expandTree, openSpec, paletteItem, ribbon, row,
  selectRow, statusText, surfaceNode, waitStatus} from '../lib.mjs';

export async function fixture(page) {
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
}

async function checkSpecMetadata(page) {
  await page.locator('.u2-palette input').first().fill('breadcrumbs');
  await page.waitForTimeout(150);
  await drag(page, paletteItem(page, 'u2-breadcrumbs'), surfaceNode(page, 'saveButton'));
  // visible segments only: the control keeps its overflow-ellipsis elements in the DOM, hidden
  const crumbs = await page.evaluate(() => [...document.querySelector(
    '.u2-designer-surface [data-u2-name="breadcrumbs1"]')?.children ?? []]
    .filter((el) => el.style.display !== 'none').map((el) => el.textContent).join(' '));
  ok('metadata-and-navigation/1a/insert-seeds-the-meta-defaults', crumbs === 'Item 1 › Item 2 › Item 3',
    `"${crumbs}"`);

  await page.locator('.u2-palette input').first().fill('text');
  await page.waitForTimeout(150);
  await drag(page, paletteItem(page, 'u2-text-input'), surfaceNode(page, 'saveButton'));
  const label = await page.evaluate(() => document.querySelector(
    '.u2-designer-surface [data-u2-name="textInput1"] .u2-input-label')?.textContent.trim() ?? '');
  ok('metadata-and-navigation/1b/insert-humanizes-the-label', label === 'Text input 1', `"${label}"`);

  await row(page, 'tabs').click({button: 'right'});
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  await page.locator('.u2-menu-item').filter({hasText: /^Add tab$/}).first().click();
  await page.waitForTimeout(300);
  const tabLabels = await page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface .u2-tabs-label')].map((l) => l.textContent));
  const status = await statusText(page);
  ok('metadata-and-navigation/1c/add-tab-appends-a-titled-selected-pane',
    tabLabels.length === 3 && tabLabels[2] === 'Tab 3' && /panel1 · \d+ nodes/.test(status),
    `tabs=${JSON.stringify(tabLabels)} status="${status}"`);

  await page.locator('.u2-palette input').first().fill('splitter');
  await page.waitForTimeout(150);
  await drag(page, paletteItem(page, 'u2-splitter'), surfaceNode(page, 'saveButton'));
  await waitStatus(page, 'splitter1');
  await selectRow(page, 'splitter1');
  // the field is addressed feature-detected: `data-u2-prop` once the panel runs on propertyForm
  // (WO-2), the propgrid row until then — the assertion is the same dropdown either way
  const choices = await page.evaluate(() => {
    const stamped = document.querySelector('[data-u2-prop="direction"]');
    const select = stamped?.querySelector('select') ??
      [...document.querySelectorAll('.u2-designer-properties .u2-propgrid-row')]
        .find((r) => r.querySelector('.u2-propgrid-name')?.textContent.trim() === 'direction')
        ?.querySelector('select');
    return select ? [...select.options].map((o) => o.value || o.textContent) : null;
  });
  await shot(page, 'metadata-and-navigation-1-spec-metadata');
  ok('metadata-and-navigation/1d/choices-render-a-dropdown', choices !== null &&
    ['horizontal', 'vertical'].every((c) => choices.includes(c)), JSON.stringify(choices));
}

async function checkCanvasNavigation(page) {
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);

  await selectRow(page, 'saveButton');
  await row(page, 'nameInput').hover();
  await page.waitForTimeout(200);
  const hover = await page.evaluate(() => {
    const round = (r) => r && {x: Math.round(r.x), y: Math.round(r.y),
      w: Math.round(r.width), h: Math.round(r.height)};
    const box = document.querySelector('.u2-designer-hover');
    return {
      display: getComputedStyle(box).display,
      box: round(box.getBoundingClientRect()),
      target: round(document.querySelector('.u2-designer-surface [data-u2-name="nameInput"]')
        ?.getBoundingClientRect()),
    };
  });
  const near = !!hover.target && ['x', 'y', 'w', 'h']
    .every((k) => Math.abs(hover.box[k] - hover.target[k]) <= 3);
  ok('metadata-and-navigation/2a/tree-hover-highlights-the-canvas-node', hover.display === 'block' && near,
    JSON.stringify(hover));
  await page.mouse.move(10, 10);

  const box = await surfaceNode(page, 'saveButton').boundingBox();
  await page.mouse.click(box.x + box.width / 2, box.y + box.height / 2, {button: 'right'});
  await page.waitForSelector('.u2-menu', {timeout: 5000});
  const menus = await page.evaluate(() => ({
    u2: document.querySelectorAll('.u2-menu').length,
    platform: [...document.querySelectorAll('*')].filter((el) =>
      el.children.length === 0 && /^Properties(\.\.\.|…)$/.test(el.textContent.trim())).length,
  }));
  await shot(page, 'metadata-and-navigation-2-single-context-menu');
  ok('metadata-and-navigation/2b/one-menu-and-no-platform-properties', menus.u2 === 1 && menus.platform === 0,
    JSON.stringify(menus));
  await page.keyboard.press('Escape');
  await page.waitForTimeout(200);

  await selectRow(page, 'nameInput');
  await page.keyboard.press('Escape');
  const parentStatus = await waitStatus(page, 'firstPane');
  ok('metadata-and-navigation/2c/esc-selects-the-parent', /› firstPane · \d+ nodes/.test(parentStatus),
    parentStatus);

  const second = await page.locator('.u2-designer-surface .u2-tabs-tab').nth(1).boundingBox();
  await page.mouse.click(second.x + second.width / 2, second.y + second.height / 2);
  const tabStatus = await waitStatus(page, 'secondPane');
  const paneShown = await page.evaluate(() => Math.round(document.querySelector(
    '.u2-designer-surface [data-u2-name="secondPane"]')?.getBoundingClientRect().height ?? 0));
  ok('metadata-and-navigation/2d/tab-header-activates-and-selects-the-pane',
    paneShown > 0 && /› secondPane · \d+ nodes/.test(tabStatus),
    `pane height=${paneShown} status="${tabStatus}"`);

  await selectRow(page, 'nameInput');
  await page.waitForTimeout(200);
  const revealed = await page.evaluate(() => {
    const round = (r) => r && {w: Math.round(r.width), h: Math.round(r.height)};
    return {
      input: round(document.querySelector('.u2-designer-surface [data-u2-name="nameInput"]')
        ?.getBoundingClientRect()),
      adorner: round(document.querySelector('.u2-designer-selected').getBoundingClientRect()),
      display: getComputedStyle(document.querySelector('.u2-designer-selected')).display,
    };
  });
  await shot(page, 'metadata-and-navigation-2-reveal-on-select');
  ok('metadata-and-navigation/2e/tree-select-reveals-the-hidden-pane', revealed.display === 'block' &&
    !!revealed.input && revealed.input.h > 0 && revealed.adorner.h > 0, JSON.stringify(revealed));

  const outline = () => page.evaluate(() => getComputedStyle(document.querySelector(
    '.u2-designer-surface [data-u2-name="saveButton"]')).outlineStyle);
  const inDesign = await outline();
  await ribbon(page, 'Control bounds').click();
  await page.waitForTimeout(100);
  const toggledOff = await outline();
  await ribbon(page, 'Control bounds').click();
  await page.waitForTimeout(100);
  await page.locator('.d4-ribbon').getByText('Run', {exact: true}).first().click();
  await page.waitForTimeout(300);
  const inRun = await outline();
  await page.locator('.d4-ribbon').getByText('Design', {exact: true}).first().click();
  await page.waitForTimeout(300);
  ok('metadata-and-navigation/2f/control-bounds-design-only-and-toggleable',
    inDesign === 'dashed' && toggledOff === 'none' && inRun === 'none' && await outline() === 'dashed',
    `design=${inDesign} off=${toggledOff} run=${inRun}`);
}

/** Walks the Open menu: a group label opens its submenu, an item label runs it. The label
 * span is the target — a group item's own text grows its submenu's once that opens. */
async function menuPick(page, ...labels) {
  await page.locator('.d4-ribbon').getByRole('button', {name: 'Open', exact: true}).first().click();
  await page.waitForSelector('.u2-menu-item', {timeout: 5000});
  for (const label of labels) {
    await page.locator('.u2-menu .u2-menu-label').filter({hasText: new RegExp(`^${label}$`)}).first().click();
    await page.waitForTimeout(150);
  }
}

async function checkRibbonAndGallery(page) {
  const names = await page.evaluate(() => [...document.querySelectorAll('.d4-ribbon button')]
    .map((b) => b.getAttribute('aria-label') || b.textContent.trim()));
  ok('metadata-and-navigation/3a/every-ribbon-button-is-addressable-by-role', names.every((n) => n !== '') &&
    ['New form', 'Open', 'Copy spec', 'Undo', 'Redo', 'Control bounds'].every((n) => names.includes(n)),
  names.join('|'));

  await page.evaluate(() => localStorage.removeItem('u2.designer.gallery'));
  await openSpec(page, EDIT_SPEC);
  await waitStatus(page, 'editRoot');
  await expandTree(page);
  await selectRow(page, 'saveButton');
  await ribbon(page, 'Delete').click();
  await page.waitForTimeout(300);
  ok('metadata-and-navigation/3b/an-edit-marks-modified', /modified/.test(await statusText(page)),
    await statusText(page));

  await ribbon(page, 'New form').click();
  const title = await confirmDiscard(page);
  const blankStatus = await statusText(page);
  const blank = await page.evaluate(() =>
    !!document.querySelector('.u2-designer-surface [data-u2-name="form1"]'));
  await shot(page, 'metadata-and-navigation-3-new-form');
  ok('metadata-and-navigation/3c/new-form-confirms-then-opens-the-blank-sample',
    title === 'Discard changes?' && blank && /^form1 · 1 node$/.test(blankStatus),
    `title="${title}" status="${blankStatus}"`);

  await menuPick(page, 'Samples', 'Settings');
  await waitStatus(page, 'settingsTabs');
  const tabs = await page.evaluate(() => [...document.querySelectorAll(
    '.u2-designer-surface .u2-tabs-label')].map((l) => l.textContent));
  ok('metadata-and-navigation/3d/the-settings-sample-loads-from-the-menu',
    JSON.stringify(tabs) === JSON.stringify(['General', 'Appearance', 'Advanced']), JSON.stringify(tabs));

  await expandTree(page);
  await selectRow(page, 'opacityInput');
  await ribbon(page, 'Delete').click();
  await page.waitForTimeout(300);
  const before = await dumpViaDialog(page);
  await menuPick(page, 'Save to gallery…');
  await page.waitForSelector('.u2-dialog input', {timeout: 5000});
  await page.locator('.u2-dialog input').first().fill('My form');
  await page.locator('.u2-dialog button').filter({hasText: /^OK$/i}).first().click();
  await page.waitForTimeout(200);
  const stored = await page.evaluate(() => localStorage.getItem('u2.designer.gallery') ?? '');
  // saving marks the document clean: 'modified' clears and New form needs no confirmation
  ok('metadata-and-navigation/3e/save-to-gallery-persists-and-marks-clean', stored.includes('My form') &&
    stored.includes('settingsTabs') && !/modified/.test(await statusText(page)),
  `${stored.length} chars status="${await statusText(page)}"`);

  await ribbon(page, 'New form').click();
  await waitStatus(page, 'form1');
  await menuPick(page, 'Gallery', 'My form');
  await waitStatus(page, 'settingsTabs');
  const after = await dumpViaDialog(page);
  await shot(page, 'metadata-and-navigation-3-gallery-round-trip');
  ok('metadata-and-navigation/3f/the-gallery-entry-loads-back-identically',
    after === before && JSON.parse(after).root.name === 'settingsTabs' && !/opacityInput/.test(after),
    `before=${before.length} after=${after.length} chars`);
  await page.evaluate(() => localStorage.removeItem('u2.designer.gallery'));
}

export const checks = [
  {id: 'metadata-and-navigation/1 spec metadata: defaults, labels, designer actions, choices', run: checkSpecMetadata},
  {id: 'metadata-and-navigation/2 canvas navigation: hover, menus, Esc, live tab headers, bounds',
    run: checkCanvasNavigation},
  {id: 'metadata-and-navigation/3 icon ribbon, New form, sample and gallery round-trip', run: checkRibbonAndGallery},
];
