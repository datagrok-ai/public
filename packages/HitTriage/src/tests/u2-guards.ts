import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {category, expect, test, waitUntil} from '@datagrok-libraries/test/src/test';
import {byKind, byName, expectOperable, state, u2test, user} from '@datagrok-libraries/u2/testing';
import {PermissionsEditor} from '../app/shared/permissions-editor';
import {openTemplateEditor} from '../app/shared/template-editor';
import {CampaignGrouping, HDCampaignsGroupingLSKey} from '../app/consts';
import type {HitDesignApp} from '../app/hit-design-app';
import type {HitDesignTemplate, TriagePermissions} from '../app/types';
import {fieldText} from '../app/shared/campaign-catalog';
import {allUsers, designCampaign, designCsv, expectSame, formDate, noComputeTemplate, Sandbox} from './fixtures';

function openDialog(name: string): HTMLElement | undefined {
  try {
    return byKind(document, 'dialog').find((d) => state(d).name === name);
  }
  catch {
    return undefined;
  }
}

async function dialogNamed(name: string): Promise<HTMLElement> {
  await waitUntil(() => openDialog(name) !== undefined, `the ${name} dialog did not open`, 15000);
  return openDialog(name)!;
}

async function designApp(sb: Sandbox): Promise<HitDesignApp> {
  const t = noComputeTemplate('ZZB3 Guards', 'ZZB3');
  await sb.writeTemplate(t);
  await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t), designCsv);
  const app = await sb.hitDesign();
  await sb.openDesign(app, 'ZZB3-1');
  return app;
}

/** The Submit dialog of a sandbox Hit Design campaign, opened as its ribbon button opens it. */
async function submitDialog(sb: Sandbox): Promise<HTMLElement> {
  const app = await designApp(sb);
  app.submitView!.render();
  return await dialogNamed('submit');
}

const openTracker = (app: HitDesignApp) =>
  (app as unknown as {_session: {openTracker(): unknown}})._session.openTracker();

const HOVER_MARK = 'data-zzb3-hover';

// no script sets :hover: each :hover rule is copied right after itself onto a marker attribute (same specificity, same
// place in the cascade), so the element carrying the marker is styled as under the pointer
function hovered<T>(el: Element, read: () => T): T {
  const added: [CSSStyleSheet | CSSGroupingRule, CSSRule][] = [];
  const copy = (owner: CSSStyleSheet | CSSGroupingRule) => {
    for (let i = owner.cssRules.length - 1; i >= 0; i--) {
      const rule = owner.cssRules[i];
      if (rule instanceof CSSGroupingRule)
        copy(rule);
      else if (rule instanceof CSSStyleRule && rule.selectorText.includes(':hover')) {
        owner.insertRule(`${rule.selectorText.replace(/:hover/g, `[${HOVER_MARK}]`)} {${rule.style.cssText}}`, i + 1);
        added.push([owner, owner.cssRules[i + 1]]);
      }
    }
  };
  try {
    for (const sheet of Array.from(document.styleSheets)) {
      try {
        copy(sheet);
      }
      catch {
        // another origin's sheet keeps its rules to itself
      }
    }
    el.setAttribute(HOVER_MARK, '');
    return read();
  } finally {
    el.removeAttribute(HOVER_MARK);
    for (const [owner, rule] of added)
      owner.deleteRule(Array.from(owner.cssRules).indexOf(rule));
  }
}

type Rgba = [number, number, number, number];

function rgba(color: string): Rgba {
  const srgb = /color\(srgb ([\d.]+) ([\d.]+) ([\d.]+)(?: \/ ([\d.]+))?\)/.exec(color);
  if (srgb)
    return [+srgb[1] * 255, +srgb[2] * 255, +srgb[3] * 255, srgb[4] === undefined ? 1 : +srgb[4]];
  const [r, g, b, a] = (/rgba?\(([^)]*)\)/.exec(color)?.[1] ?? '0, 0, 0, 0').split(/[\s,/]+/).map(Number);
  return [r, g, b, a ?? 1];
}

const over = (top: Rgba, under: Rgba): Rgba =>
  [0, 1, 2].map((i) => top[i] * top[3] + under[i] * (1 - top[3])).concat(1) as Rgba;

function painted(el: Element): string {
  let color: Rgba = [255, 255, 255, 1];
  const chain: Element[] = [];
  for (let at: Element | null = el; at !== null; at = at.parentElement)
    chain.unshift(at);
  for (const at of chain) {
    const style = getComputedStyle(at);
    color = over(rgba(style.backgroundColor), color);
    const layers = Array.from(style.backgroundImage.matchAll(/gradient\((?:[^,(]*,\s*)?((?:rgba?|color)\([^)]*\))/g));
    for (const layer of layers.reverse())
      color = over(rgba(layer[1]), color);
  }
  return `rgb(${color.slice(0, 3).map(Math.round).join(', ')})`;
}

const textLook = (e: Element): string => `${getComputedStyle(e).color} ${getComputedStyle(e).fontWeight}`;

function breadcrumbsLinkLook(): string {
  const host = ui.div([ui.breadcrumbs(['ZZB3 app', 'ZZB3 campaign']).root], 'd4-ribbon-name');
  document.body.append(host);
  try {
    return textLook(host.querySelector('.ui-breadcrumbs-element a')!);
  } finally {
    host.remove();
  }
}

const isOpen = (v: DG.ViewBase) => Array.from(grok.shell.views).some((x) => x.dart === v.dart);

/** What stopped `gesture`, or '' when nothing did. */
function refusal(gesture: () => void): string {
  try {
    gesture();
    return '';
  }
  catch (e) {
    return (e as Error).message;
  }
}

async function closeWithCancel(dialog: HTMLElement): Promise<void> {
  user.click(byName(dialog, 'Cancel'));
  await waitUntil(() => !dialog.isConnected, 'CANCEL left the dialog open');
}

category('u2 guards', () => {
  test('Submit: every control the dialog offers can be operated, the read-only rows\' pencils and download included',
    u2test(async () => {
      const sb = await Sandbox.create('guards');
      try {
        const dialog = await submitDialog(sb);
        expectOperable(dialog);
        await closeWithCancel(dialog);
      } finally {
        await sb.dispose();
      }
    }), {timeout: 120000});

  test('Submit: the Name and File Path pencils edit in place, Escape gives the row back and keeps the dialog',
    u2test(async () => {
      const sb = await Sandbox.create('guards');
      try {
        const dialog = await submitDialog(sb);
        const name = byName(dialog, 'Name');
        expect(refusal(() => user.type(name, 'x')).includes('is read-only'), true,
          'Name takes no keys before its pencil');
        user.click(byName(dialog, 'Edit Campaign Name'));
        expect(state(name).readOnly, false, 'Name after its pencil');
        user.type(name, 'ZZB3 renamed', {replace: true});
        expect(state(name).value, 'ZZB3 renamed', 'Name while editing');
        user.press(name, 'Escape');
        expect(state(name).readOnly, true, 'Name after Escape');
        expect(state(name).value, 'Friendly ZZB3-1', 'Escape keeps the saved name');
        expect(dialog.isConnected, true, 'Escape in the row leaves the dialog open');
        const path = byName(dialog, 'File Path');
        user.click(byName(dialog, 'Edit file path'));
        expect(state(path).readOnly, false, 'File Path after its pencil');
        user.press(path, 'Escape');
        expect(state(path).readOnly, true, 'File Path after Escape');
        await closeWithCancel(dialog);
      } finally {
        await sb.dispose();
      }
    }), {timeout: 120000});

  test('Submit: a Date field reads as the campaigns table shows it', u2test(async () => {
    const sb = await Sandbox.create('guards');
    try {
      const dialog = await submitDialog(sb);
      expect(state(byName(dialog, 'Deadline')).value, fieldText(formDate), 'the Deadline row');
      await closeWithCancel(dialog);
    } finally {
      await sb.dispose();
    }
  }), {timeout: 120000});

  test('Edit Permissions: every row control works; Can Edit checks Can View; OK saves both', u2test(async () => {
    const saved: TriagePermissions[] = [];
    void PermissionsEditor.show({view: [allUsers], edit: []}, {subject: 'ZZB3-1', save: (p) => saved.push(p)});
    const dialog = await dialogNamed('edit-permissions');
    expectOperable(dialog);
    const [grid] = byKind(dialog, 'access-grid');
    const view = byName(grid, 'Can View for All users') as HTMLInputElement;
    const edit = byName(grid, 'Can Edit for All users') as HTMLInputElement;
    expectSame([view.checked, edit.checked], [true, false], 'Can View, Can Edit as stored');
    user.click(edit);
    expectSame([view.checked, edit.checked], [true, true], 'Can View, Can Edit after checking Can Edit');
    user.click(byName(dialog, 'OK'));
    await waitUntil(() => !dialog.isConnected, 'OK left Edit Permissions open', 10000);
    expectSame(saved, [{edit: [allUsers], view: [allUsers]}], 'what OK saved');
  }), {timeout: 60000});

  test('Edit Permissions of a campaign view closes with the view when another campaign opens',
    u2test(async () => {
      const sb = await Sandbox.create('guards');
      try {
        const t = noComputeTemplate('ZZB3 Guards', 'ZZB3');
        await sb.writeTemplate(t);
        await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t), designCsv);
        await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', t), designCsv);
        const app = await sb.hitDesign();
        await sb.openDesign(app, 'ZZB3-1');
        const icon = document.querySelector<HTMLElement>('[data-u2-item="Edit campaign permissions"]');
        expect(icon !== null, true, 'the ribbon\'s Edit campaign permissions');
        icon!.click();
        const dialog = await dialogNamed('edit-permissions');
        await sb.openDesign(app, 'ZZB3-2');
        await waitUntil(() => !dialog.isConnected, 'Edit Permissions of ZZB3-1 stayed open over ZZB3-2', 10000);
      } finally {
        await sb.dispose();
      }
    }), {timeout: 120000});

  test('Submit: Enter in the edited Name saves it and keeps the dialog; Download gives the campaign\'s table',
    u2test(async () => {
      const sb = await Sandbox.create('guards');
      const click = HTMLAnchorElement.prototype.click;
      const downloads: string[] = [];
      HTMLAnchorElement.prototype.click = function(this: HTMLAnchorElement) {
        if (this.hasAttribute('download'))
          downloads.push(this.download);
        else
          click.call(this);
      };
      try {
        const dialog = await submitDialog(sb);
        const name = byName(dialog, 'Name');
        user.click(byName(dialog, 'Edit Campaign Name'));
        user.type(name, 'ZZB3 renamed', {replace: true});
        user.press(name, 'Enter');
        await waitUntil(() => !dialog.isConnected || state(name).readOnly === true, 'Enter left the Name editing',
          30000);
        expect(dialog.isConnected, true, 'Enter in the row leaves the dialog open');
        await sb.savesDone();
        expect((await sb.readCampaign('ZZB3-1')).friendlyName, 'ZZB3 renamed', 'the name Enter saved');
        user.click(byName(dialog, 'Download'));
        expectSame(downloads, ['ZZB3-1.csv'], 'what Download saved');
        await closeWithCancel(dialog);
      } finally {
        HTMLAnchorElement.prototype.click = click;
        await sb.dispose();
      }
    }), {timeout: 120000});

  test('Edit Permissions: each row starts with its user\'s picture or its group\'s icon', u2test(async () => {
    const me = DG.User.current().group;
    void PermissionsEditor.show({view: [allUsers, me.id], edit: [me.id]}, {subject: 'ZZB3-1', save: () => true});
    const dialog = await dialogNamed('edit-permissions');
    const [grid] = byKind(dialog, 'access-grid');
    const pictureOf = (principal: string): Element | null => {
      const cell = Array.from(grid.querySelectorAll<HTMLTableRowElement>('tbody tr'))
        .find((r) => r.dataset.u2Key === principal)?.cells[0];
      const first = cell?.firstElementChild ?? null;
      return first !== null && first !== cell!.lastElementChild ? first : null;
    };
    const icons = (el: Element | null) => {
      try {
        return el === null ? 0 : byKind(el, 'icon').length;
      }
      catch {
        return 0;
      }
    };
    const group = pictureOf(allUsers);
    expect(group !== null && icons(group) === 1, true, 'the All users row starts with the group icon');
    const person = pictureOf(me.id);
    expect(person !== null && person.firstElementChild !== null && icons(person) === 0, true,
      `the ${me.friendlyName} row starts with the user's picture or initials`);
    await closeWithCancel(dialog);
  }), {timeout: 60000});

  test('Progress tracker: Escape and the close button keep the moves, and save them', u2test(async () => {
    const sb = await Sandbox.create('guards');
    try {
      const app = await designApp(sb);
      const df = app.dataFrame!;
      for (const [row, how] of [[0, 'Escape'], [1, 'the close button']] as const) {
        openTracker(app);
        const tracker = await dialogNamed('progress-tracker');
        const version = app.campaign!.version ?? 0;
        df.set('Stage', row, 'Tested');
        if (how === 'Escape')
          user.press(tracker, 'Escape');
        else
          user.click(Array.from(tracker.querySelectorAll('button[aria-label="Close"]')).find((b) => !b.textContent)!);
        await waitUntil(() => !tracker.isConnected, `${how} left the tracker open`, 10000);
        await waitUntil(() => app.campaign?.version === version + 1 && !app.isSaving, `${how}: the move was not saved`,
          30000);
        expect(df.get('Stage', row), 'Tested', `${how}: the move stays in the table`);
      }
      const saved = await grok.dapi.files.readCsv(sb.savePath('ZZB3-1'));
      expectSame(saved.col('Stage')!.toList(), ['Tested', 'Tested'], 'the saved stages');
    } finally {
      await sb.dispose();
    }
  }), {timeout: 120000});

  test('Template editor: a stage is free text: one typed in a new row is the template\'s stage', u2test(async () => {
    const sb = await Sandbox.create('guards');
    let created: HitDesignTemplate | undefined;
    try {
      const app = await sb.hitDesign();
      await openTemplateEditor<HitDesignTemplate>({app, onCreated: (t) => created = t});
      await waitUntil(() => grok.shell.v?.name === 'New Template', 'the template editor did not open', 15000);
      const editor = byName(grok.shell.v.root, 'template-editor');
      user.type(byName(editor, 'name'), 'ZZB3 Stages');
      user.type(byName(editor, 'key'), 'ZZBS');
      const stages = byName(editor, 'stages');
      user.click(byName(stages, 'Add stage'));
      user.type(byName(stages, 'Name'), 'ZZB3 assay review');
      const create = byName(editor, 'Create');
      await waitUntil(() => state(create).enabled, () => `CREATE stayed off: ${state(create).disabledReason}`, 30000);
      user.click(create);
      await waitUntil(() => created !== undefined, 'CREATE made no template', 30000);
      expectSame(created!.stages, ['ZZB3 assay review'], 'the stages the template keeps');
    } finally {
      if (grok.shell.v?.name === 'New Template')
        grok.shell.v.close();
      await sb.dispose();
    }
  }), {timeout: 120000});

  test('Campaigns table: the hovered row is marked, and not as a section header is', u2test(async () => {
    const sb = await Sandbox.create('guards');
    const grouping = localStorage.getItem(HDCampaignsGroupingLSKey);
    try {
      localStorage.setItem(HDCampaignsGroupingLSKey, CampaignGrouping.Status);
      const t = noComputeTemplate('ZZB3 Guards', 'ZZB3');
      await sb.writeTemplate(t);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-1', t), designCsv);
      await sb.writeCampaign(designCampaign(sb, 'ZZB3-2', t, {status: 'Submitted'}), designCsv);
      const app = await sb.hitDesign();
      const table = byName(app.infoView.root, 'campaigns');
      const row = () => table.querySelector('tr[role="row"]');
      const section = () => table.querySelector('th[role="rowheader"]');
      await waitUntil(() => row() !== null && section() !== null,
        'the campaigns table shows no section with a campaign', 30000);
      const resting = painted(row()!);
      const header = painted(section()!);
      const under = hovered(row()!, () => painted(row()!));
      expect(under !== resting, true, `the row under the pointer looks as it rests: ${under}`);
      expect(under !== header, true, `the row under the pointer looks like a section header: ${under}`);
    } finally {
      grouping === null ? localStorage.removeItem(HDCampaignsGroupingLSKey) :
        localStorage.setItem(HDCampaignsGroupingLSKey, grouping);
      await sb.dispose();
    }
  }), {timeout: 120000});

  test('Campaign view: the crumb is a breadcrumbs link; a press on it closes the campaign', u2test(async () => {
    const sb = await Sandbox.create('guards');
    try {
      const app = await designApp(sb);
      const view = app.designView!;
      const crumb = () => document.querySelector('.d4-view-ribbon > span.d4-ribbon-name');
      await waitUntil(() => crumb()?.textContent === app.mainView.name,
        'the campaign view shows no crumb to the landing', 10000);
      // the chrome hears a press on the crumb from the frame after it placed its group in the ribbon
      await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
      const link = crumb()!.querySelector('a');
      expect(link !== null && textLook(link) === breadcrumbsLinkLook(), true,
        `the crumb is not drawn as a ui.breadcrumbs link: ${textLook(link ?? crumb()!)} vs ${breadcrumbsLinkLook()}`);
      user.click(link!);
      await waitUntil(() => !isOpen(view), 'the press on the crumb left the campaign open', 10000);
      expect(isOpen(app.mainView), true, 'the landing stays open');
    } finally {
      await sb.dispose();
    }
  }), {timeout: 120000});
});
