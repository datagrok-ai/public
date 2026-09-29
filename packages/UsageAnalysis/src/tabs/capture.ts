import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './ua';
import {TimelineView} from './timeline';
import {UaToolbox} from '../ua-toolbox';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {onRowContextMenu} from '../utils';
import '../../css/usage_analysis.css';

const CAPTURE_ITEMS = ['clicks', 'inputs', 'requests', 'calls', 'errors'];
const SUBJECTS = ['user', 'group', 'package', 'everyone'];
const ALL_ACTIVITY = 'all activity';
const SCOPES = [ALL_ACTIVITY, 'view', 'element', 'function', 'error'];
const NO_LEVEL = 'none';
const LEVELS = [NO_LEVEL, 'error', 'warning', 'info', 'debug'];
const DURATIONS: {[name: string]: number} = {'30 min': 30, '2 h': 120, '1 d': 1440, '2 d': 2880, '7 d': 10080};
const HIDDEN_COLUMNS = ['status', 'capture', 'anonymous', 'name', 'max_events', 'window_minutes', 'max_sessions',
  'created_at', 'expires_at', 'ended_at', 'id'];

/** Capture rules (`capture_rules`): who is captured, by whom and why; new rules and stops go through the server. */
export class CaptureView extends UaView {
  rulesViewer?: UaFilterableQueryViewer;
  /** The rule to show again in the context panel once the rules reload. */
  reshow?: string;

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Capture';
  }

  async initViewers(path?: string): Promise<void> {
    this.rulesViewer = new UaFilterableQueryViewer({
      filterSubscription: this.uaToolbox.filterStream,
      name: 'Capture rules',
      queryName: 'CaptureRules',
      processDataFrame: (t: DG.DataFrame) => {
        t.onCurrentRowChanged.subscribe(() => this.showRule(t, t.currentRowIdx));
        const i = this.reshow ? t.col('rule')!.toList().indexOf(this.reshow) : -1;
        this.reshow = undefined;
        if (i >= 0)
          t.currentRowIdx = i;
        return t;
      },
      createViewer: (t: DG.DataFrame) => {
        const grid = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
        grid.columns.setOrder(['rule', 'author', 'subject', 'scope', 'reason', 'active', 'events']);
        for (const name of HIDDEN_COLUMNS)
          grid.col(name)!.visible = false;
        grid.col('reason')!.width = 250;
        onRowContextMenu(grid, (menu, i) => {
          if (t.get('status', i) === 'active')
            menu.item('Stop...', () => this.stopDialog(t.get('rule', i)));
          menu.item('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'rule', t.get('rule', i)));
        });
        return grid;
      },
    });
    this.viewers.push(this.rulesViewer);
    const toolbar = ui.divH([ui.button('New rule...', () => this.newRuleDialog())], 'ua-toolbar');
    this.root.append(ui.divV([toolbar, ui.box(this.rulesViewer.root)], 'ui-box'));
  }

  showRule(t: DG.DataFrame, i: number): void {
    if (i < 0)
      return;
    const rule: string = t.get('rule', i);
    const active = t.get('status', i) === 'active';
    const buttons = [ui.button('Timeline', () => TimelineView.open(this.uaToolbox.viewHandler, 'rule', rule))];
    if (active)
      buttons.unshift(ui.button('Stop...', () => this.stopDialog(rule)));
    const acc = DG.Accordion.create();
    acc.addPane(rule, () => ui.divV([
      ui.tableFromMap({
        'Name': t.get('name', i) ?? '',
        'Author': t.get('author', i) ?? '',
        'Subject': `${t.get('subject', i)}${t.get('anonymous', i) ? ' (anonymous)' : ''}`,
        'Scope': t.get('scope', i),
        'Capture': t.get('capture', i),
        'Reason': t.get('reason', i),
        'Status': t.get('status', i),
        'Created': t.col('created_at')!.getString(i),
        [active ? 'Expires' : 'Ended']: t.col(active ? 'expires_at' : 'ended_at')!.getString(i) ||
          t.col('expires_at')!.getString(i),
        'Events': `${t.get('events', i)} of ${t.get('max_events', i)}`,
        'Window': `${t.get('window_minutes', i)} min, at most ${t.get('max_sessions', i)} sessions`,
      }),
      ui.buttonsInput(buttons),
    ]), true);
    grok.shell.o = acc.root;
  }

  stopDialog(rule: string): void {
    const reason = ui.input.string('Reason', {tooltipText: 'Why the rule stops; required'});
    const dialog = ui.dialog(`Stop ${rule}`)
      .add(reason)
      .onOK(async () => {
        try {
          await grok.functions.call('CaptureRuleStop', {id: rule, reason: reason.value.trim()});
          grok.shell.info(`Stopped ${rule}`);
          this.reshow = rule;
          this.rulesViewer?.reloadViewer();
        }
        catch (e: any) {
          grok.shell.error(`${rule}: ${e?.message ?? e}`);
        }
      })
      .show();
    const ok = dialog.getButton('OK');
    ok.disabled = true;
    ui.tooltip.bind(ok, () => reason.value?.trim() ? 'Stop the rule' : 'Enter the reason');
    reason.onChanged.subscribe(() => ok.disabled = !reason.value?.trim());
  }

  newRuleDialog(): void {
    const subject = ui.input.choice('Subject', {value: 'user', items: SUBJECTS, nullable: false});
    const who = ui.input.string('Who', {tooltipText: 'The login, the group name or the package name'});
    const scope = ui.input.choice('Scope', {value: ALL_ACTIVITY, items: SCOPES, nullable: false});
    const scopeValue = ui.input.string('Scope value',
      {tooltipText: 'A view name, a part of an element path, a function nqName or an error signature'});
    const capture = ui.input.multiChoice('Capture', {value: ['clicks', 'requests', 'errors'], items: CAPTURE_ITEMS});
    const level = ui.input.choice('Server level', {value: NO_LEVEL, items: LEVELS, nullable: false});
    const flags = ui.input.string('Debug flags', {tooltipText: 'Comma-separated server debug flags, e.g. query,storage'});
    const duration = ui.input.choice('For', {value: '1 d', items: Object.keys(DURATIONS), nullable: false});
    const maxEvents = ui.input.int('Max events', {value: 10000});
    const anonymous = ui.input.bool('Anonymous', {tooltipText: 'Group and everyone rules only: no user, session or IP'});
    const name = ui.input.string('Name');
    const reason = ui.input.string('Reason', {tooltipText: 'Why this activity is captured; required'});

    const problem = (): string | null => {
      if (subject.value !== 'everyone' && !who.value?.trim())
        return `Enter the ${subject.value === 'user' ? 'login' : `${subject.value} name`}`;
      if (scope.value !== ALL_ACTIVITY && !scopeValue.value?.trim())
        return `Enter the ${scope.value}`;
      if (subject.value === 'everyone' && scope.value === ALL_ACTIVITY)
        return 'A rule for everyone needs a scope';
      if (!capture.value?.length && level.value === NO_LEVEL && !flags.value?.trim())
        return 'Choose what to capture';
      if (!reason.value?.trim())
        return 'Enter the reason';
      return null;
    };
    const refresh = () => {
      who.root.hidden = subject.value === 'everyone';
      scopeValue.root.hidden = scope.value === ALL_ACTIVITY;
      anonymous.enabled = subject.value === 'group' || subject.value === 'everyone';
      if (!anonymous.enabled && anonymous.value)
        anonymous.value = false;
      dialog.getButton('OK').disabled = problem() != null;
    };

    const dialog = ui.dialog('New capture rule');
    const inputs: DG.InputBase[] = [subject, who, scope, scopeValue, capture, level, flags, duration, maxEvents,
      anonymous, name, reason];
    for (const input of inputs) {
      dialog.add(input);
      input.onChanged.subscribe(() => refresh());
    }
    dialog.onOK(async () => {
      const body = {
        name: name.value?.trim() || undefined,
        subject: {type: subject.value, value: subject.value === 'everyone' ? undefined : who.value.trim()},
        scope: scope.value === ALL_ACTIVITY ? undefined : {type: scope.value, value: scopeValue.value.trim()},
        capture: {
          ...Object.fromEntries(CAPTURE_ITEMS.map((item) => [item, capture.value?.includes(item) ?? false])),
          serverLevel: level.value === NO_LEVEL ? null : level.value,
          debugFlags: (flags.value ?? '').split(',').map((f) => f.trim()).filter((f) => f),
        },
        anonymous: anonymous.value,
        maxEvents: maxEvents.value,
        forMinutes: DURATIONS[duration.value!],
        reason: reason.value.trim(),
      };
      try {
        const rule = JSON.parse(await grok.functions.call('CaptureRuleAdd', {rule: JSON.stringify(body)}));
        grok.shell.info(`Created cap-${rule.number}`);
        this.rulesViewer?.reloadViewer();
      }
      catch (e: any) {
        grok.shell.error(`Capture rule: ${e?.message ?? e}`);
      }
    });
    dialog.show();
    ui.tooltip.bind(dialog.getButton('OK'), () => problem() ?? 'Create the rule');
    refresh();
  }
}
