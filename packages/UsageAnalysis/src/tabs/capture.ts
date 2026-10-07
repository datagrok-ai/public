import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import {UaFilterableQueryViewer} from '../viewers/ua-filterable-query-viewer';
import {formatTime, onRowContextMenu, problemLine, scrollToStartOnFirstDraw, showProblem} from '../utils';
import '../../css/usage_analysis.css';

const CAPTURE_ITEMS = ['clicks', 'inputs', 'requests', 'calls', 'errors'];
const SUBJECTS = ['user', 'group', 'package', 'everyone'];
const ALL_ACTIVITY = 'all activity';
const SCOPES = [ALL_ACTIVITY, 'view', 'element', 'function', 'error'];
const TRIGGERS: {[scope: string]: string} = {view: 'the view opens', element: 'the element is clicked',
  function: 'the function runs', error: 'the error occurs'};
const NO_LEVEL = 'none';
const LEVELS = [NO_LEVEL, 'error', 'warning', 'info', 'debug'];
const CREDENTIALS_FLAG = 'credentials';
const DURATIONS: {[name: string]: number} = {'30 min': 30, '2 h': 120, '1 d': 1440, '2 d': 2880, '7 d': 10080};
/** The server's field names in a refusal, as the New rule dialog labels them. */
const LABELS: [RegExp, string][] = [[/\bsubject\.type\b/g, 'Subject'], [/\bsubject\.value\b/g, 'Who'],
  [/\bscope\.type\b/g, 'Scope'], [/\bscope\.value\b/g, 'Scope value'], [/\b(capture\.)?serverLevel\b/g, 'Server level'],
  [/\b(capture\.)?debugFlags\b/g, 'Debug flags'], [/\bmaxEvents\b/g, 'Max events'],
  [/\b(forMinutes|expiresAt)\b/g, 'For']];

/** Capture rules (`capture_rules`): who is captured, by whom and why; new rules and stops go through the server. */
export class CaptureView extends UaView {
  rulesViewer?: UaFilterableQueryViewer;
  /** The rule to show again in the context panel once the rules reload. */
  reshow?: string;
  /** Reasons of the rules stopped here: the server saves its stop record a moment after the stop. */
  stopReasons: {[rule: string]: string} = {};
  /** The context panel of the rule shown last. */
  shown?: HTMLElement;
  /** The New rule dialog's debug flags, read once. */
  debugFlags?: Promise<string[]>;

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
        const rules: string[] = t.col('rule')!.toList();
        const stopReasons = t.col('stop_reason')!;
        for (let i = 0; i < t.rowCount; i++) {
          if (!stopReasons.get(i) && this.stopReasons[rules[i]] && t.get('status', i) === 'stopped')
            stopReasons.set(i, this.stopReasons[rules[i]], false);
        }
        t.onCurrentRowChanged.subscribe(() => this.showRule(t, t.currentRowIdx));
        const i = this.reshow ? rules.indexOf(this.reshow) : -1;
        this.reshow = undefined;
        if (i >= 0) {
          t.currentRowIdx = i;
          this.showRule(t, i);
        }
        else if (this.shown && grok.shell.o === this.shown)
          grok.shell.setCurrentObject(null, false, true);
        return t;
      },
      createViewer: (t: DG.DataFrame) => {
        const grid = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
        const order = ['rule', 'name', 'author', 'subject', 'scope', 'reason', 'active', 'stop_reason', 'events'];
        grid.columns.setVisible(order);
        grid.columns.setOrder(order);
        grid.col('name')!.width = 150;
        grid.col('reason')!.width = 250;
        grid.col('active')!.width = 170;
        grid.col('stop_reason')!.width = 200;
        grid.col('stop_reason')!.name = 'stop reason';
        grid.onCellPrepare((gc) => {
          if (!gc.isTableCell || gc.gridColumn.column?.name !== 'active' || gc.cell.value !== 'active')
            return;
          const ends = formatTime(t.get('expires_at', gc.cell.rowIndex));
          const today = ends.startsWith(formatTime(new Date()).substring(0, 10));
          gc.customText = `active · ends ${ends.substring(today ? 11 : 5, 16)} UTC`;
        });
        onRowContextMenu(grid, (menu, i) => {
          if (t.get('status', i) === 'active')
            menu.item('Stop...', () => this.stopDialog(t.get('rule', i)));
          menu.item('Timeline', () => this.openTimeline('rule', t.get('rule', i)));
        });
        scrollToStartOnFirstDraw(grid);
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
    const status: string = t.get('status', i);
    const active = status === 'active';
    const buttons = [ui.button('Timeline', () => this.openTimeline('rule', rule))];
    if (active)
      buttons.unshift(ui.button('Stop...', () => this.stopDialog(rule)));
    const details: {[key: string]: string} = {
      'Name': t.get('name', i) ?? '',
      'Author': t.get('author', i) ?? '',
      'Subject': t.get('subject', i),
      'Scope': t.get('scope', i),
      'Capture': t.get('capture', i),
      'Reason': t.get('reason', i),
      'Status': status,
      'Created': formatTime(t.get('created_at', i)),
      [active ? 'Expires' : 'Ended']:
        formatTime(t.get(active ? 'expires_at' : 'ended_at', i) ?? t.get('expires_at', i)),
    };
    if (t.get('stopped_by', i))
      details['Stopped by'] = t.get('stopped_by', i);
    if (t.get('stop_reason', i))
      details['Stop reason'] = t.get('stop_reason', i);
    details['Events'] = `${t.get('events', i)} of ${t.get('max_events', i)}`;
    const perSession: {[key: string]: HTMLElement} = {'Per session': CaptureView.perSession(t.get('scope', i),
      t.get('window_minutes', i), t.get('max_sessions', i))};
    const acc = DG.Accordion.create();
    acc.addPane(rule, () => ui.divV([ui.tableFromMap({...details, ...perSession}), ui.buttonsInput(buttons)]), true);
    this.shown = acc.root;
    // `grok.shell.o = ...` ignores a change within a second of the last one, so a quick second click showed the first row.
    grok.shell.setCurrentObject(acc.root, true, true);
  }

  /** What a rule of [scope] captures of each session, with the server's semantics of the window and the limit. */
  static perSession(scope: string, windowMinutes: number, maxSessions: number): HTMLElement {
    const type = (scope ?? '').split(' ')[0];
    const trigger = TRIGGERS[type];
    const text = trigger ? `${windowMinutes} min from when ${trigger}, up to ${maxSessions} sessions` :
      'every session, until the rule ends';
    return ui.tooltip.bind(ui.divText(text), trigger ?
      `When ${trigger} (${scope}) in a session, that session is captured for the next ${windowMinutes} min, and ` +
      `each repeat restarts the ${windowMinutes} min; never past the rule's end. Only the first ${maxSessions} ` +
      'sessions it happens in are captured.' :
      `A rule without a scope captures every session of its subject until the rule ends; the ${windowMinutes}-min ` +
      `window and the ${maxSessions}-session limit apply to scoped rules only.`);
  }

  stopDialog(rule: string): void {
    const reason = ui.input.string('Reason', {tooltipText: 'Why the rule stops; required'});
    const line = problemLine();
    const dialog = ui.dialog(`Stop ${rule}`)
      .add(reason)
      .add(line)
      .onOK(async () => {
        try {
          await grok.dapi.log.stopCaptureRule(rule, reason.value.trim());
          grok.shell.info(`Stopped ${rule}`);
          this.stopReasons[rule] = reason.value.trim();
          this.reshow = rule;
          this.rulesViewer?.reloadViewer();
        }
        catch (e: any) {
          grok.shell.error(`${rule}: ${e?.message ?? e}`);
        }
      })
      .show();
    const ok = dialog.getButton('OK');
    const refresh = () => showProblem(ok, line, reason.value?.trim() ? null : 'Enter the reason');
    ui.tooltip.bind(ok, 'Stop the rule');
    reason.onChanged.subscribe(() => refresh());
    refresh();
    reason.input.focus();
  }

  /** [message] of a refused rule with the server's field names replaced by the dialog's labels. */
  static refusal(message: string): string {
    for (const [name, label] of LABELS)
      message = message.replace(name, label);
    return message;
  }

  /** The server's debug flags (`grok.dapi.log.getLoggingPolicy`), in their order, but `credentials`. */
  static async loadDebugFlags(): Promise<string[]> {
    const policy = await grok.dapi.log.getLoggingPolicy();
    return (policy.debugFlags ?? []).filter((f: string) => f !== CREDENTIALS_FLAG);
  }

  async newRuleDialog(): Promise<void> {
    this.debugFlags ??= CaptureView.loadDebugFlags().catch(() => {
      this.debugFlags = undefined;
      return [];
    });
    const debugFlags = await this.debugFlags;
    const subject = ui.input.choice('Subject', {value: 'user', items: SUBJECTS, nullable: false});
    const names: {[subject: string]: Promise<string[]>} = {};
    const namesOf = (type: string): Promise<string[]> => names[type] ??= (type === 'user' ?
      grok.dapi.users.list().then((users) => users.map((u) => u.login)) : type === 'group' ?
        grok.dapi.groups.list().then((groups) => groups.filter((g) => !g.personal).map((g) => g.friendlyName)) :
        grok.dapi.packages.list().then((packages) => packages.map((p) => p.name)))
      .then((list) => [...new Set(list)].sort());
    const who = ui.typeAhead('Who', {minLength: 1, limit: 30, source: async (query: string) => {
      const q = query.toLowerCase();
      return (await namesOf(subject.value!)).filter((n) => n.toLowerCase().includes(q));
    }});
    who.setTooltip('The login, the group name or the package name');
    const scope = ui.input.choice('Scope', {value: ALL_ACTIVITY, items: SCOPES, nullable: false});
    const scopeValue = ui.input.string('Scope value',
      {tooltipText: 'A view name, a part of an element path, a function nqName or an error signature'});
    const capture = ui.input.multiChoice('Capture', {value: ['clicks', 'requests', 'errors'], items: CAPTURE_ITEMS});
    const level = ui.input.choice('Server level', {value: NO_LEVEL, items: LEVELS, nullable: false});
    const flags = ui.input.multiChoice('Debug flags', {value: [], items: debugFlags,
      tooltipText: 'The server debug categories this rule turns on, as in Settings > Logger'});
    ui.setDisplay(flags.root, debugFlags.length > 0);
    const duration = ui.input.choice('For', {value: '1 d', items: Object.keys(DURATIONS), nullable: false});
    const maxEvents = ui.input.int('Max events', {value: 10000});
    const anonymous = ui.input.bool('Anonymous', {tooltipText: 'No user, session or IP'});
    const anonymousNote = ui.divText('Only for group and everyone rules', 'ua-note');
    const name = ui.input.string('Name');
    const reason = ui.input.string('Reason', {tooltipText: 'Why this activity is captured; required'});
    const line = problemLine();

    const problem = (): string | null => {
      if (subject.value !== 'everyone' && !who.value?.trim())
        return `Enter the ${subject.value === 'user' ? 'login' : `${subject.value} name`} in Who`;
      if (scope.value !== ALL_ACTIVITY && !scopeValue.value?.trim())
        return `Enter the ${scope.value} in Scope value`;
      if (!reason.value?.trim())
        return 'Enter the reason';
      return null;
    };
    const refresh = () => {
      ui.setDisplay(who.root, subject.value !== 'everyone');
      ui.setDisplay(scopeValue.root, scope.value !== ALL_ACTIVITY);
      const canBeAnonymous = subject.value === 'group' || subject.value === 'everyone';
      anonymous.enabled = canBeAnonymous;
      ui.setDisplay(anonymousNote, !canBeAnonymous);
      if (!canBeAnonymous && anonymous.value)
        anonymous.value = false;
      showProblem(ok, line, problem());
    };

    const dialog = ui.dialog('New capture rule');
    const inputs: DG.InputBase[] = [subject, who, scope, scopeValue, capture, level, flags, duration, maxEvents,
      anonymous, name, reason];
    for (const input of inputs) {
      dialog.add(input);
      if (input === anonymous)
        dialog.add(anonymousNote);
      input.onChanged.subscribe(() => refresh());
    }
    dialog.add(line);
    dialog.addButton('OK', async () => {
      const body = {
        name: name.value?.trim() || undefined,
        subject: {type: subject.value, value: subject.value === 'everyone' ? undefined : who.value.trim()},
        scope: scope.value === ALL_ACTIVITY ? undefined : {type: scope.value, value: scopeValue.value.trim()},
        capture: {
          ...Object.fromEntries(CAPTURE_ITEMS.map((item) => [item, capture.value?.includes(item) ?? false])),
          serverLevel: level.value === NO_LEVEL ? null : level.value,
          debugFlags: flags.value ?? [],
        },
        anonymous: anonymous.value,
        maxEvents: maxEvents.value,
        forMinutes: DURATIONS[duration.value!],
        reason: reason.value.trim(),
      };
      ok.disabled = true;
      try {
        const rule = await grok.dapi.log.addCaptureRule(body);
        dialog.close();
        grok.shell.info(`Created cap-${rule.number}`);
        this.reshow = `cap-${rule.number}`;
        this.rulesViewer?.reloadViewer();
      }
      catch (e: any) {
        ok.disabled = false;
        line.textContent = CaptureView.refusal(`${e?.message ?? e}`);
        ui.setDisplay(line, true);
      }
    });
    const ok = dialog.getButton('OK');
    ui.tooltip.bind(ok, 'Create the rule');
    dialog.show();
    refresh();
    who.input.focus();
  }
}
