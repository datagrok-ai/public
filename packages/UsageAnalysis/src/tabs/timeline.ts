import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import {emptyState, scrollToStartOnFirstDraw} from '../utils';
import {ErrorsView} from './errors';
import '../../css/usage_analysis.css';

export const TIMELINE_KEYS = ['action', 'request', 'session', 'report', 'rule'];

/** One action, request, session, report or capture rule in time order, from `grok.dapi.log.getTimeline`. */
export class TimelineView extends UaView {
  keyInput = ui.input.choice('By', {value: 'action', items: TIMELINE_KEYS, nullable: false});
  idInput = ui.input.string('Id');
  host: HTMLDivElement = ui.box();
  /** `<key>=<id>` of the record shown, for the tab's path. */
  shownParam?: string;
  /** The app's `?<key>=<id>` parameters: the platform passes them to the app function, not in the URL. */
  static urlParams: {[key: string]: string | undefined} = {};

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Timeline';
  }

  async initViewers(path?: string): Promise<void> {
    this.idInput.input.addEventListener('keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this.load();
    });
    this.idInput.setTooltip('An action id (Clicks), a request id (Errors), a session id, a report number or cap-<n>');
    const toolbar = ui.divH([this.keyInput.root, this.idInput.root, ui.button('Show', () => this.load())],
      'ua-toolbar ua-timeline-toolbar');
    this.root.append(ui.divV([toolbar, this.host], 'ui-box'));
    if (this.shownParam)
      return;
    const key = TIMELINE_KEYS.find((k) => TimelineView.urlParams[k]) ?? 'action';
    this.show(key, TimelineView.urlParams[key] ?? '');
  }

  static noId(): HTMLElement {
    return emptyState('Enter an id and press Show',
      'Or open Timeline from a click, an error occurrence, a session or a capture rule');
  }

  show(key: string, value: string): void {
    this.keyInput.value = key;
    this.idInput.value = value;
    this.load();
  }

  load(): void {
    const key = this.keyInput.value!;
    const value = (this.idInput.value ?? '').trim();
    ui.empty(this.host);
    if (!value) {
      this.host.append(TimelineView.noId());
      return;
    }
    this.shownParam = `${key}=${encodeURIComponent(value)}`;
    if (this.uaToolbox?.viewHandler.getCurrentView() === this)
      this.uaToolbox.viewHandler.updatePath();
    this.host.append(ui.waitBox(async () => {
      try {
        const t = ErrorsView.frame(await grok.dapi.log.getTimeline({[key]: value}));
        if (t.rowCount === 0)
          return emptyState(`No records for ${key} ${value}`, 'Check the id, or choose another By');
        t.name = `Timeline of ${key} ${value}`;
        const grid = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
        grid.columns.setOrder(['time', 'source', 'server', 'kind', 'summary', 'status', 'ms', 'requestId', 'user']);
        grid.col('time')!.format = 'yyyy-MM-dd HH:mm:ss.fff UTC';
        grid.col('time')!.width = 190;
        grid.col('summary')!.width = 500;
        grid.col('requestId')!.name = 'request id';
        scrollToStartOnFirstDraw(grid);
        return grid.root;
      }
      catch (e: any) {
        return ui.divText(`Timeline: ${e?.message ?? e}`, 'd4-viewer-error');
      }
    }));
  }
}
