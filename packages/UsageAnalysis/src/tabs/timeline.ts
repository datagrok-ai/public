import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {UaView} from './ua';
import {UaToolbox} from '../ua-toolbox';
import {ViewHandler} from '../view-handler';
import {scrollToStartOnFirstDraw} from '../utils';
import '../../css/usage_analysis.css';

export const TIMELINE_KEYS = ['action', 'request', 'session', 'report', 'rule'];

/** One action, request, session, report or capture rule in time order, from the server's `Timeline` function. */
export class TimelineView extends UaView {
  keyInput = ui.input.choice('By', {value: 'action', items: TIMELINE_KEYS, nullable: false});
  idInput = ui.input.string('Id');
  host: HTMLDivElement = ui.box();
  /** `<key>=<id>` of the record shown, for the tab's path. */
  shownParam?: string;
  private fromUrl: {key: string, value: string} | null = null;
  /** The app's `?<key>=<id>` parameters: the platform passes them to the app function, not in the URL. */
  static urlParams: {[key: string]: string | undefined} = {};

  constructor(uaToolbox?: UaToolbox) {
    super(uaToolbox);
    this.name = 'Timeline';
    const key = TIMELINE_KEYS.find((k) => TimelineView.urlParams[k]);
    if (key)
      this.fromUrl = {key, value: TimelineView.urlParams[key]!};
  }

  static open(handler: ViewHandler, key: string, value: string): void {
    handler.changeTab('Timeline');
    (handler.getCurrentView() as TimelineView).show(key, value);
  }

  async initViewers(path?: string): Promise<void> {
    this.idInput.input.addEventListener('keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this.load();
    });
    this.idInput.setTooltip('An action or request id (x-request-id), a session id, a report number or cap-<n>');
    const form = ui.form([this.keyInput, this.idInput], {classes: 'ua-toolbar'});
    form.append(ui.buttonsInput([ui.button('Show', () => this.load())]));
    this.root.append(ui.divV([ui.div([form], 'ua-toolbar'), this.host], 'ui-box'));
    if (!this.shownParam && this.fromUrl)
      this.show(this.fromUrl.key, this.fromUrl.value);
  }

  show(key: string, value: string): void {
    this.keyInput.value = key;
    this.idInput.value = value;
    this.load();
  }

  load(): void {
    const key = this.keyInput.value!;
    const value = (this.idInput.value ?? '').trim();
    if (!value)
      return;
    this.shownParam = `${key}=${encodeURIComponent(value)}`;
    if (this.uaToolbox?.viewHandler.getCurrentView() === this)
      this.uaToolbox.viewHandler.updatePath();
    ui.empty(this.host);
    this.host.append(ui.waitBox(async () => {
      try {
        const t: DG.DataFrame = await grok.functions.call('Timeline', {spec: JSON.stringify({[key]: value})});
        if (t.rowCount === 0)
          return ui.divText(`No records for ${key} ${value}`);
        t.name = `Timeline of ${key} ${value}`;
        const grid = DG.Viewer.grid(t, {showRowHeader: false, allowRowSelection: false, allowBlockSelection: false});
        grid.col('time')!.format = 'yyyy-MM-dd HH:mm:ss.fff';
        grid.col('summary')!.width = 500;
        scrollToStartOnFirstDraw(grid);
        return grid.root;
      }
      catch (e: any) {
        return ui.divText(`Timeline: ${e?.message ?? e}`, 'd4-viewer-error');
      }
    }));
  }
}
