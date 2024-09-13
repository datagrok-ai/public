/* Do not change these import lines. Datagrok will import API library in exactly the same manner */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const WIDGETS_STORAGE = 'widgets';

export interface UserWidgetSettings {
  factoryName?: string;
  caption?: string;
  ignored?: boolean;
}

export interface UserWidgetsSettings {
  [index: string]: UserWidgetSettings;
}

export let settings: UserWidgetsSettings;

export function getSettings(): UserWidgetsSettings {
  if (!settings) {
    const savedSettings: {[key: string]: any} = grok.userSettings.get(WIDGETS_STORAGE) ?? {};
    for (const key of Object.keys(savedSettings))
      savedSettings[key] = JSON.parse(savedSettings[key]);
    settings = savedSettings;
  }
  return settings;
}

export function saveSettings(): void {
  console.log(settings);
  let s: {[key: string]: any} = {};
  for (const key of Object.keys(settings))
    s[key] = JSON.stringify(settings[key]);
  grok.userSettings.addAll(WIDGETS_STORAGE, s);
}

export function widgetHost(w: DG.Widget, widgetHeader?: HTMLDivElement): HTMLElement {
  const host = ui.box(null, 'power-pack-widget-host');
  widgetHeader ??= ui.div();
  function remove(): void {
    host.remove();
    if (w.factory?.name) {
      const widgetSettings = settings[w.factory.name] ?? (settings[w.factory.name] = { });
      widgetSettings.ignored = true;
      saveSettings();
      grok.shell.info('To control widget visibility, go to Tools | Widgets');
    }
  }

  if (w.props.hasProperty('order'))
    host.style.order = w.props.order;


  const header = ui.div([
    ui.divText(w.props.hasProperty('caption') ? w.props.caption : '', 'd4-dialog-title'),
    widgetHeader,
    ui.icons.settings(() => {grok.shell.o = w;}, 'Edit settings'),
    ui.icons.close(remove, 'Remove'),
  ], 'd4-dialog-header');

  if (w.root.classList.contains('widget-narrow'))
    host.classList.add('widget-narrow');
  if (w.root.classList.contains('widget-wide'))
    host.classList.add('widget-wide');

  host.appendChild(header);
  host.appendChild(ui.box(w.root, 'power-pack-widget-content'));
  ui.tools.setHoverVisibility(host, Array.from(host.querySelectorAll('i')));

  return host;
}
