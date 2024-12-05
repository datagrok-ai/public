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


function initWidgetHost(host: HTMLDivElement, w: DG.Widget) {
  function remove(): void {
    host.remove();
    if (w.factory?.name) {
      const widgetSettings = settings[w.factory.name] ?? (settings[w.factory.name] = { });
      widgetSettings.ignored = true;
      saveSettings();
    }
  }

  if (w.props.hasProperty('order'))
    host.style.order = w.props.order;

  const header = host.querySelector('.d4-dialog-header')!;
  header.appendChild(ui.icons.settings(() => {grok.shell.o = w;}, 'Edit settings'));
  header.appendChild(ui.icons.close(remove, 'Remove'));

  if (w.root.classList.contains('widget-narrow'))
    host.classList.add('widget-narrow');
  if (w.root.classList.contains('widget-wide'))
    host.classList.add('widget-wide');

  host.querySelector('.power-pack-widget-content')!.appendChild(w.root);
  ui.tools.setHoverVisibility(host, Array.from(host.querySelectorAll('i')));
}

function createWidgetHost(title: string): HTMLDivElement {
  const header = ui.div([ui.divText(title, 'd4-dialog-title'),], 'd4-dialog-header');
  const host = ui.box(null, 'power-pack-widget-host');
  host.appendChild(header);
  host.appendChild(ui.box(null, 'power-pack-widget-content'));
  return host;
}

export function widgetHostFromFunc(f: DG.Func) {
  const host: HTMLDivElement = createWidgetHost(f.friendlyName);
  const contentDiv: HTMLElement = (host.querySelector('.power-pack-widget-content')!) as HTMLElement;

  f.apply().then(function(w: DG.Widget) {
      if (w) {
        w.factory = f;
        initWidgetHost(host, w);
      }
      else
        host.remove();
  })
  .catch((e) => {
    host.style.display = 'none';
    host.remove();
    console.error(`Error creating widget ${f.name}`, e);
  })
 .finally(() => ui.setUpdateIndicator(contentDiv, false, ''));

  setTimeout(() => {
    if (contentDiv!.children.length == 0)
      ui.setUpdateIndicator(contentDiv, true, '');
  }, 1000);

  if (f.options['order'] !== null)
    host.style.order = f.options['order'];

  return host;
}


export function widgetHost(w: DG.Widget/*, widgetHeader?: HTMLDivElement*/): HTMLElement {
  const host = createWidgetHost(w.props.caption ?? '');
  initWidgetHost(host, w);
  //widgetHeader ??= ui.div();

  return host;
}
