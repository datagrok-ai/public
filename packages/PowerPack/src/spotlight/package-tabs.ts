import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';

/** Spotlight tabs that packages add: a function with `meta.role: spotlightTab` returning a widget is a tab, named
 * by its friendly name, after the built-in ones. `meta.spotlightDefault: true` asks for the tab to be the one the
 * Spotlight opens on; it is opened the first time it appears, and from then on the tab the user last picked wins.
 * A tab keeps the Spotlight's height; what it finds goes to the preview below the Spotlight, through
 * PowerPack:showSpotlightPreview, as the Workspace tab's previews do. */
export const SPOTLIGHT_TAB_ROLE = 'spotlightTab';

export function spotlightTabFunctions(): DG.Func[] {
  return DG.Func.find({meta: {role: SPOTLIGHT_TAB_ROLE}});
}

export async function spotlightTabContent(func: DG.Func): Promise<HTMLElement> {
  try {
    const widget: DG.Widget = await func.apply();
    return widget?.root ?? ui.divText(`${func.friendlyName} returned nothing to show`, 'pp-empty-message');
  }
  catch (e) {
    console.error(e);
    return ui.divText(`${func.friendlyName} could not be shown: ${e instanceof Error ? e.message : e}`, 'pp-empty-message');
  }
}

const SETTINGS = 'power-pack-spotlight-tabs';
const OPENED_KEY = 'openedDefaultTabs';

/** `funcs` are the functions that got a tab. */
export function openDefaultTabOnce(tabs: DG.TabControl, funcs: DG.Func[]): void {
  const func = funcs.find((f) => `${f.options['spotlightDefault']}` === 'true');
  const opened: string[] = JSON.parse(grok.userSettings.getValue(SETTINGS, OPENED_KEY) ?? '[]');
  if (!func || opened.includes(func.nqName))
    return;
  tabs.currentPane = tabs.getPane(func.friendlyName);
  grok.userSettings.add(SETTINGS, OPENED_KEY, JSON.stringify([...opened, func.nqName]));
}
