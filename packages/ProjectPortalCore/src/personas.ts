/** Mode and persona against the platform: remembered per user in the platform's user settings
 * (`grok.userSettings`, saved on the server, so they follow the user to another browser), and the
 * user's groups for a portfolio's persona mapping. The rules are in ./model/personas. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {PortfolioRow} from './generated/db';
import {funcOf} from './catalogue';
import {asMode, asPersona, defaultPersona, Persona, PersonaConfig, PortalMode} from './model/personas';

/** The user-settings entry the portal's preferences live under. */
export const PREFERENCES = 'ProjectPortalCore';
/** Advanced until the user chooses: the program page shows everything, as it did before modes. */
export const DEFAULT_MODE: PortalMode = 'advanced';

export const readPreference = (key: string): string | null => {
  try {
    return grok.userSettings.getValue(PREFERENCES, key) ?? null;
  } catch {
    return null;
  }
};

/** Writes a preference and saves it to the server in the background; a failed save is retried
 * by the platform, so the page never waits for it. */
export const writePreference = (key: string, value: string | null): void => {
  if (value === null) {
    // deleting a key that is not there throws
    if (readPreference(key) !== null)
      grok.userSettings.delete(PREFERENCES, key);
  } else
    grok.userSettings.add(PREFERENCES, key, value);
  void grok.userSettings.flush().catch((e) => grok.log.warning(`Project Portal: preference ${key} not saved: ${e?.message ?? e}`));
};

export const savedMode = (): PortalMode => asMode(readPreference('mode')) ?? DEFAULT_MODE;
export const saveMode = (mode: PortalMode): void => writePreference('mode', mode);

/** The persona the user chose, or null when they never did. */
export const savedPersona = (): Persona | null => asPersona(readPreference('persona'));
/** Remembers the user's persona; null forgets it, so the portfolio's mapping decides again. */
export const savePersona = (persona: Persona | null): void => writePreference('persona', persona);

let groupNames: Promise<string[]> | null = null;

/** The names of the groups the current user is in, transitively; read once per session. */
export function myGroupNames(): Promise<string[]> {
  groupNames ??= grok.dapi.groups.currentUserGroups()
    .then((gs) => [...new Set(gs.flatMap((g) => [g.friendlyName, g.name]).filter((n): n is string => !!n))])
    .catch(() => []);
  return groupNames;
}

/** The persona the program page uses in [portfolio]: the user's choice, else the portfolio's
 * mapping of the user's groups, else the portfolio's default, else by setting. */
export async function personaIn(portfolio: PortfolioRow): Promise<{persona: Persona, chosen: boolean}> {
  const chosen = savedPersona();
  if (chosen)
    return {persona: chosen, chosen: true};
  return {persona: defaultPersona(portfolio.personas as PersonaConfig | undefined, await myGroupNames(), portfolio.setting),
    chosen: false};
}

/** Runs a portfolio-declared card's function for a program: it gets `portfolio`, `program` and
 * `stage` when it has inputs by those names. A view it returns opens; a table opens as a table
 * view; anything else is shown as a message. */
export async function runFunctionCard(name: string, ctx: {portfolio: string, program: string, stage?: string | null}):
  Promise<void> {
  const func = funcOf(name);
  if (!func)
    throw new Error(`${name} is not installed here`);
  const params: Record<string, unknown> = {};
  for (const p of func.inputs) {
    const v = (ctx as Record<string, unknown>)[p.name];
    if (v !== undefined && v !== null)
      params[p.name] = v;
  }
  const result = await func.apply(params);
  if (result instanceof DG.ViewBase)
    grok.shell.addView(result);
  else if (result instanceof DG.DataFrame)
    grok.shell.addTableView(result);
  else if (result !== undefined && result !== null)
    grok.shell.info(String(result));
}
