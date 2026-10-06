/** Finding a program across the portfolios the user can see, and the programs they opened last,
 * kept in their user settings (they follow the user to another browser). The rules are in
 * ./model/find. */
import {portalDb} from './generated/db';
import {addRecent, FindProgram, parseRecents, RecentProgram} from './model/find';
import {readPreference, writePreference} from './personas';
import {ROW_LIMIT} from './rows';

const RECENT_KEY = 'recent';

/** Every program the user can see, with its portfolio, in the shape the search reads. */
export async function findablePrograms(): Promise<FindProgram[]> {
  const [portfolios, programs] = await Promise.all([portalDb.portfolios.query({limit: ROW_LIMIT}),
    portalDb.programs.query({limit: ROW_LIMIT})]);
  const byId = new Map(portfolios.map((p) => [p.id, p]));
  return programs.filter((p) => byId.has(p.portfolio_id)).map((p) => ({id: p.id, code: p.code, name: p.name,
    portfolio: byId.get(p.portfolio_id)!.code, portfolioName: byId.get(p.portfolio_id)!.name, status: p.status ?? null,
    target: p.target ?? null, details: (p.details ?? null) as Record<string, unknown> | null}));
}

/** The programs the user opened last, newest first, as saved. */
export const recentPrograms = (): RecentProgram[] => parseRecents(readPreference(RECENT_KEY));

/** Puts the program first in the user's recents. */
export function rememberProgram(portfolio: string, id: string): void {
  writePreference(RECENT_KEY, JSON.stringify(addRecent(recentPrograms(), {id, portfolio, at: new Date().toISOString()})));
}

/** Forgets the user's recents. */
export const forgetRecents = (): void => writePreference(RECENT_KEY, null);
