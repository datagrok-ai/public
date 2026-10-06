import {DEFAULT_MONOMERS_PER_ROW} from '../constants';

/**
 * Read the `MonomersPerRow` package property into a usable row width.
 *
 * Datagrok hands an `int` property back as a number once it has round-tripped
 * through the server, but a `defaultValue` straight out of package.json arrives
 * as the string written there, and a package whose settings have not loaded yet
 * gives `undefined`. All three have to end up at a sane integer — a `NaN` here
 * would reach the layout and produce a drawing with no coordinates.
 *
 * @param {unknown} raw The raw property value.
 * @return {number} A row width of at least 1, or {@link DEFAULT_MONOMERS_PER_ROW}.
 */
export function resolveMonomersPerRow(raw: unknown): number {
  const parsed = typeof raw === 'number' ? raw : Number.parseInt(String(raw ?? ''), 10);
  return Number.isFinite(parsed) && parsed >= 1 ? Math.round(parsed) : DEFAULT_MONOMERS_PER_ROW;
}
