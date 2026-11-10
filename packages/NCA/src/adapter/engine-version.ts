/**
 * The installed NCA engine version — `@datagrok-libraries/sci-comp`'s own
 * `package.json` version, read at BUILD time (webpack inlines the JSON; the
 * upstream package exports the `./package.json` subpath for exactly this).
 *
 * Why it exists (GROK-20960, F-N.6): a study's cached per-analysis results
 * (`summaryCache`) are invalidated only on SOURCE drift. Nothing recorded which
 * engine computed them, so a numbers-changing sci-comp release (0.11.0: BLQ
 * substitutions reaching the integrator; the route-aware IV-bolus dose-time
 * gate) would leave every open study silently showing pre-fix numbers until
 * each analysis happened to be reopened — thesis #4's forbidden silent state.
 * Every compute now stamps this string on the analysis entry
 * (`AnalysisEntry.engineVersion`), and `isStale` (`nca-repository.ts`)
 * treats a stamp that differs from it (or is missing) as stale — the
 * existing "stale — refreshes on open" cue, no new UI.
 *
 * A dependency-free leaf module (guardrails: module-graph footguns) so any
 * layer may import it at module-init time.
 */
import sciCompPkg from '@datagrok-libraries/sci-comp/package.json';

export const INSTALLED_ENGINE_VERSION: string = sciCompPkg.version;
