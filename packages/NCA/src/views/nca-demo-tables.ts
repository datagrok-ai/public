/**
 * The bundled demo-table catalogue — pure data, and deliberately a LEAF MODULE
 * (it imports nothing, not even the Datagrok API).
 *
 * **Why this is its own file.** Two reasons, one structural and one about
 * cohesion:
 *
 *  1. **It breaks an import cycle rather than working around one.** `views/`
 *     contains a pre-existing 4-module cycle — `nca-launch` →
 *     `nca-mapping-pane` (`nca-launch.ts:21`) → `nca-add-to-study`
 *     (`nca-mapping-pane.ts:47`) → `nca-analysis-open`
 *     (`nca-add-to-study.ts:35`) → back to `nca-launch`
 *     (`nca-analysis-open.ts:20`). That cycle has never crashed the bundle only
 *     because every export crossing it is a *function*, evaluated at call time.
 *     The ★ balloon copy needs these display names, and having
 *     `nca-hub-help.ts` import them from `nca-launch.ts` pulled the help module
 *     into the cycle: an eager `export const` there read `DEMO_TABLE_FILES`
 *     during module evaluation, observed an uninitialized binding, threw, and
 *     took the WHOLE package down — every function silently unregistered, with
 *     the only symptom being `Unable to find JS function "test"` on a green
 *     webpack build. Homing the data in a dependency-free module removes the
 *     edge outright, so consumers may hold it in plain `const`s again.
 *  2. **It is where this data belongs.** `nca-launch.ts` documents itself as
 *     "pure functions that compose the analysis flow"; a catalogue of bundled
 *     fixtures with scientist-facing display copy is not that.
 *
 * **The convention this file protects, worth knowing before adding to `views/`:**
 * inside a cyclic module group, a top-level `const` that reads another cyclic
 * module's binding is a bundle-killer, while a function that reads it at call
 * time is safe. Prefer importing shared data from a leaf like this one.
 */

/**
 * Frame-level tag carrying the demo entry's `file` — the STABLE identity key for
 * `ensureDemoTables`' idempotence. `df.name` is user-mutable (rename the table
 * and the next ★ click would add a duplicate); the file id is not.
 *
 * Dot-namespaced to avoid colliding with a platform-reserved frame tag, whose
 * failure mode (ours silently ignored, or a platform tag overwritten) would be
 * invisible. Note the package's COLUMN-level tags use dashes
 * (`nca-passthrough`, `nca-param`); these are the first FRAME-level tags here,
 * so there is no in-package precedent either way.
 */
export const DEMO_ID_TAG = 'nca.demo.id';

/**
 * Frame-level tag stating where a demo table's numbers came from. Set on BOTH
 * demo frames, not only the synthetic one: a tag present only on the simulated
 * table would make its ABSENCE the signal, which is a worse design than stating
 * provenance positively in both cases.
 *
 * **This is machine-readable bookkeeping, NOT user-facing disclosure.** Measured
 * on a live server 2026-07-28: a frame tag renders in no platform surface at all
 * (the properties panel was confirmed rendered and contained neither the value
 * nor the word "provenance"). It also does not travel with an export or a
 * screenshot. See `files/REGEN.md` before describing the demo as
 * "disclosed" anywhere.
 */
export const PROVENANCE_TAG = 'nca.demo.provenance';

/** One bundled demo table: where it loads from, what it is called, and where its
 *  numbers came from. */
export interface DemoTableEntry {
  /** Path under the package AppData root, e.g. `demo/indomethacin_iv.csv`. */
  readonly file: string;
  /** Display name stamped onto `df.name` — deliberately DECOUPLED from the
   *  filename so the table can read as a scientist would name it. */
  readonly name: string;
  /** Value for {@link PROVENANCE_TAG}: real lineage, or simulated + how. */
  readonly provenance: string;
  /** One clause saying what this table lets the user DO — the ★ balloon's
   *  reason to name it at all. Lives here, beside `name`, so the balloon is a
   *  pure projection of this constant and no copy is typed twice. */
  readonly unlocks: string;
  /**
   * Scientist-facing caveat, rendered italic in the ★ balloon. Present ONLY
   * where the data warrants one — absent means "nothing to disclose", never
   * "we forgot".
   *
   * **Hard budget: 4 words maximum** (user decision, 2026-07-28). The first
   * version of this caveat was a full sentence plus a meta-clause about its own
   * persistence; the user reviewed it on screen and asked for one short italic
   * line instead. Four words is the whole budget, so spend it on the two facts
   * that matter — the numbers and the structure — not on hedging. A test
   * enforces the count.
   */
  readonly qualifier?: string;
}

/**
 * The demo tables the Hub's ★ icon opens (topic: nca-hub-onboarding). Two
 * shapes, deliberately: an IV-bolus single-dose set and a multi-route set, so a
 * cold-start user can reach both the plain-NCA and the bioavailability paths
 * without hunting for data. Both carry SMILES, so the structure-passthrough
 * column path is exercisable straight from the demo.
 *
 * Order is load order AND seeding order — `ensureDemoTables` returns in this
 * order and the Hub seeds its picker with the first entry.
 *
 * **Display name is not derived from the filename.** The two are separate fields
 * on purpose: the file is the stable id, the name is scientist-facing copy
 * carrying compound + route + SPECIES (user decision, 2026-07-27 — species is
 * the highest-value discriminator for a tool whose two fixtures are one human
 * clinical dataset and one rat study). A parallel name array would let the two
 * drift by index; the balloon and the tests both read this constant.
 */
export const DEMO_TABLE_FILES: readonly DemoTableEntry[] = [
  {
    file: 'demo/indomethacin_iv.csv',
    name: 'Indomethacin IV (human)',
    provenance:
      'Real data. Indomethacin IV bolus, 6 subjects — Kwan et al. 1976, distributed as ' +
      'R datasets::Indometh; the bundled copy is bit-exact against the PKNCA 0.12.1 ' +
      'reference suite. The indomethacin structure is the actual compound studied.',
    unlocks: '6 subjects, IV bolus. Compute basic PK parameters (CL, Vz, t½).',
  },
  {
    file: 'demo/propranolol_multiroute.csv',
    name: 'Propranolol IV/PO/SC (rat)',
    provenance:
      'SIMULATED data — a seeded one-compartment simulation (SEED=20260526; ' +
      'CL=0.5 L/h/kg, Vd=2.5 L/kg, t½≈3.5 h, F_PO=0.35, F_SC=0.85), not a real ' +
      'study. The propranolol structure is ILLUSTRATIVE ONLY: these numbers are ' +
      'dose-linear and low-clearance, whereas real rat propranolol is saturable and ' +
      'flow-limited (CL 43-94 mL/min/kg, Iwamoto 1985). Do not read the compound ' +
      'identity as a claim about the kinetics.',
    unlocks: '64 animals, three routes. Compare routes and compute bioavailability.',
    // Comma, not a semicolon: these are two noun phrases, not two independent
    // clauses, and the comma matches the marks already used in `unlocks`.
    // Terminal period, because the two sentences it follows both carry one.
    qualifier: 'Simulated data, illustrative structure.',
  },
];
