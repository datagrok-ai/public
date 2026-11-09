/**
 * Minimum in-memory audit log for UC-01 — scoped to PRD IR-12 content only.
 *
 * **Regulatory posture:** ship v1.0 as RUO. The full ALCOA+ /
 * 21 CFR Part 11 §11.10(e) 12-field schema, before/after vectors, cryptographic
 * record hashing, persistent backend, and electronic signatures are
 * **Phase 3 P2** per PRD §10. This module deliberately implements only the
 * PRD IR-12 scope — "log of all user-modified parameters (λz point selection,
 * BLQ strategy, excluded subjects)" — as a flat in-memory record with five
 * named fields plus an opaque `payload` slot.
 *
 * Migration target: when a GxP track starts, `AuditEntry` is replaced with the full
 * ALCOA+ / Part 11 12-field schema and the in-memory array with a persistent backend.
 * Keeping the IR-12 scope now avoids data-model rework later.
 */

/**
 * Compute-result status shared with sci-comp's `nca.ComputeResult.status`.
 * Mirrored here (not imported) so this module stays usable for non-compute
 * audit entries (e.g. a UI action like 'subject-excluded' that has no
 * associated ComputeResult yet).
 */
export type AuditStatus = 'ok' | 'partial' | 'failed';

/**
 * One audit-log entry — the five PRD IR-12 fields plus an opaque payload.
 *
 * Fields are exactly those listed in the UC-01 spec §3.1 row 15 ("timestamp,
 * subjectKey, action, ruleHash, status") and the cycle-state Phase 2 plan.
 * The 12-field ALCOA+ extension (`user_id`, `signature_hash`, before/after
 * vectors, record hash chain, ...) is out of scope per AD-04.
 *
 * - `timestamp` — milliseconds since epoch, captured at write time.
 * - `subjectKey` — composite profile identifier. For UC-01 single-compound
 *   slice this is the subjectId string; later phases extend to
 *   `subject:analyte:occasion` per OQ-2.
 * - `action` — open-string action discriminator. The three IR-12-scoped
 *   action categories are `'lambdaZ-points-changed'`, `'blq-strategy-changed'`,
 *   `'subject-excluded'`. Compute itself logs `'compute'`. New action strings
 *   may be added without a breaking change.
 * - `ruleHash` — stable string fingerprint of the `nca.NcaRules` that were
 *   in effect for this entry. The hashing algorithm is owned by the (Phase 4)
 *   controller; this module only carries the value. An empty string is valid
 *   for actions that have no associated rules.
 * - `status` — matches `ComputeResult.status` when `action === 'compute'`;
 *   `'ok'` is the conventional default for UI actions that succeed without a
 *   compute (e.g. a successful BLQ-strategy change).
 * - `payload` — optional opaque structured detail. Examples per action:
 *   `'lambdaZ-points-changed'` → `{before: Int32Array, after: Int32Array}`;
 *   `'subject-excluded'` → `{reason: string}`. Consumers are responsible for
 *   their own payload shapes; nothing in this module inspects `payload`.
 */
export interface AuditEntry {
  readonly timestamp: number;
  readonly subjectKey: string;
  readonly action: string;
  readonly ruleHash: string;
  readonly status: AuditStatus;
  readonly payload?: unknown;
  /**
   * NCA Studies tagging (AD-7). Optional + additive: coordinator-driven computes
   * set them so audit entries are attributable to a project/analysis; existing
   * single-analysis computes leave them undefined (byte-compatible with every
   * prior `append()` call). First-class fields rather than buried in `payload`
   * so a future audit view can query by project/analysis without parsing the
   * opaque slot.
   */
  readonly studyId?: string;
  readonly analysisId?: string;
}

/**
 * In-memory append-only audit log.
 *
 * One instance per package session; the (Phase 4) controller owns the
 * lifetime. No persistence, no eviction, no cap — UC-01 expected workload is
 * a few dozen subjects per compute and a handful of computes per session.
 *
 * Thread safety: single-threaded by construction (the Datagrok platform
 * dispatches user actions on the main thread; the worker pool — Phase 2+ in
 * the architecture, not UC-01 — never touches the audit log directly).
 */
export class AuditLog {
  private readonly _entries: AuditEntry[] = [];

  /** Append one entry. The entry is stored by reference; callers should not mutate it. */
  append(entry: AuditEntry): void {
    this._entries.push(entry);
  }

  /** Read-only snapshot of all entries in insertion order. */
  entries(): readonly AuditEntry[] {
    return this._entries;
  }
}

/**
 * Package-level singleton — "one instance per package session" per the module
 * head, and "the package's in-memory AuditLog" per UC-01 spec §9 row 14.
 * Production code paths (Phase 6 grok-function entry point, Phase 5 dialog
 * "Compute" button) read this via {@link getAuditLog} and never construct
 * their own. Tests replace it via {@link setAuditLog} for isolation.
 *
 * Owning the lifetime here (not at call sites) guarantees DoD §9 row 14 —
 * every `computeNca` invocation appends to *the package's* log — at the
 * module boundary, independent of whether the caller threaded a log through.
 */
let packageAuditLog: AuditLog = new AuditLog();

/** Read the package-level audit log. */
export function getAuditLog(): AuditLog {
  return packageAuditLog;
}

/**
 * Replace the package-level audit log. **Tests-only API.** Production code
 * must not call this; the singleton is the package's session-long instance.
 */
export function setAuditLog(log: AuditLog): void {
  packageAuditLog = log;
}
