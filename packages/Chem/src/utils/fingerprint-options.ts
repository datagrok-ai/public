/**
 * Fingerprint options: the optional parameters of `getFingerprints` / `chemGetFingerprints` and the RDKit details they
 * become (RDKit's MinimalLib names, `Code/MinimalLib/common.h`: `radius`, `nBits`, `useChirality`). The defaults:
 * Morgan radius 2 (ECFP4), 2048 bits, no chirality, no features. No `datagrok-api` here: the
 * RDKit and Crux workers import it.
 */
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {defaultMorganFpLength, defaultMorganFpRadius, Fingerprint} from './chem-common';

export interface FingerprintOptions {
  /** Morgan: the radius (2: ECFP4, 3: ECFP6). */
  radius?: number | null;
  /** The number of bits (`nBits`): every type but MACCS (167 keys). */
  size?: number | null;
  /** Morgan, atom pairs, topological torsions: stereocentres take part (as RDKit encodes them); Crux computes
   * them in both modes ({@link needsCrux}). */
  useChirality?: boolean | null;
  /** Morgan: feature invariants (FCFP) instead of atom invariants (ECFP); Crux computes them in both modes. */
  useFeatures?: boolean | null;
}

/** The options each type reads; the others are ignored for it. */
const READS: {[type in Fingerprint]: (keyof FingerprintOptions)[]} = {
  [Fingerprint.Morgan]: ['radius', 'size', 'useChirality', 'useFeatures'],
  [Fingerprint.RDKit]: ['size'],
  [Fingerprint.Pattern]: ['size'],
  [Fingerprint.MACCS]: [],
  [Fingerprint.AtomPair]: ['size', 'useChirality'],
  [Fingerprint.TopologicalTorsion]: ['size', 'useChirality'],
};

/** The whole-number options: their defaults and ranges. */
const NUMBERS = {
  radius: {default: defaultMorganFpRadius, min: 0, max: 10},
  size: {default: defaultMorganFpLength, min: 1, max: 65536},
} as const;

/**
 * The options that differ from the defaults and that `type` reads (null when none: the default fingerprint). Values may
 * come as text (a function call's parameters: '3', 'true'); a radius or size that is not a whole number in range
 * throws.
 */
export function fingerprintOptions(type: Fingerprint, options?: FingerprintOptions | null):
  FingerprintOptions | null {
  if (!options)
    return null;
  const set: FingerprintOptions = {};
  for (const key of READS[type] ?? []) {
    const raw: unknown = options[key];
    if (raw === null || raw === undefined || raw === '')
      continue;
    if (key === 'radius' || key === 'size') {
      const {default: byDefault, min, max} = NUMBERS[key];
      const value = Number(raw);
      if (!Number.isInteger(value) || value < min || value > max)
        throw new Error(`Fingerprint ${key} must be a whole number from ${min} to ${max}, not ${raw}`);
      if (value !== byDefault)
        set[key] = value;
    } else if (raw === true || raw === 'true')
      set[key] = true;
  }
  return Object.keys(set).length ? set : null;
}

/**
 * Whether only Crux computes these fingerprints, in both modes, with no RDKit behind it (as molzip):
 * - FCFP (`useFeatures`): RDKit.js reads no `useFeatures` (MinimalLib's Morgan details) and gives ECFP; Crux gives
 *   RDKit's feature-invariant Morgan (`GetMorganFeatureAtomInvGen`) bit for bit;
 * - `useChirality`: Chem's RDKit.js 2024.09 under the new stereo perception (Chem's) gives an enantiomer pair the same
 *   Morgan bits and ignores it for atom pairs and torsions (its generators read the legacy perception's CIP codes);
 *   Crux gives RDKit 2026.03's, which label the molecule by the new CIP labeller first.
 */
export function needsCrux(options: FingerprintOptions | null): boolean {
  return !!options?.useFeatures || !!options?.useChirality;
}

/**
 * The part of a cached fingerprint column's name for these options (`fingerprintOptions`' result): none for the
 * defaults, so saved projects' caches (`~smiles.Morgan`) stay valid; e.g. `(radius=3, size=1024, chirality, features)`.
 */
export function fingerprintOptionsKey(options: FingerprintOptions | null): string {
  if (!options)
    return '';
  const parts: string[] = [];
  if (options.radius !== undefined)
    parts.push(`radius=${options.radius}`);
  if (options.size !== undefined)
    parts.push(`size=${options.size}`);
  if (options.useChirality)
    parts.push('chirality');
  if (options.useFeatures)
    parts.push('features');
  return `(${parts.join(', ')})`;
}

/**
 * The details JSON of RDKit's `get_<type>_fp_as_uint8array(details)` for these options (`fingerprintOptions`' result),
 * or null for the call without details. Morgan always has its details, as before (`{"radius":2,"nBits":2048}`).
 */
export function fingerprintDetails(type: Fingerprint, options: FingerprintOptions | null): string | null {
  if (type === Fingerprint.Morgan) {
    return JSON.stringify({radius: options?.radius ?? defaultMorganFpRadius,
      nBits: options?.size ?? defaultMorganFpLength, ...(options?.useChirality ? {useChirality: true} : {}),
      ...(options?.useFeatures ? {useFeatures: true} : {})});
  }
  if (!options)
    return null;
  return JSON.stringify({...(options.size !== undefined ? {nBits: options.size} : {}),
    ...(options.useChirality ? {useChirality: true} : {})});
}

/** The module's `get_<type>_fp_as_uint8array` method for a fingerprint type. */
const FP_METHOD: {[type in Fingerprint]: string} = {
  [Fingerprint.Morgan]: 'get_morgan_fp_as_uint8array',
  [Fingerprint.RDKit]: 'get_rdkit_fp_as_uint8array',
  [Fingerprint.Pattern]: 'get_pattern_fp_as_uint8array',
  [Fingerprint.MACCS]: 'get_maccs_fp_as_uint8array',
  [Fingerprint.AtomPair]: 'get_atom_pair_fp_as_uint8array',
  [Fingerprint.TopologicalTorsion]: 'get_topological_torsion_fp_as_uint8array',
};

/**
 * Throws the module's refusal of these options (a fingerprint of ethanol with them): the worker code catches errors
 * per molecule, so a refused option would give every molecule no fingerprint instead of the caller's fallback.
 */
export function assertFingerprintOptions(module: RDModule, type: Fingerprint, options?: FingerprintOptions | null):
  void {
  const opts = fingerprintOptions(type, options);
  if (!opts) // the defaults, which every module computes
    return;
  const details = fingerprintDetails(type, opts)!;
  const mol = module.get_mol('CCO');
  try {
    (mol as any)[FP_METHOD[type]](details);
  } finally {
    mol.delete();
  }
}
