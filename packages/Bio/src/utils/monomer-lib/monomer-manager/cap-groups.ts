/* The cap groups that close a monomer's R groups, and which one a library's definition means.
   Definitions arrive with fields missing, the SMILES written half a dozen ways ([*:1][H], [H][*:1],
   [O:2], O[R2], C[*:1]=C) or a word where the SMILES belongs ([*:3]Azide) — so the SMILES is read
   first when it is chemistry, and the name and the alternate ID only when it is not. */
import {RGroup} from '@datagrok-libraries/bio/src/types/monomer-library';
import {HELM_RGROUP_FIELDS} from '@datagrok-libraries/bio/src/utils/const';

export interface CapGroup {
  /** The name a library writes: `capGroupName` and the tail of `alternateId`. */
  name: string;
  /** What the cap is, in a word or two. */
  description: string;
  /** The cap with `[*]` where it attaches. */
  smiles: string;
  /** Other names for it, compared lower-case without punctuation. */
  aliases?: string[];
  /** Other ways a library writes the cap without its attachment point. */
  fragments?: string[];
}

/** The fifteen caps offered first, most used first: the order is the count over the HELM core
 * library and two production libraries (H, OH, then allyl, bromo, chloro, alkyne, azide, iodo),
 * then the other small caps a chemist reaches for. Allyl is `C=C–`, as those libraries write it. */
export const COMMON_CAP_GROUPS: CapGroup[] = [
  {name: 'H', description: 'hydrogen', smiles: '[*][H]', aliases: ['hydrogen', 'proton'], fragments: ['H', '[H]']},
  {name: 'OH', description: 'hydroxyl', smiles: 'O[*]',
    aliases: ['hydroxyl', 'hydroxy'], fragments: ['O', 'OH', '[OH]']},
  {name: 'Allyl', description: 'allyl', smiles: 'C=C[*]',
    aliases: ['alkene', 'alkenyl'], fragments: ['C=C', '[CH2]=[CH]']},
  {name: 'Br', description: 'bromo', smiles: 'Br[*]',
    aliases: ['bromo', 'bromine', 'bromide'], fragments: ['Br', '[Br]']},
  {name: 'Cl', description: 'chloro', smiles: 'Cl[*]',
    aliases: ['chloro', 'chlorine', 'chloride'], fragments: ['Cl', '[Cl]']},
  {name: 'Alkyne', description: 'ethynyl', smiles: 'C#C[*]',
    aliases: ['alkynyl', 'ethynyl', 'acetylene'], fragments: ['C#C', '[CH]#C']},
  {name: 'Azide', description: 'azide', smiles: '[N-]=[N+]=N[*]', aliases: ['azido', 'n3'],
    fragments: ['N3', '[N-]=[N+]=N', 'N=[N+]=[N-]', 'N=N=N']},
  {name: 'I', description: 'iodo', smiles: 'I[*]', aliases: ['iodo', 'iodine', 'iodide'], fragments: ['I', '[I]']},
  {name: 'F', description: 'fluoro', smiles: 'F[*]',
    aliases: ['fluoro', 'fluorine', 'fluoride'], fragments: ['F', '[F]']},
  {name: 'NH2', description: 'amino', smiles: 'N[*]', aliases: ['amino', 'amine'], fragments: ['N', 'NH2', '[NH2]']},
  {name: 'Me', description: 'methyl', smiles: 'C[*]', aliases: ['methyl', 'ch3'], fragments: ['C', 'CH3', '[CH3]']},
  {name: 'Ac', description: 'acetyl', smiles: 'CC(=O)[*]',
    aliases: ['acetyl'], fragments: ['CC(=O)', 'C(C)=O', 'C(=O)C', 'CC=O']},
  {name: 'OMe', description: 'methoxy', smiles: 'CO[*]', aliases: ['methoxy', 'och3'], fragments: ['CO', 'OC']},
  {name: 'SH', description: 'thiol', smiles: 'S[*]',
    aliases: ['thiol', 'sulfhydryl', 'mercapto'], fragments: ['S', 'SH', '[SH]']},
  {name: 'Et', description: 'ethyl', smiles: 'CC[*]', aliases: ['ethyl', 'c2h5'], fragments: ['CC']},
];

export interface CapInference {
  cap: CapGroup;
  /** One of {@link COMMON_CAP_GROUPS}; false for a cap the definition names that is not among them. */
  known: boolean;
  /** The field the cap was read from. */
  from: 'smiles' | 'name' | 'alternateId';
  /** A field that says otherwise, worth showing next to the choice. */
  note?: string;
}

/** Canonical SMILES of a valid SMILES, null for anything else; the caller's chemistry engine. */
export type Canonicalize = (smiles: string) => string | null;

/** Every way an attachment point is written: [*:1], [*1], [1*], [1*:1], [R1], [R], [*], *. */
const ATTACHMENT = /\[\*:?\d*\]|\[\d+\*(?::\d+)?\]|\[R\d*\]|\*/g;
/** A cap atom carrying the R number as its atom map: [O:2], [OH:2], [H:1]. */
const MAPPED_ATOM = /\[([A-Z][a-z]?)(H\d?)?:\d+\]/g;

/** A field of a definition read from a library file, whatever case its key was written in. */
const fieldOf = (rg: Partial<RGroup>, field: string): string => {
  const record: {[key: string]: string | undefined} = rg;
  const key = Object.keys(record).find((k) => k.toLowerCase() === field.toLowerCase());
  const value = key === undefined ? undefined : record[key];
  return typeof value === 'string' ? value.trim() : '';
};

const plain = (s: string): string => s.toLowerCase().replace(/[^a-z0-9]/g, '');

export class CapGroups {
  private static _canonical: Map<string, CapGroup> | undefined;

  /** The R number of a label (`R2` → 2), or null. */
  static numberOf(label: string): number | null {
    const m = /^R(\d+)$/i.exec(label.trim());
    return m ? Number.parseInt(m[1]) : null;
  }

  /** X: what some libraries write where the cap is not decided — a placeholder, not a cap group. */
  static isPlaceholderName(name: string): boolean {
    return name.trim().toUpperCase() === 'X';
  }

  static isPlaceholder(cap: CapGroup): boolean {
    return CapGroups.isPlaceholderName(cap.name) || CapGroups.isPlaceholderName(CapGroups.fragmentOf(cap.smiles));
  }

  /** The monomer with every R group closed by its cap: each attachment of the monomer and of its cap
   * becomes one ring-closure bond (`N[*:1]` + `[*:1][H]` → `N%91.[H]%91`), so a cap of any size joins
   * at the right atom. An R group with no usable cap is closed with H. Canonical when `canonicalize`
   * is given. */
  static capped(smiles: string, rgroups: Partial<RGroup>[], canonicalize?: Canonicalize): string {
    let monomer = smiles;
    const caps: string[] = [];
    const numbers = new Set(Array.from(smiles.matchAll(/\[\*:(\d+)\]/g), (m) => Number.parseInt(m[1])));
    for (const n of numbers) {
      const def = rgroups.find((rg) => CapGroups.labelOf(rg) === `R${n}`);
      const inferred = def === undefined ? null : CapGroups.infer(def, canonicalize);
      const usable = inferred !== null && (inferred.known || (canonicalize?.(inferred.cap.smiles) ?? null) !== null);
      const cap = CapGroups.numbered(usable ? inferred!.cap.smiles : COMMON_CAP_GROUPS[0].smiles, n);
      const closure = `%${n < 10 ? 90 + n : n}`;
      const joined = [CapGroups._close(monomer, `[*:${n}]`, closure), CapGroups._close(cap, `[*:${n}]`, closure)];
      if (joined[0] === null || joined[1] === null)
        continue;
      monomer = joined[0];
      caps.push(joined[1]);
    }
    const result = [monomer, ...caps].join('.');
    return canonicalize?.(result) ?? result;
  }

  /** `token` (an attachment) replaced by a ring-closure bond on the atom it hangs from; null where
   * that atom cannot be told. */
  private static _close(smiles: string, token: string, closure: string): string | null {
    if (smiles.startsWith(token)) {
      const rest = smiles.slice(token.length);
      const atom = /^(\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnops])/.exec(rest);
      return atom === null ? null : atom[0] + closure + rest.slice(atom[0].length);
    }
    if (smiles.includes(`(${token})`))
      return smiles.replace(`(${token})`, closure);
    return smiles.includes(token) ? smiles.replace(token, closure) : null;
  }

  static byName(name: string): CapGroup | null {
    const p = plain(name);
    if (p === '')
      return null;
    return COMMON_CAP_GROUPS.find((c) => plain(c.name) === p || (c.aliases ?? []).includes(p)) ?? null;
  }

  /** The cap SMILES without its attachment point: `O[*:2]` → `O`, `[O:2]` → `O`, `[*:3]Azide` → `Azide`. */
  static fragmentOf(smiles: string): string {
    return smiles.replace(ATTACHMENT, '').replace(MAPPED_ATOM, (_, el, h) => `${el}${h ?? ''}`)
      .replace(/\(\)/g, '').trim();
  }

  /** The cap SMILES with its attachment numbered `n`, whatever form it was written in. */
  static numbered(smiles: string, n: number): string {
    return smiles.replace(ATTACHMENT, `[*:${n}]`);
  }

  /** The R group a known cap makes on R`n`, in the form the HELM libraries write. */
  static rgroup(cap: CapGroup, n: number): RGroup {
    return {
      [HELM_RGROUP_FIELDS.CAP_GROUP_SMILES_UPPERCASE]: CapGroups.numbered(cap.smiles, n),
      [HELM_RGROUP_FIELDS.ALTERNATE_ID]: `R${n}-${cap.name}`,
      [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: cap.name,
      [HELM_RGROUP_FIELDS.LABEL]: `R${n}`,
    };
  }

  /** An R group nothing defines yet: the label alone. */
  static undefinedRGroup(label: string): RGroup {
    return {[HELM_RGROUP_FIELDS.ALTERNATE_ID]: '', [HELM_RGROUP_FIELDS.CAP_GROUP_NAME]: '',
      [HELM_RGROUP_FIELDS.LABEL]: label};
  }

  static smilesOf(rg: Partial<RGroup>): string {
    return fieldOf(rg, HELM_RGROUP_FIELDS.CAP_GROUP_SMILES);
  }

  static nameOf(rg: Partial<RGroup>): string {
    return fieldOf(rg, HELM_RGROUP_FIELDS.CAP_GROUP_NAME);
  }

  static alternateIdOf(rg: Partial<RGroup>): string {
    return fieldOf(rg, HELM_RGROUP_FIELDS.ALTERNATE_ID);
  }

  static labelOf(rg: Partial<RGroup>): string {
    return fieldOf(rg, HELM_RGROUP_FIELDS.LABEL);
  }

  /** Which cap `rg` means, or null when nothing in it names one. The SMILES decides while it is
   * chemistry; a known cap's name over different chemistry stays that chemistry, as a cap of its
   * own, so picking it never rewrites the structure. */
  static infer(rg: Partial<RGroup>, canonicalize?: Canonicalize): CapInference | null {
    const smiles = CapGroups.smilesOf(rg);
    const name = CapGroups.nameOf(rg);
    const alternateId = CapGroups.alternateIdOf(rg);
    const altName = alternateId.replace(/^R\d+-?/i, '');
    const named = CapGroups.byName(name) ?? CapGroups.byName(altName);
    const disagree = (cap: CapGroup): string | undefined =>
      named !== null && named !== cap ? `the name says ${named.name}` : undefined;

    const fragment = CapGroups.fragmentOf(smiles);
    if (fragment !== '') {
      const byFragment = COMMON_CAP_GROUPS.find((c) => (c.fragments ?? []).includes(fragment)) ??
        CapGroups.byName(fragment);
      if (byFragment !== null)
        return {cap: byFragment, known: true, from: 'smiles', note: disagree(byFragment)};
      const canonical = canonicalize?.(smiles.replace(ATTACHMENT, '[*]')) ?? null;
      if (canonical !== null) {
        const byChemistry = CapGroups._byCanonical(canonicalize!).get(canonical);
        if (byChemistry !== undefined)
          return {cap: byChemistry, known: true, from: 'smiles', note: disagree(byChemistry)};
        const own = name || altName || fragment;
        return {cap: {name: own, description: 'as defined', smiles: smiles.replace(ATTACHMENT, '[*]')},
          known: false, from: 'smiles', note: named !== null ? `not the usual ${named.name}` : undefined};
      }
      // a SMILES that is no chemistry and no cap's name is what will be saved, so it stays the cap —
      // falling back to a name beside it would drop it silently
      // X is the cap left undecided: kept as written, so the row can say so, not swapped for a name beside it
      if (CapGroups.isPlaceholderName(fragment)) {
        const cap = {name: 'X', description: 'placeholder', smiles: smiles.replace(ATTACHMENT, '[*]')};
        return {cap, known: false, from: 'smiles'};
      }
      if (canonicalize !== undefined) {
        const cap = {name: name || altName || fragment, description: 'as defined',
          smiles: smiles.replace(ATTACHMENT, '[*]')};
        return {cap, known: false, from: 'smiles', note: 'the SMILES is not valid'};
      }
    }
    if (named !== null)
      return {cap: named, known: true, from: CapGroups.byName(name) !== null ? 'name' : 'alternateId'};
    const own = name || altName;
    if (own === '')
      return null;
    return {cap: {name: own, description: 'as the library defines it', smiles: smiles.replace(ATTACHMENT, '[*]')},
      known: false, from: name ? 'name' : 'alternateId'};
  }

  private static _byCanonical(canonicalize: Canonicalize): Map<string, CapGroup> {
    if (CapGroups._canonical === undefined) {
      const map = new Map<string, CapGroup>();
      for (const cap of COMMON_CAP_GROUPS) {
        const canonical = canonicalize(cap.smiles);
        if (canonical !== null)
          map.set(canonical, cap);
      }
      CapGroups._canonical = map;
    }
    return CapGroups._canonical;
  }
}
