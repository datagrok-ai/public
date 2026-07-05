/* eslint-disable max-len */
/* eslint-disable camelcase */
import {RDModule, RDMol, RDReaction, MolList} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {BRANCH_DELIMITER} from '../../rendering/rdkit-reaction-renderer';
import {EnumeratorConfig} from './config';
import {applyProductFilters, buildExclusionQmols, computeMolStats, MolStats, stripIsotopesFromSmiles} from './filters';

export interface TemplateInput {
  smarts: string;
  blockingSmartsList: string[];
  reactionName: string;
}

export interface RouteStep {
  reactants: string[];
  product: string;
  templateSmarts: string;
  reactionName: string;
}

export type Route = RouteStep[];

export interface ProductRecord {
  smiles: string;
  routes: Route[];
  firstRound: number;
}

export interface EnumerationResult {
  productsByRound: ProductRecord[][];
  warnings: string[];
}

export interface EnumerationProgress {
  round: number;
  numRounds: number;
  templateIndex: number;
  numTemplates: number;
  productsSoFar: number;
  combosDone?: number;
  combosTotal?: number;
}

// Per-round narrowing of the global pools, indexed by round-1. An undefined field falls back to the
// global list for that round. `templates` is matched against the global list by SMARTS string.
export interface PerRoundOverride {
  templates?: TemplateInput[];
  buildingBlocks?: string[];
  reagents?: string[];
}

export interface EnumerateOptions {
  rdkit: RDModule;
  config: EnumeratorConfig;
  templates: TemplateInput[];
  buildingBlocks: string[];
  exclusionSmarts: string[];
  // Optional library of "reagents". When non-empty, switches to reagents mode: every step uses
  // EXACTLY ONE reactant from the main pool (BBs in round 1, round-(R-1) products in round R > 1)
  // and fills every remaining slot from this reagents library. Used to grow derivatives of each
  // BB across rounds (P1 = BB + reagents, P2 = P1 + reagents, …). Overrides depth_first /
  // breadth_first slot logic when present.
  reagents?: string[];
  // Round r (1-based) narrows templates/BBs/reagents to perRoundOverrides[r-1]; reagents MODE itself
  // stays global (driven by `reagents` presence) — a per-round `reagents` list only narrows the pool.
  perRoundOverrides?: PerRoundOverride[];
  onProgress?: (p: EnumerationProgress) => void;
  isCancelled?: () => boolean;
}

export function splitSmartsByReactants(lhs: string): string[] {
  const parts: string[] = [];
  let depth = 0;
  let start = 0;
  for (let i = 0; i < lhs.length; i++) {
    const c = lhs[i];
    if (c === '(' || c === '[') depth++;
    else if (c === ')' || c === ']') depth--;
    else if (c === '.' && depth === 0) {
      parts.push(lhs.slice(start, i));
      start = i + 1;
    }
  }
  parts.push(lhs.slice(start));
  return parts.map((p) => p.trim()).filter((p) => p.length > 0);
}

// A reaction-SMARTS reactant slot may be wrapped in `(...)` to denote a tied group of fragments
// that must come from a single connected molecule (multi-fragment SMARTS). Strip the outer parens
// so the inner SMARTS can be parsed by get_qmol.
export function stripOuterParens(s: string): string {
  if (s.length < 2 || s[0] !== '(' || s[s.length - 1] !== ')') return s;
  let depth = 0;
  for (let i = 0; i < s.length; i++) {
    if (s[i] === '(' || s[i] === '[') depth++;
    else if (s[i] === ')' || s[i] === ']') depth--;
    if (depth === 0 && i < s.length - 1) return s;
  }
  return s.slice(1, -1).trim();
}

function tryGetMol(rdkit: RDModule, smiles: string): RDMol | null {
  try {
    const m = rdkit.get_mol(smiles);
    if (!m) return null;
    if (!m.is_valid()) {m.delete(); return null;}
    return m;
  } catch {
    return null;
  }
}

function tryGetQmol(rdkit: RDModule, smarts: string): RDMol | null {
  try {
    const q = rdkit.get_qmol(smarts);
    if (!q) return null;
    if (!q.is_valid()) {q.delete(); return null;}
    return q;
  } catch {
    return null;
  }
}

function tryGetRxn(rdkit: RDModule, smarts: string): RDReaction | null {
  try {
    const r = rdkit.get_rxn(smarts);
    return r ?? null;
  } catch {
    return null;
  }
}

interface ParsedTemplate {
  index: number;
  smarts: string;
  reactionName: string;
  rxn: RDReaction;
  // [slot][fragment] — tied groups (parenthesized in the SMARTS) yield multi-fragment slots.
  // A BB matches a slot only if it matches every fragment.
  reactantQmols: RDMol[][];
  blockingQmols: RDMol[];
  numReactants: number;
}

function getSlotFragments(rs: string): string[] {
  const stripped = stripOuterParens(rs);
  return splitSmartsByReactants(stripped);
}

function parseTemplate(rdkit: RDModule, t: TemplateInput, idx: number, warnings: string[]): ParsedTemplate | null {
  const arrowIdx = t.smarts.indexOf('>>');
  if (arrowIdx < 0) {
    warnings.push(`Template ${idx + 1} (${t.reactionName || 'unnamed'}) has no '>>' separator; skipping.`);
    return null;
  }
  const lhs = t.smarts.slice(0, arrowIdx);
  const slotSmarts = splitSmartsByReactants(lhs);
  const rxn = tryGetRxn(rdkit, t.smarts);
  if (!rxn) {
    warnings.push(`Template ${idx + 1} (${t.reactionName || 'unnamed'}): get_rxn failed; skipping.`);
    return null;
  }
  const reactantQmols: RDMol[][] = [];
  const blockingQmols: RDMol[] = [];
  for (const slot of slotSmarts) {
    const fragments = getSlotFragments(slot);
    const fragQmols: RDMol[] = [];
    let ok = true;
    for (const frag of fragments) {
      const q = tryGetQmol(rdkit, frag);
      if (!q) {
        warnings.push(`Template ${idx + 1} (${t.reactionName || 'unnamed'}): invalid reactant SMARTS fragment '${frag}'; skipping template.`);
        ok = false;
        break;
      }
      fragQmols.push(q);
    }
    if (!ok) {
      rxn.delete();
      for (const slotQs of reactantQmols) for (const q of slotQs) q.delete();
      for (const q of fragQmols) q.delete();
      return null;
    }
    reactantQmols.push(fragQmols);
  }
  for (const bs of t.blockingSmartsList) {
    const trimmed = (bs ?? '').trim();
    if (!trimmed) continue;
    const q = tryGetQmol(rdkit, trimmed);
    if (q) blockingQmols.push(q);
  }
  return {
    index: idx, smarts: t.smarts, reactionName: t.reactionName,
    rxn, reactantQmols, blockingQmols, numReactants: slotSmarts.length,
  };
}

function disposeTemplate(t: ParsedTemplate): void {
  if ((t as any)._disposed) return;
  (t as any)._disposed = true;
  try {t.rxn.delete();} catch (e) {
    console.warn(`Template ${t.index + 1} rxn delete failed: ${e instanceof Error ? e.message : String(e)}`);
  }
  for (let i = 0; i < t.reactantQmols.length; i++) {
    const slotQs = t.reactantQmols[i];
    for (let j = 0; j < slotQs.length; j++) {
      try {slotQs[j].delete();} catch (e) {
        console.warn(`Template ${t.index + 1} reactant qmol [${i}][${j}] delete failed: ${e instanceof Error ? e.message : String(e)}`);
      }
    }
    slotQs.length = 0;
  }
  t.reactantQmols.length = 0;
  for (let i = 0; i < t.blockingQmols.length; i++) {
    try {t.blockingQmols[i].delete();} catch (e) {
      console.warn(`Template ${t.index + 1} blocking qmol [${i}] delete failed: ${e instanceof Error ? e.message : String(e)}`);
    }
  }
  t.blockingQmols.length = 0;
}

function bbMatchesSlot(mol: RDMol, fragQmols: RDMol[]): boolean {
  for (const q of fragQmols) if (!hasMatch(mol, q)) return false;
  return true;
}

function hasMatch(mol: RDMol, qmol: RDMol): boolean {
  let raw: string;
  try {
    raw = mol.get_substruct_match(qmol);
  } catch {
    return false;
  }
  if (!raw || raw === '' || raw === '{}') return false;
  try {
    const obj = JSON.parse(raw);
    return Array.isArray(obj.atoms) && obj.atoms.length > 0;
  } catch {
    return raw.length > 2;
  }
}

class MolCache {
  private map = new Map<string, RDMol>();
  constructor(private rdkit: RDModule) {}
  get(smiles: string): RDMol | null {
    const cached = this.map.get(smiles);
    if (cached) return cached;
    const m = tryGetMol(this.rdkit, smiles);
    if (!m) return null;
    this.map.set(smiles, m);
    return m;
  }
  dispose(): void {
    for (const m of this.map.values())
      try {m.delete();} catch {/* ignore */}

    this.map.clear();
  }
}

function canonicalize(rdkit: RDModule, smiles: string): string | null {
  const mol = tryGetMol(rdkit, smiles);
  if (!mol) return null;
  try {
    return mol.get_smiles();
  } catch {
    return null;
  } finally {
    try {mol.delete();} catch {/* ignore */}
  }
}

function* cartesian<T>(slots: T[][]): Generator<T[]> {
  const n = slots.length;
  if (n === 0) {yield []; return;}
  const idx = new Array(n).fill(0);
  const lengths = slots.map((s) => s.length);
  if (lengths.some((l) => l === 0)) return;
  while (true) {
    yield idx.map((i, k) => slots[k][i]);
    let k = n - 1;
    while (k >= 0) {idx[k]++; if (idx[k] < lengths[k]) break; idx[k] = 0; k--;}
    if (k < 0) return;
  }
}

// Format a route as a single multi-step reaction string parseable by the Chem package's
// reaction renderer (parseMultiStepReaction in rdkit-reaction-renderer.ts).
//
// Each step is emitted as its own self-contained `reactants>>product` segment, joined by the
// renderer's BRANCH_DELIMITER (`--**--`). The renderer splits on the delimiter so each step is
// recovered exactly as it was, with no fake "prev_product → next_reactants" intermediates at
// the seams. If a step uses a product from a previous step (linear chain, depth-first), that
// product simply appears in this step's reactants list — no elision.
//
// Example route:   bb1+bb2→p1, bb3+bb4→p2, p2+bb5→p3, p3+p1→p4
//   formatted:     "bb1.bb2>>p1--**--bb3.bb4>>p2--**--p2.bb5>>p3--**--p3.p1>>p4"
//   renderer sees: 4 real steps in order, no fake seams.
export function formatRoute(route: Route): string {
  if (route.length === 0) return '';
  return route.map((s) => `${s.reactants.join('.')}>>${s.product}`).join(BRANCH_DELIMITER);
}

export interface OutputRow {
  product: string;
  route: string;
  template: string;
  reaction_name: string;
  round: number;
  n_routes: number;
}

export async function enumerate(opts: EnumerateOptions): Promise<{rows: OutputRow[]; warnings: string[]}> {
  const {rdkit, config, templates, buildingBlocks, exclusionSmarts, reagents,
    perRoundOverrides, onProgress, isCancelled} = opts;
  const warnings: string[] = [];

  // Shared canonicalize + dedup, used for both the global BB/reagent pools and per-round overrides,
  // so an override list matches the canonical SMILES used in the product pools.
  const canonUnique = (list: string[], label: string): string[] => {
    const out: string[] = [];
    for (const s of list) {
      const c = canonicalize(rdkit, s);
      if (c) out.push(c);
      else warnings.push(`Skipped invalid ${label} SMILES: ${s}`);
    }
    return Array.from(new Set(out));
  };

  const parsedTemplates: ParsedTemplate[] = [];
  for (let i = 0; i < templates.length; i++) {
    const p = parseTemplate(rdkit, templates[i], i, warnings);
    if (p) parsedTemplates.push(p);
  }
  if (parsedTemplates.length === 0) {
    warnings.push('No valid templates found.');
    return {rows: [], warnings};
  }

  const exclusion = buildExclusionQmols(rdkit, exclusionSmarts);

  const uniqueBBs = canonUnique(buildingBlocks, 'BB');

  const reagentsMode = !!(reagents && reagents.length > 0);
  const uniqueReagents: string[] = reagentsMode ? canonUnique(reagents!, 'reagent') : [];
  if (reagentsMode && uniqueReagents.length === 0)
    warnings.push('Reagents mode: no valid reagents after canonicalization; falling back to BB-only mode.');
  const useReagents = reagentsMode && uniqueReagents.length > 0;

  // Precompute per-round override sets, indexed by round-1. An empty/non-matching override falls
  // back to the global pool (rather than zeroing the round silently) and warns.
  const roundAllowedSmarts: (Set<string> | null)[] = [];
  const roundBBs: (string[] | null)[] = [];
  const roundReagents: (string[] | null)[] = [];
  for (let r = 0; r < config.enumeration.num_rounds; r++) {
    const o = perRoundOverrides?.[r];
    if (o?.templates && o.templates.length > 0) {
      const allowed = new Set(o.templates.map((t) => t.smarts));
      const intersects = parsedTemplates.some((t) => allowed.has(t.smarts));
      if (intersects)
        roundAllowedSmarts.push(allowed);
      else {
        warnings.push(`Round ${r + 1}: template override matches none of the global reaction templates; ` +
          `using the global template set instead.`);
        roundAllowedSmarts.push(null);
      }
    } else {
      if (o?.templates) warnings.push(`Round ${r + 1}: empty template override; using the global template set instead.`);
      roundAllowedSmarts.push(null);
    }
    const bbs = o?.buildingBlocks ? canonUnique(o.buildingBlocks, 'round BB') : null;
    roundBBs.push(bbs && bbs.length > 0 ? bbs : null);
    const rgs = o?.reagents ? canonUnique(o.reagents, 'round reagent') : null;
    roundReagents.push(rgs && rgs.length > 0 ? rgs : null);
  }

  const molCache = new MolCache(rdkit);

  // Round 0 ("keep building blocks in output") must mirror round 1's ACTIVE BB pool, but only in the
  // modes where a round-1 override actually narrows it (useReagents / depth_first — see
  // `eligibleSmiles` below). In plain breadth-first mode the override is a no-op, so seeding round 0
  // from it anyway would both misreport the kept-BB output AND leak the override into round 1's own
  // pool via `allPriorPool` (which folds in every prior round, including round 0).
  const round1OverrideActive = useReagents || config.enumeration.depth_first;
  const productPools: ProductRecord[][] = [];
  productPools.push((round1OverrideActive ? (roundBBs[0] ?? uniqueBBs) : uniqueBBs)
    .map((s) => ({smiles: s, routes: [], firstRound: 0})));

  const {max_num_components, max_num_combinations_per_template, max_num_routes_per_compound,
    max_num_routes_per_compound: _maxRoutes} = config;
  void _maxRoutes;

  const buildBlockingFilter = (mol: RDMol, blocking: RDMol[]): boolean => {
    for (const q of blocking) if (hasMatch(mol, q)) return true;
    return false;
  };

  // Time-based yield: surrender to the event loop when we've been computing continuously for more
  // than YIELD_INTERVAL_MS. Use setTimeout(0), NOT requestAnimationFrame — RAF is fully suspended
  // (not just throttled) in a backgrounded tab, which would hang the run (and Cancel with it).
  const YIELD_INTERVAL_MS = 50;
  let lastYield = performance.now();
  type ProgressCtx = {round: number; ti: number; total: number; combosTotal: number; combosDone: number; productsSoFar: number};
  let progressContext: ProgressCtx | null = null;
  const yieldIfNeeded = async (): Promise<void> => {
    const now = performance.now();
    if (now - lastYield < YIELD_INTERVAL_MS) return;
    lastYield = now;
    if (progressContext && onProgress) {
      onProgress({
        round: progressContext.round, numRounds: config.enumeration.num_rounds,
        templateIndex: progressContext.ti, numTemplates: progressContext.total,
        productsSoFar: progressContext.productsSoFar,
        combosDone: progressContext.combosDone, combosTotal: progressContext.combosTotal,
      });
    }
    await new Promise<void>((resolve) => setTimeout(resolve, 0));
    lastYield = performance.now();
  };

  try {
    for (let round = 1; round <= config.enumeration.num_rounds; round++) {
      if (isCancelled?.()) break;

      const newPool = new Map<string, ProductRecord>();
      const prevRoundProducts = round === 1 ? [] : productPools[round - 1].map((p) => p.smiles);
      const allPriorPool = new Set<string>();
      for (let r = 0; r < round; r++) for (const p of productPools[r]) allPriorPool.add(p.smiles);

      // Per-round narrowing; each falls back to the global list when this round has no override.
      const allowedSmarts = roundAllowedSmarts[round - 1];
      const activeBBs = roundBBs[round - 1] ?? uniqueBBs;
      const activeReagents = roundReagents[round - 1] ?? uniqueReagents;

      // Pool of SMILES that can fill any reactant slot this round.
      // - depth_first round 1: only original BBs (BBs combine with BBs).
      // - depth_first round R > 1: BBs ∪ R(R-1) products. The combo filter below then enforces
      //   exactly one R(R-1) product + (N-1) BBs per combo (linear chain extension).
      // - breadth_first: any product produced so far (rounds 0..R-1).
      // - reagents mode reuses the depth_first pool for the "main" slot only; other slots draw
      //   from the reagents library and the depth-first combo filter is skipped (the slot pools
      //   already enforce "exactly one from main, rest from reagents").
      const eligibleSmiles = useReagents || config.enumeration.depth_first ?
        (round === 1 ? activeBBs : Array.from(new Set([...activeBBs, ...prevRoundProducts]))) :
        Array.from(allPriorPool);
      // A per-round BB override is a no-op in breadth-first mode (eligibleSmiles above never reads
      // activeBBs there) — surface that at run time, not just via the UI's dot indicator.
      if (!useReagents && !config.enumeration.depth_first && roundBBs[round - 1] != null) {
        warnings.push(`Round ${round}: building-block override has no effect in breadth-first ` +
          `mode — a round draws from all earlier products regardless of the per-step BB subset.`);
      }

      for (let ti = 0; ti < parsedTemplates.length; ti++) {
        if (isCancelled?.()) break;
        const t = parsedTemplates[ti];
        if (t.numReactants > max_num_components) continue;
        // Per-step template subset: skip any global template not chosen for this round.
        if (allowedSmarts && !allowedSmarts.has(t.smarts)) continue;

        progressContext = {
          round, ti, total: parsedTemplates.length,
          combosTotal: 0, combosDone: 0,
          productsSoFar: newPool.size,
        };
        onProgress?.({round, numRounds: config.enumeration.num_rounds,
          templateIndex: ti, numTemplates: parsedTemplates.length,
          productsSoFar: newPool.size});

        // Build one or more slot configurations to enumerate. Each configuration is a
        // slots[i] = SMILES candidates for reactant slot i; cartesian(slots) generates combos.
        //
        // - Standard mode: a single configuration where every slot draws from eligibleSmiles.
        // - Reagents mode: one configuration per "main slot" index. For each i in 0..N-1,
        //   slot i draws from eligibleSmiles (BBs ∪ round-(R-1) products) and every other
        //   slot draws from the reagents library. Different "main slot" choices can produce
        //   duplicate combos for symmetric templates; the product-pool dedup handles that.
        const slotConfigs: string[][][] = [];

        if (useReagents) {
          for (let mainSlot = 0; mainSlot < t.numReactants; mainSlot++) {
            if (isCancelled?.()) break;
            const slots: string[][] = [];
            for (let i = 0; i < t.numReactants; i++) slots.push([]);

            for (const smi of eligibleSmiles) {
              if (isCancelled?.()) break;
              await yieldIfNeeded();
              const mol = molCache.get(smi);
              if (!mol) continue;
              if (buildBlockingFilter(mol, t.blockingQmols)) continue;
              if (bbMatchesSlot(mol, t.reactantQmols[mainSlot])) slots[mainSlot].push(smi);
            }

            if (slots[mainSlot].length === 0) continue;

            let otherSlotsOk = true;
            for (let i = 0; i < t.numReactants; i++) {
              if (i === mainSlot) continue;
              for (const smi of activeReagents) {
                if (isCancelled?.()) break;
                await yieldIfNeeded();
                const mol = molCache.get(smi);
                if (!mol) continue;
                if (buildBlockingFilter(mol, t.blockingQmols)) continue;
                if (bbMatchesSlot(mol, t.reactantQmols[i])) slots[i].push(smi);
              }
              if (slots[i].length === 0) {otherSlotsOk = false; break;}
            }
            if (otherSlotsOk) slotConfigs.push(slots);
          }
        } else {
          const slots: string[][] = [];
          for (let i = 0; i < t.numReactants; i++) slots.push([]);
          for (const smi of eligibleSmiles) {
            if (isCancelled?.()) break;
            await yieldIfNeeded();
            const mol = molCache.get(smi);
            if (!mol) continue;
            if (buildBlockingFilter(mol, t.blockingQmols)) continue;
            for (let i = 0; i < t.numReactants; i++)
              if (bbMatchesSlot(mol, t.reactantQmols[i])) slots[i].push(smi);
          }
          if (!slots.some((s) => s.length === 0)) slotConfigs.push(slots);
        }

        if (slotConfigs.length === 0) continue;

        let totalCombos = 0;
        for (const cfg of slotConfigs) {
          let c = 1;
          for (const s of cfg) c *= s.length;
          totalCombos += c;
        }
        const comboCap = max_num_combinations_per_template;
        let executed = 0;
        let truncated = false;
        progressContext.combosTotal = totalCombos;

        configLoop:
        for (const slots of slotConfigs) {
          for (const combo of cartesian(slots)) {
            if (isCancelled?.()) break configLoop;
            if (comboCap >= 0 && executed >= comboCap) {
              truncated = true;
              break configLoop;
            }
            executed++;
            progressContext.combosDone = executed;
            progressContext.productsSoFar = newPool.size;
            await yieldIfNeeded();

            // In reagents mode the slot pools already enforce "exactly one main + reagents
            // elsewhere", so skip the depth-first combo filter (it would reject everything
            // because reagents aren't in the BB or prev-round sets).
            if (!useReagents && config.enumeration.depth_first && round > 1) {
            // Strict depth-first = linear chain extension. Each step takes EXACTLY ONE
            // round-(R-1) product and reacts it with original BBs only — no merging two
            // complex products in a single step (that would be convergent/breadth-first).
              const prevSet = new Set(prevRoundProducts);
              const bbSet = new Set(activeBBs);
              let prevCount = 0;
              let bbCount = 0;
              for (const c of combo) {
                if (prevSet.has(c)) prevCount++;
                else if (bbSet.has(c)) bbCount++;
              }
              if (prevCount !== 1 || prevCount + bbCount !== combo.length) continue;
            }

            let molList: MolList | null = null;
            let result: ReturnType<RDReaction['run_reactants']> | null = null;
            // Always create fresh mols from SMILES for the MolList — never share with the BB cache.
            // copy() appears to be shallow at the WASM level: a copy still aliases the original's
            // managed memory, and molList.delete() reaches into both, corrupting the heap and
            // surfacing as "null function" errors during dispose. Fresh-parse-per-combo matches the
            // Chem package's reactions.ts pattern.
            const inputMols: RDMol[] = [];
            try {
              molList = new rdkit.MolList();
              let allValid = true;
              for (const smi of combo) {
                const fresh = tryGetMol(rdkit, smi);
                if (!fresh) {allValid = false; break;}
                inputMols.push(fresh);
                molList.append(fresh);
              }
              if (!allValid) continue;

              try {
              // FIX 1: Changed cap from 1 → 0 (unlimited) so reactions that fire at
              // multiple sites return all products instead of only the first one.
                result = t.rxn.run_reactants(molList, 0);
              } catch (e) {
                warnings.push(`run_reactants failed for template ${t.index + 1}: ${e instanceof Error ? e.message : String(e)}`);
                continue;
              }
              if (!result || result.size() === 0) continue;

              // Match the example in rdkit-api.ts exactly: result.get(i).next() once, no reset/at_end.
              // Calling reset() or at_end() before next() can leave the iterator in a corrupt state.
              // FIX 2: Was result.get(0) — only the first product set was ever processed.
              // Now loop over ALL result sets so every reaction site's product is collected.
              const producedSmilesSet = new Set<string>();
              for (let ri = 0; ri < result.size(); ri++) {
                const productSet = result.get(ri);
                try {
                  const productMol = productSet.next();
                  try {
                    if (!productMol || !productMol.is_valid()) continue;

                    let productSmiles: string;
                    try {productSmiles = productMol.get_smiles();} catch {continue;}
                    if (config.products_specs.remove_isotope_information && /\[\d+[A-Za-z]/.test(productSmiles)) {
                      const stripped = stripIsotopesFromSmiles(productSmiles);
                      const re = canonicalize(rdkit, stripped);
                      if (!re) continue;
                      productSmiles = re;
                    }
                    if (producedSmilesSet.has(productSmiles))
                      continue;
                    producedSmilesSet.add(productSmiles);

                    // Always evaluate against a fresh mol parsed from the canonical SMILES — never touch
                    // productMol directly (it can be in a fragile post-reaction state, and downstream
                    // ops like get_substruct_match on it appear to corrupt the WASM heap).
                    const evalMol = tryGetMol(rdkit, productSmiles);
                    if (!evalMol) continue;
                    try {
                      let stats: MolStats;
                      try {
                        stats = computeMolStats(evalMol);
                      } catch (e) {
                        warnings.push(`computeMolStats failed: ${e instanceof Error ? e.message : String(e)}`);
                        continue;
                      }
                      const fr = applyProductFilters(stats, config.products_specs, exclusion.qmols, evalMol);
                      if (!fr.pass) continue;
                    } finally {
                      try {evalMol.delete();} catch {/* ignore */}
                    }

                    const step: RouteStep = {
                      reactants: combo.slice(),
                      product: productSmiles,
                      templateSmarts: t.smarts,
                      reactionName: t.reactionName,
                    };

                    // Build the base route(s) — the synthesis history that must precede this step.
                    // For each combo component that's a previously-synthesized product (from ANY
                    // earlier round, not just round-(r-1)), splice in its known route. BBs and any
                    // unknown reactant are treated as starting materials with no prior history.
                    //
                    // Looking only at round-(r-1) was the original behavior, but it dropped the
                    // synthesis steps of round-(r-2) (or earlier) products that appear in the combo.
                    // In breadth-first mode (or whenever a step legitimately consumes products from
                    // multiple earlier rounds), the resulting route would look like it started from
                    // pre-formed peptides instead of bare BBs. Now we look across all prior rounds.
                    let baseRoutes: Route[];
                    if (round === 1)
                      baseRoutes = [[]];
                    else {
                      const allPrevProducts = new Map<string, ProductRecord>();
                      for (let r = 1; r < round; r++) {
                        for (const p of productPools[r]) {
                          // First occurrence wins (typically the shorter route).
                          if (!allPrevProducts.has(p.smiles)) allPrevProducts.set(p.smiles, p);
                        }
                      }
                      const prevComponents = combo.filter((c) => allPrevProducts.has(c));
                      if (prevComponents.length === 0) baseRoutes = [[]];
                      else {
                        const prevRouteLists = prevComponents.map((pc) => {
                          const rec = allPrevProducts.get(pc);
                          return rec && rec.routes.length > 0 ? rec.routes : [[]];
                        });
                        baseRoutes = [];
                        for (const combo2 of cartesian(prevRouteLists)) {
                          const merged: Route = [];
                          for (const r of combo2) for (const s of r) merged.push(s);
                          baseRoutes.push(merged);
                        }
                      }
                    }

                    let rec = newPool.get(productSmiles);
                    if (!rec) {
                      rec = {smiles: productSmiles, routes: [], firstRound: round};
                      newPool.set(productSmiles, rec);
                    }
                    for (const base of baseRoutes) {
                      if (max_num_routes_per_compound >= 0 && rec.routes.length >= max_num_routes_per_compound) break;
                      rec.routes.push([...base, step]);
                    }
                  } finally {
                    try {productMol?.delete();} catch {/* ignore */}
                  }
                } finally {
                  try {productSet.delete();} catch {/* ignore */}
                }
              } // end FIX 2: for ri (product sets)
            } catch (e) {
              warnings.push(`Combo execution failed for template ${t.index + 1}: ${e instanceof Error ? e.message : String(e)}`);
            } finally {
              try {result?.delete();} catch {/* ignore */}
              try {molList?.delete();} catch {/* ignore */}
              for (const m of inputMols) try {m.delete();} catch {/* ignore */}
              inputMols.length = 0;
            }
          }
        }
        if (truncated)
          warnings.push(`Template ${t.index + 1} (${t.reactionName || ''}): truncated at ${executed}/${totalCombos} combinations (cap=${comboCap}).`);
      }

      productPools.push(Array.from(newPool.values()));
    }
  } catch (e) {
    console.error('Enumeration failed:', e instanceof Error ? e.message : String(e));
  } finally {
    for (const t of parsedTemplates) disposeTemplate(t);
    exclusion.dispose();
    molCache.dispose();
  }

  const finalProducts = new Map<string, ProductRecord>();
  const startRound = config.keep_building_blocks_in_final_output ? 0 : 1;
  for (let r = startRound; r < productPools.length; r++) {
    for (const p of productPools[r]) {
      const ex = finalProducts.get(p.smiles);
      if (!ex) finalProducts.set(p.smiles, {...p, routes: p.routes.slice()});
      else for (const route of p.routes) ex.routes.push(route);
    }
  }

  const rows: OutputRow[] = [];
  for (const rec of finalProducts.values()) {
    if (rec.routes.length === 0) {
      rows.push({product: rec.smiles, route: '', template: '', reaction_name: '',
        round: rec.firstRound, n_routes: 0});
    } else {
      const cap = max_num_routes_per_compound;
      const limited = cap >= 0 ? rec.routes.slice(0, cap) : rec.routes;
      for (const route of limited) {
        const last = route[route.length - 1];
        rows.push({
          product: rec.smiles,
          route: formatRoute(route),
          template: last.templateSmarts,
          reaction_name: last.reactionName,
          round: rec.firstRound,
          n_routes: rec.routes.length,
        });
      }
    }
  }

  return {rows, warnings};
}
