/* eslint-disable max-len */
import {cleanupHelmSymbol} from '@datagrok-libraries/bio/src/helm/utils';

import {PolyToolEnumeratorParams, PolyToolEnumeratorTypes,
  PolyToolPlaceholder, PolyToolBreadthPlaceholder} from './types';

// For example keep monomers presented in HELMCoreLibrary.json only (not [NH2])
export const PT_HELM_EXAMPLE = 'PEPTIDE1{R.[Aca].T.G.H.F.G.A.A.Y.P.E.[meI]}$$$$';

// hwe migration (Phase 6): enumeration used to parse the HELM into a mutable
// JSDraw2 / HELM Web Editor mol (`org.helm.webeditor.IO.parseHelm`), clone it,
// overwrite `atom.elem` at a position and re-serialize (`...IO.getHelm`). The
// net effect on a peptide HELM is simply swapping the monomer token at a given
// 0-based position. The standalone hwe model is immutable, so instead of
// rebuilding a mutable mol we perform the swap directly on the HELM string —
// bracket-aware so inline-SMILES monomers (`[... |$_R1;..$|]`) and multi-char
// symbols are tokenized correctly. This drops the `JSDraw2`/`org` globals and,
// because untouched positions are copied verbatim, preserves the input's exact
// formatting (e.g. the `$$$$V2.0` / `[Ac(1)]` suffix the tests assert).

/** One monomer token's `[start, end)` char span inside the HELM string. */
type MonomerSpan = {start: number, end: number};

/** A single enumerated variant: position→new-symbol swaps + accumulated name. */
type EnumState = {swaps: Map<number, string>, name: string};

/** HELM-id form of a monomer symbol: multi-char symbols are bracketed. */
function idSymbol(sym: string): string {
  return sym != null && sym.length > 1 ? `[${sym}]` : sym;
}

/**
 * Bracket-aware extraction of every monomer-token span across all polymer
 * bodies (`PEPTIDE1{m.m.m}|PEPTIDE2{m.m}`), in source order — the same order
 * the legacy parser laid atoms out for a peptide. `[`/`]` depth shields the
 * `.`/`}`/`$` that appear inside inline-SMILES monomers.
 */
function getMonomerSpans(helm: string): MonomerSpan[] {
  const spans: MonomerSpan[] = [];
  let depth = 0;
  let inBody = false;
  let tokenStart = -1;
  for (let i = 0; i < helm.length; i++) {
    const c = helm[i];
    if (c === '[') depth++;
    else if (c === ']') depth--;

    if (depth === 0 && c === '$' && !inBody)
      break; // end of the polymer section (section-1 separator)

    if (!inBody) {
      if (depth === 0 && c === '{') {
        inBody = true;
        tokenStart = i + 1;
      }
    } else {
      if (depth === 0 && c === '.') {
        spans.push({start: tokenStart, end: i});
        tokenStart = i + 1;
      } else if (depth === 0 && c === '}') {
        spans.push({start: tokenStart, end: i});
        inBody = false;
      }
    }
  }
  return spans;
}

/** Materialize one variant: copy the source HELM, replacing swapped tokens. */
function materialize(helm: string, spans: MonomerSpan[], state: EnumState): [string, string] {
  if (state.swaps.size === 0)
    return [helm, state.name];

  let res = '';
  let cursor = 0;
  for (let pos = 0; pos < spans.length; pos++) {
    if (state.swaps.has(pos)) {
      res += helm.slice(cursor, spans[pos].start);
      res += idSymbol(state.swaps.get(pos)!);
      cursor = spans[pos].end;
    }
  }
  res += helm.slice(cursor);
  return [res, state.name];
}

/**
 * Swap each position in `[start, end]` with each monomer in `monomers`, chained
 * onto every input state (monomer outer, position inner — matching the legacy
 * `polyToolEnumeratorCore` loop nesting). `oldSymbols` carries the cleaned
 * source symbol per position for the name diff.
 */
function coreRange(
  states: EnumState[], start: number, end: number, monomers: string[], oldSymbols: string[]
): EnumState[] {
  const res: EnumState[] = [];
  for (const st of states) {
    for (let monI = 0; monI < monomers.length; monI++) {
      for (let posI = 0; posI <= end - start; posI++) {
        const pos = start + posI;
        const newSymbol = monomers[monI];
        // The old symbol is the CURRENT value at the position (the legacy code
        // read `atoms[pos].elem`), so re-swapping a position diffs against the
        // value the previous round wrote, not the original.
        const oldSymbol = st.swaps.has(pos) ? st.swaps.get(pos)! : oldSymbols[pos];
        const swaps = new Map(st.swaps);
        swaps.set(pos, newSymbol);
        const name = `${st.name}-${idSymbol(oldSymbol)}${pos + 1}${idSymbol(newSymbol)}`;
        res.push({swaps, name});
      }
    }
  }
  return res;
}

function getPtEnumeratorSingle(
  initial: EnumState, placeholders: PolyToolPlaceholder[], oldSymbols: string[]
): EnumState[] {
  return placeholders
    .map((ph) => coreRange([initial], ph.position, ph.position, ph.monomers, oldSymbols))
    .reduce((acc, l) => acc.concat(l), []);
}

function getPtEnumeratorMatrix(
  initial: EnumState, placeholders: PolyToolPlaceholder[], oldSymbols: string[]
): EnumState[] {
  let states = [initial];
  for (const ph of placeholders)
    states = coreRange(states, ph.position, ph.position, ph.monomers, oldSymbols);
  return states;
}

/** Parallel (zip) enumeration: the i-th result takes the i-th monomer from each placeholder position.
 * All placeholders must have the same number of monomers (validated upstream).
 * With K positions and N monomers each, produces exactly N results. */
function getPtEnumeratorParallel(
  initial: EnumState, placeholders: PolyToolPlaceholder[], oldSymbols: string[]
): EnumState[] {
  if (placeholders.length === 0)
    return [];

  const monomerCount = placeholders[0].monomers.length;
  for (const ph of placeholders) {
    if (ph.monomers.length !== monomerCount)
      throw new Error(`Parallel enumeration requires all positions to have the same number of monomers`);
  }

  const res: EnumState[] = new Array<EnumState>(monomerCount);
  for (let i = 0; i < monomerCount; i++) {
    const swaps = new Map(initial.swaps);
    let name = initial.name;
    for (const ph of placeholders) {
      const pos = ph.position;
      const newSymbol = ph.monomers[i];
      const oldSymbol = swaps.has(pos) ? swaps.get(pos)! : oldSymbols[pos];
      name += `-${idSymbol(oldSymbol)}${pos + 1}${idSymbol(newSymbol)}`;
      swaps.set(pos, newSymbol);
    }
    res[i] = {swaps, name};
  }
  return res;
}

function getPtEnumeratorBreadth(
  initial: EnumState, placeholdersBreadth: PolyToolBreadthPlaceholder[], oldSymbols: string[]
): EnumState[] {
  if (placeholdersBreadth.length == 0)
    return [];

  let states = [initial];
  for (const phb of placeholdersBreadth)
    states = coreRange(states, phb.start, phb.end, phb.monomers, oldSymbols);
  return states;
}

/**
 * @param {string} helm Molecule string in HELM format
 * @param {string} id Base name for the enumerated variants
 * @param  params Enumeration parameters (placeholders by 0-based position)
 * @returns {[string, string][]} List of `[helm, id]` enumerated molecules. Covered with tests.
 */
export function doPolyToolEnumerateHelm(
  helm: string, id: string, params: PolyToolEnumeratorParams
): [ /* helm */ string, /* id */ string][] {
  const spans = getMonomerSpans(helm);
  const oldSymbols = spans.map((s) => cleanupHelmSymbol(helm.slice(s.start, s.end)));
  const initial: EnumState = {swaps: new Map<number, string>(), name: id};

  let resStates: EnumState[] = [];
  if (params.placeholders) {
    switch (params.type) {
    case PolyToolEnumeratorTypes.Single: {
      resStates = getPtEnumeratorSingle(initial, params.placeholders, oldSymbols);
      break;
    }
    case PolyToolEnumeratorTypes.Parallel: {
      resStates = getPtEnumeratorParallel(initial, params.placeholders, oldSymbols);
      break;
    }
    case PolyToolEnumeratorTypes.Matrix: {
      resStates = getPtEnumeratorMatrix(initial, params.placeholders, oldSymbols);
      break;
    }
    }
  }

  if (params.breadthPlaceholders)
    resStates = resStates.concat(getPtEnumeratorBreadth(initial, params.breadthPlaceholders, oldSymbols));

  if (params.keepOriginal)
    resStates = [initial, ...resStates];

  return resStates.map((st) => materialize(helm, spans, st));
}
