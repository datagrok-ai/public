/**
 * Profile-based global alignment for antibody numbering.
 *
 * Aligns an input antibody sequence against a consensus profile
 * using a modified Needleman-Wunsch algorithm. The profile defines
 * which positions exist in the numbering scheme and what amino acids
 * are expected at each position.
 *
 * Gap positions in the profile (empty amino acid lists) represent
 * CDR/variable regions where the input sequence can have insertions.
 */

import type { ChainType, AlignmentResult, Scheme } from './types';
import type { ProfileData } from './data/consensus';
import { getConsensusProfile, getCtermProfile } from './data/consensus';
import { blosum62Score } from './data/blosum62';

// Scoring parameters
const MATCH_CONSENSUS = 6;   // bonus for matching a consensus amino acid
const MATCH_BLOSUM_SCALE = 1; // scale factor for BLOSUM scores on mismatch
const GAP_OPEN = -10;         // penalty for opening a gap
const GAP_EXTEND = -2;        // penalty for extending a gap
const GAP_POSITION_BONUS = 4; // reduce gap penalty at gap positions (CDRs)

/**
 * Score a residue against a profile position.
 * Gap positions (empty consensus) get a flat score since anything goes there.
 */
function scorePosition(aa: string, consensusAAs: string[]): number {
  if (consensusAAs.length === 0) {
    // This is a gap/insertion position in the profile - mild positive score
    return 1;
  }
  // Check if aa is in the consensus set
  if (consensusAAs.includes(aa)) {
    return MATCH_CONSENSUS;
  }
  // Use BLOSUM62 against the best consensus AA
  let best = -10;
  for (const caa of consensusAAs) {
    const s = blosum62Score(aa, caa) * MATCH_BLOSUM_SCALE;
    if (s > best) best = s;
  }
  return best;
}

/** Direction for traceback */
const enum Dir {
  DIAG = 0,
  UP = 1,   // gap in profile (insertion in sequence)
  LEFT = 2, // gap in sequence (deletion from profile)
}

interface AlignmentCell {
  score: number;
  dir: Dir;
}

/**
 * Perform profile-based global alignment.
 * Returns the position code for each residue in the input sequence,
 * or '-' if the residue is outside the numbered region.
 */
function profileAlign(
  seq: string,
  profile: ProfileData,
): { positionCodes: string[]; matchedPositions: number; totalProfilePositions: number } {
  const n = seq.length;       // sequence length
  const m = profile.length;   // profile length

  // DP matrix: (n+1) x (m+1)
  // We use affine gap penalties with separate match/insert/delete states
  const INF = -1e9;

  // M[i][j] = best score aligning seq[0..i-1] to profile[0..j-1] ending in match
  // X[i][j] = ... ending in gap in profile (insertion in seq)
  // Y[i][j] = ... ending in gap in sequence (deletion from profile)
  const M: number[][] = Array.from({ length: n + 1 }, () => new Float64Array(m + 1) as unknown as number[]);
  const X: number[][] = Array.from({ length: n + 1 }, () => new Float64Array(m + 1) as unknown as number[]);
  const Y: number[][] = Array.from({ length: n + 1 }, () => new Float64Array(m + 1) as unknown as number[]);
  const trM: Dir[][] = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1) as unknown as Dir[]);
  const trX: Dir[][] = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1) as unknown as Dir[]);
  const trY: Dir[][] = Array.from({ length: n + 1 }, () => new Uint8Array(m + 1) as unknown as Dir[]);

  // Initialize
  M[0][0] = 0;
  X[0][0] = INF;
  Y[0][0] = INF;

  // Allow free leading gaps in sequence (seq can start after profile begins)
  // and free leading gaps in profile (profile can start after seq begins)
  for (let j = 1; j <= m; j++) {
    // Deletion from profile (skip profile positions at start)
    const isGap = profile[j - 1][1].length === 0;
    M[0][j] = INF;
    X[0][j] = INF;
    // Mild penalty for skipping non-gap positions, free for gap positions
    Y[0][j] = isGap ? Y[0][j - 1] : Y[0][j - 1] + GAP_EXTEND * 0.3;
    if (j === 1) Y[0][j] = isGap ? 0 : GAP_OPEN * 0.3;
    trY[0][j] = Dir.LEFT;
  }

  for (let i = 1; i <= n; i++) {
    // Insertion in sequence (residues before profile starts)
    M[i][0] = INF;
    Y[i][0] = INF;
    X[i][0] = X[i - 1][0] + GAP_EXTEND * 0.3;
    if (i === 1) X[i][0] = GAP_OPEN * 0.3;
    trX[i][0] = Dir.UP;
  }

  // Fill DP matrix
  for (let i = 1; i <= n; i++) {
    const aa = seq[i - 1];
    for (let j = 1; j <= m; j++) {
      const [, consensusAAs] = profile[j - 1];
      const isGapPos = consensusAAs.length === 0;

      // Match/mismatch score
      const matchScore = scorePosition(aa, consensusAAs);

      // M[i][j]: match state - came from any state diagonally
      const mFromM = M[i - 1][j - 1] + matchScore;
      const mFromX = X[i - 1][j - 1] + matchScore;
      const mFromY = Y[i - 1][j - 1] + matchScore;
      if (mFromM >= mFromX && mFromM >= mFromY) {
        M[i][j] = mFromM; trM[i][j] = Dir.DIAG;
      } else if (mFromX >= mFromY) {
        M[i][j] = mFromX; trM[i][j] = Dir.UP;
      } else {
        M[i][j] = mFromY; trM[i][j] = Dir.LEFT;
      }

      // X[i][j]: gap in profile (insertion in sequence)
      // Reduced penalty at positions that are naturally gap positions
      const gapOpenAdj = isGapPos ? GAP_OPEN + GAP_POSITION_BONUS : GAP_OPEN;
      const gapExtAdj = isGapPos ? GAP_EXTEND + GAP_POSITION_BONUS * 0.5 : GAP_EXTEND;

      const xFromM = M[i - 1][j] + gapOpenAdj;
      const xFromX = X[i - 1][j] + gapExtAdj;
      if (xFromM >= xFromX) {
        X[i][j] = xFromM; trX[i][j] = Dir.DIAG;
      } else {
        X[i][j] = xFromX; trX[i][j] = Dir.UP;
      }

      // Y[i][j]: gap in sequence (deletion from profile)
      // More permissive at gap positions (CDR gaps are expected)
      const delOpen = isGapPos ? 0 : GAP_OPEN;
      const delExt = isGapPos ? 0 : GAP_EXTEND;

      const yFromM = M[i][j - 1] + delOpen;
      const yFromY = Y[i][j - 1] + delExt;
      if (yFromM >= yFromY) {
        Y[i][j] = yFromM; trY[i][j] = Dir.DIAG;
      } else {
        Y[i][j] = yFromY; trY[i][j] = Dir.LEFT;
      }
    }
  }

  // Find best end score.
  // We prefer endpoints that consume the full profile (j = m) to ensure
  // all scheme positions are assigned. We allow trailing sequence residues
  // (i < n) for sequences extending beyond the variable region.
  // We also allow j < m for truncated sequences missing C-terminal residues,
  // but apply a penalty for unused non-gap profile positions.
  let bestScore = INF;
  let bestI = n, bestJ = m;
  let bestState: 'M' | 'X' | 'Y' = 'M';

  // Primary: full profile consumed, any amount of sequence consumed
  for (let i = 0; i <= n; i++) {
    for (const [state, mat] of [['M', M], ['X', X], ['Y', Y]] as const) {
      if (mat[i][m] > bestScore) {
        bestScore = mat[i][m]; bestState = state; bestI = i; bestJ = m;
      }
    }
  }

  // Secondary: full sequence consumed, partial profile (for truncated sequences).
  // Apply a small penalty per skipped non-gap profile position.
  for (let j = 0; j < m; j++) {
    // Count skipped non-gap positions at the end
    let skippedNonGap = 0;
    for (let k = j; k < m; k++) {
      if (profile[k][1].length > 0) skippedNonGap++;
    }
    const penalty = skippedNonGap * 3; // mild penalty per skipped position
    for (const [state, mat] of [['M', M], ['X', X], ['Y', Y]] as const) {
      if (mat[n][j] - penalty > bestScore) {
        bestScore = mat[n][j] - penalty; bestState = state; bestI = n; bestJ = j;
      }
    }
  }

  // Traceback
  const alignment: Array<[seqIdx: number, profIdx: number]> = [];
  let ci = bestI, cj = bestJ;
  let curState = bestState;

  // Trailing unaligned residues in sequence
  for (let i = n; i > bestI; i--) {
    alignment.push([i - 1, -1]);
  }

  while (ci > 0 || cj > 0) {
    if (curState === 'M') {
      if (ci === 0 && cj === 0) break;
      if (ci === 0 || cj === 0) {
        // Edge case
        if (ci > 0) { alignment.push([ci - 1, -1]); ci--; curState = 'X'; }
        else { cj--; curState = 'Y'; }
        continue;
      }
      alignment.push([ci - 1, cj - 1]);
      const tr = trM[ci][cj];
      ci--; cj--;
      if (tr === Dir.DIAG) curState = 'M';
      else if (tr === Dir.UP) curState = 'X';
      else curState = 'Y';
    } else if (curState === 'X') {
      if (ci === 0) break;
      alignment.push([ci - 1, -1]); // seq residue not aligned to profile
      const tr = trX[ci][cj];
      ci--;
      if (tr === Dir.DIAG) curState = 'M';
      else curState = 'X';
    } else { // Y
      if (cj === 0) break;
      // Skip profile position (deletion)
      const tr = trY[ci][cj];
      cj--;
      if (tr === Dir.DIAG) curState = 'M';
      else curState = 'Y';
    }
  }

  // Leading unaligned residues
  while (ci > 0) {
    alignment.push([ci - 1, -1]);
    ci--;
  }

  alignment.reverse();

  // Build position codes from alignment
  const positionCodes: string[] = new Array(n).fill('-');
  let matchedPositions = 0;
  const totalProfilePositions = profile.filter(p => p[1].length > 0).length;

  // Track which profile positions were used (for insertion labeling)
  const usedProfilePositions = new Set<number>();
  // First pass: assign direct matches
  for (const [seqIdx, profIdx] of alignment) {
    if (profIdx >= 0 && seqIdx >= 0) {
      const [posNum] = profile[profIdx];
      positionCodes[seqIdx] = String(posNum);
      usedProfilePositions.add(profIdx);
      if (profile[profIdx][1].length > 0) {
        matchedPositions++;
      }
    }
  }

  // Determine the alignment span: only create insertions between the first
  // and last directly-matched profile positions. Residues outside this span
  // (leader peptide, constant region, etc.) should remain as '-'.
  let firstMatchedSeqIdx = -1;
  let lastMatchedSeqIdx = -1;
  for (const [seqIdx, profIdx] of alignment) {
    if (profIdx >= 0 && seqIdx >= 0) {
      if (firstMatchedSeqIdx === -1) firstMatchedSeqIdx = seqIdx;
      lastMatchedSeqIdx = seqIdx;
    }
  }

  // Second pass: handle insertions (residues within the alignment span
  // that weren't directly matched to a profile position)
  // Insertion codes use letters: A, B, C, ...
  const insertionCounters = new Map<string, number>();
  for (let si = 0; si < n; si++) {
    if (positionCodes[si] !== '-') continue;
    // Only create insertions within the alignment span
    if (si < firstMatchedSeqIdx || si > lastMatchedSeqIdx) continue;
    // Find the nearest assigned position before this one
    let prevPos = '';
    for (let k = si - 1; k >= 0; k--) {
      if (positionCodes[k] !== '-') {
        prevPos = positionCodes[k];
        break;
      }
    }
    if (prevPos) {
      // Extract base position number (strip any existing insertion letter)
      const basePos = prevPos.replace(/[A-Z]$/, '');
      const count = (insertionCounters.get(basePos) ?? 0) + 1;
      insertionCounters.set(basePos, count);
      const insertionLetter = String.fromCharCode(64 + count); // A=1, B=2, ...
      positionCodes[si] = basePos + insertionLetter;
    }
  }

  return { positionCodes, matchedPositions, totalProfilePositions };
}

/**
 * Compute percent identity between aligned sequence and consensus.
 */
function computeIdentity(
  seq: string,
  positionCodes: string[],
  profile: ProfileData,
): number {
  const profileMap = new Map<number, string[]>();
  for (const [pos, aas] of profile) {
    if (aas.length > 0) profileMap.set(pos, aas);
  }

  let matches = 0;
  let total = 0;

  for (let i = 0; i < seq.length; i++) {
    const code = positionCodes[i];
    if (code === '-') continue;
    const posNum = parseInt(code, 10);
    if (isNaN(posNum)) continue;
    const consensusAAs = profileMap.get(posNum);
    if (!consensusAAs || consensusAAs.length === 0) continue;
    total++;
    if (consensusAAs.includes(seq[i])) {
      matches++;
    }
  }

  return total > 0 ? matches / total : 0;
}

/**
 * Validate conserved residues for an antibody sequence.
 * Checks for conserved cysteines (positions 23, 104 in IMGT-like schemes)
 * and conserved tryptophan (position 41 in IMGT).
 */
function validateConserved(
  seq: string,
  positionCodes: string[],
  profile: ProfileData,
): string {
  // Build position-to-residue map
  const posToAA = new Map<number, string>();
  for (let i = 0; i < seq.length; i++) {
    const code = positionCodes[i];
    if (code === '-') continue;
    const posNum = parseInt(code, 10);
    if (!isNaN(posNum)) {
      posToAA.set(posNum, seq[i]);
    }
  }

  // Find the conserved cysteine positions by checking profile
  // In most schemes, Cys is at specific positions where profile says ['C']
  const cysPositions: number[] = [];
  for (const [pos, aas] of profile) {
    if (aas.length === 1 && aas[0] === 'C') {
      cysPositions.push(pos);
    }
  }

  for (const pos of cysPositions) {
    const aa = posToAA.get(pos);
    if (aa && aa !== 'C') {
      return `Expected conserved Cys at position ${pos}, found ${aa}`;
    }
  }

  return '';
}

/**
 * Try to find the C-terminal of the variable region.
 * Returns the index in the sequence where the C-terminal motif starts,
 * or -1 if not found.
 */
function findCTerminal(seq: string, chain: ChainType): number {
  const ctermProfile = getCtermProfile(chain);
  const motifLen = ctermProfile.length;

  if (seq.length < motifLen) return -1;

  let bestScore = -Infinity;
  let bestPos = -1;

  // Scan for best match across the sequence (variable region FW4 can be
  // far from the C-terminus in sequences with constant regions)
  const searchStart = Math.max(0, Math.floor(seq.length * 0.3));
  for (let start = searchStart; start <= seq.length - motifLen; start++) {
    let score = 0;
    for (let j = 0; j < motifLen; j++) {
      const aa = seq[start + j];
      const consensusAAs = ctermProfile[j][1];
      if (consensusAAs.includes(aa)) {
        score += 3;
      } else {
        score -= 1;
      }
    }
    if (score > bestScore) {
      bestScore = score;
      bestPos = start;
    }
  }

  // Require a reasonable match
  return bestScore >= motifLen * 1.5 ? bestPos : -1;
}

/**
 * Align a sequence to all chain type profiles and pick the best match.
 */
export function alignSequence(
  seq: string,
  scheme: Scheme,
  chains: ChainType[] = ['H', 'K', 'L'],
): AlignmentResult {
  if (!seq || seq.length < 10) {
    return { numbering: [], percentIdentity: 0, chainType: 'H', error: 'Sequence too short' };
  }

  // Pre-compute CTERM positions for all chain types
  const ctermPositions = new Map<ChainType, number>();
  for (const chain of chains) {
    const idx = findCTerminal(seq, chain);
    if (idx >= 0) ctermPositions.set(chain, idx);
  }

  // Detect scFv: if Heavy and Light CTERM positions are both found
  // and far apart, this is a single-chain variable fragment (VH + VL)
  const ctermH = ctermPositions.get('H') ?? -1;
  const ctermK = ctermPositions.get('K') ?? -1;
  const ctermL = ctermPositions.get('L') ?? -1;
  const ctermLight = Math.max(ctermK, ctermL);
  const lightChain: ChainType = ctermK >= ctermL ? 'K' : 'L';

  const scfvDomainStart = new Map<ChainType, number>();
  if (ctermH >= 0 && ctermLight >= 0 && Math.abs(ctermH - ctermLight) > 80) {
    if (ctermLight < ctermH) {
      // VL-linker-VH: Light domain first, Heavy domain second.
      // Trim Heavy to skip past VL domain so it aligns cleanly to VH.
      const ctLen = getCtermProfile(lightChain).length;
      scfvDomainStart.set('H', ctermLight + ctLen);
    }
    // For VH-VL: Heavy is already at position 0 and aligns well naturally.
    // Don't trim Light - leaving it untrimmed means VH residues confuse
    // the Light alignment, keeping Heavy as the best match.
  }

  const isScfv = ctermH >= 0 && ctermLight >= 0 && Math.abs(ctermH - ctermLight) > 80;

  let bestResult: AlignmentResult | null = null;
  let bestIdentity = -1;
  let heavyResult: AlignmentResult | null = null;
  let heavyIdentity = -1;

  for (const chain of chains) {
    const profile = getConsensusProfile(scheme, chain);

    let trimmedSeq = seq;
    let seqOffset = 0;
    const ctermIdx = ctermPositions.get(chain) ?? -1;
    const ctermLen = getCtermProfile(chain).length;

    // C-terminal trimming
    if (ctermIdx >= 0) {
      const ctermEnd = Math.min(seq.length, ctermIdx + ctermLen + 5);
      trimmedSeq = seq.substring(0, ctermEnd);
    }

    // N-terminal trimming
    const domainStart = scfvDomainStart.get(chain);
    if (domainStart !== undefined && domainStart > 0 && domainStart < trimmedSeq.length - 80) {
      // scFv: skip past the other domain
      seqOffset = domainStart;
      trimmedSeq = trimmedSeq.substring(domainStart);
    } else if (trimmedSeq.length > 180) {
      // Long sequence (with constant region): try a few N-terminal start positions
      const starts = [0];
      const firstAAs = profile.slice(0, 5).map(p => p[1]).flat();
      for (let i = 1; i < Math.min(30, trimmedSeq.length - 80); i++) {
        if (firstAAs.includes(trimmedSeq[i])) {
          starts.push(i);
          break;
        }
      }

      let bestLocalIdentity = -1;
      let bestOffset = 0;

      for (const start of starts) {
        const subSeq = trimmedSeq.substring(start);
        const result = profileAlign(subSeq, profile);
        const identity = computeIdentity(subSeq, result.positionCodes, profile);
        if (identity > bestLocalIdentity) {
          bestLocalIdentity = identity;
          bestOffset = start;
        }
      }

      if (bestOffset > 0) {
        seqOffset = bestOffset;
        trimmedSeq = trimmedSeq.substring(seqOffset);
      }
    }

    const result = profileAlign(trimmedSeq, profile);
    const identity = computeIdentity(trimmedSeq, result.positionCodes, profile);

    // Build full numbering array
    const fullNumbering: string[] = new Array(seq.length).fill('-');
    for (let i = 0; i < trimmedSeq.length; i++) {
      fullNumbering[seqOffset + i] = result.positionCodes[i];
    }

    // KABAT CDR-H1 insertion placement: insertions should go after position 35
    // (not left-aligned after 31). Redistribute if needed.
    if (scheme === 'kabat' && chain === 'H') {
      const cdr1Indices: number[] = [];
      for (let i = 0; i < fullNumbering.length; i++) {
        if (fullNumbering[i] === '-') continue;
        const baseNum = parseInt(fullNumbering[i], 10);
        if (!isNaN(baseNum) && baseNum >= 31 && baseNum <= 35) {
          cdr1Indices.push(i);
        }
      }
      if (cdr1Indices.length > 5) {
        // First 5 get base positions 31-35, rest get 35A, 35B, etc.
        for (let j = 0; j < cdr1Indices.length; j++) {
          if (j < 5) {
            fullNumbering[cdr1Indices[j]] = String(31 + j);
          } else {
            fullNumbering[cdr1Indices[j]] = '35' + String.fromCharCode(64 + j - 4);
          }
        }
      }
    }

    const validationError = validateConserved(trimmedSeq, result.positionCodes, profile);

    const chainResult: AlignmentResult = {
      numbering: fullNumbering,
      percentIdentity: identity,
      chainType: chain,
      error: identity < 0.3 ? 'Low sequence identity; may not be an antibody variable region' :
             validationError || '',
    };

    // Save Heavy result for scFv preference
    if (chain === 'H') {
      heavyResult = chainResult;
      heavyIdentity = identity;
    }

    if (identity > bestIdentity) {
      bestIdentity = identity;
      bestResult = chainResult;
    }
  }

  // For scFv sequences, prefer Heavy chain when identities are close.
  // In scFv constructs, both VH and VL domains match their respective
  // profiles well. AntPack (HMMER-based) naturally picks VH; we replicate
  // this by preferring Heavy when its identity is within 0.05 of the best.
  if (isScfv && bestResult && bestResult.chainType !== 'H' &&
      heavyResult && heavyIdentity > 0.9 && bestIdentity - heavyIdentity < 0.05) {
    bestResult = heavyResult;
  }

  return bestResult!;
}
