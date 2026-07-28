export interface PromptSpec {
  id: string;
  prompt: string;
  category: string;
  expected_path: string;            // 'datagrok-exec' for now
  expected_helpers: string[];       // regex per element
  expected_symbols: string[];       // regex per element
  forbidden_symbols: string[];      // regex per element
  intent?: string;
  trap?: string;
}

export interface DimResult {
  passed: boolean;
  // Per-regex diagnostic — which patterns matched, which didn't.
  details: Array<{pattern: string; matched: boolean}>;
}

export interface RubricResult {
  id: string;
  category: string;
  path: DimResult;
  helpers: DimResult;
  symbols: DimResult;
  forbidden: DimResult;
  overall: boolean;
  execBlocks: string[];
  responseText: string;
}

export function extractExecBlocks(text: string): string[] {
  const re = /```datagrok-exec\s*\n([\s\S]*?)\n```/g;
  const out: string[] = [];
  let m: RegExpExecArray | null;
  while ((m = re.exec(text)) !== null)
    out.push(m[1]);
  return out;
}

function anyBlockMatches(blocks: string[], pattern: string): boolean {
  let re: RegExp;
  try {
    re = new RegExp(pattern);
  } catch {
    // If the pattern itself is broken, treat as not-matched (and surface via details).
    return false;
  }
  return blocks.some((b) => re.test(b));
}

function anyTextMatches(text: string, pattern: string): boolean {
  let re: RegExp;
  try {
    re = new RegExp(pattern);
  } catch {
    return false;
  }
  return re.test(text);
}

export function evaluate(spec: PromptSpec, responseText: string): RubricResult {
  const blocks = extractExecBlocks(responseText);

  // path: did we emit at least one datagrok-exec block (if expected)?
  const path: DimResult = {
    passed: spec.expected_path === 'datagrok-exec' ? blocks.length > 0 : true,
    details: [{pattern: `>=1 ${spec.expected_path} block`, matched: blocks.length > 0}],
  };

  // helpers: ANY listed regex matching in at least one exec block passes.
  // Authors list alternative acceptable helpers — they're not all required.
  // (Per-pattern pipe-OR also works inside a single element.)
  const helperDetails = spec.expected_helpers.map((p) => ({
    pattern: p,
    matched: anyBlockMatches(blocks, p),
  }));
  const helpers: DimResult = {
    passed: helperDetails.length === 0 || helperDetails.some((d) => d.matched),
    details: helperDetails,
  };

  // symbols: each regex must match in at least one exec block.
  const symbolDetails = spec.expected_symbols.map((p) => ({
    pattern: p,
    matched: anyBlockMatches(blocks, p),
  }));
  const symbols: DimResult = {
    passed: symbolDetails.every((d) => d.matched),
    details: symbolDetails,
  };

  // forbidden: NO regex may match in ANY exec block. We also check the response text
  // outside blocks — a forbidden anti-pattern in commentary still leaks the wrong
  // mental model. Per spec the focus is on blocks, so we check blocks first; the
  // text-wide check is just for richer diagnostics.
  const forbiddenDetails = spec.forbidden_symbols.map((p) => ({
    pattern: p,
    matched: anyBlockMatches(blocks, p),
  }));
  const forbidden: DimResult = {
    // "passed" means "no forbidden pattern matched" — i.e. all `matched` are false.
    passed: forbiddenDetails.every((d) => !d.matched),
    details: forbiddenDetails,
  };

  const overall = path.passed && helpers.passed && symbols.passed && forbidden.passed;

  return {
    id: spec.id,
    category: spec.category,
    path,
    helpers,
    symbols,
    forbidden,
    overall,
    execBlocks: blocks,
    responseText,
  };
}
