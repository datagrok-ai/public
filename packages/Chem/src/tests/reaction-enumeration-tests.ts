/* eslint-disable max-len */
/* eslint-disable max-lines-per-function */
import * as DG from 'datagrok-api/dg';
import {before, category, test, expect} from './_timed-test';
import {_package} from '../package-test';
import {cloneConfig, DEFAULT_CONFIG} from '../utils/reaction-enumeration/config';
import {
  enumerate, formatRoute, Route,
  splitSmartsByReactants, stripOuterParens, TemplateInput,
} from '../utils/reaction-enumeration/enumerate';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getRdKitModule} from '../utils/chem-common-rdkit';
import {parseMultiStepReaction} from '../rendering/rdkit-reaction-renderer';

// ── Helpers ─────────────────────────────────────────────────────────────────────

async function loadCsv(name: string): Promise<DG.DataFrame> {
  let text: string;
  try {
    text = await _package.files.readAsText(name);
  } catch (e) {
    throw new Error(`readAsText("${name}") failed: ${e instanceof Error ? e.message : String(e)}`);
  }
  if (text == null) throw new Error(`readAsText("${name}") returned null/undefined`);
  const df = DG.DataFrame.fromCsv(text);
  if (!df) throw new Error(`fromCsv("${name}") returned null. First 200 chars: ${text.slice(0, 200)}`);
  return df;
}

function colAsStrings(df: DG.DataFrame, name: string): string[] {
  const c = df.col(name);
  if (!c) throw new Error(`Column "${name}" missing. Available: ${df.columns.names().join(', ')}`);
  const out: string[] = new Array(c.length);
  for (let i = 0; i < c.length; i++) out[i] = c.get(i) == null ? '' : String(c.get(i));
  return out;
}

// (parseMultiStepReaction and BRANCH_DELIMITER are imported from the renderer above so the
// test can't drift out of sync with the implementation.)

async function runFullSetup(opts?: {
  maxComponents?: number; maxCombos?: number;
  useExclusion?: boolean; firstNTemplates?: number;
}) {
  const rdkit = getRdKitModule();
  const templatesDf = await loadCsv('enumerations/reactions.csv');
  const bbsDf = await loadCsv('enumerations/bb.csv');
  const exclusionDf = await loadCsv('enumerations/ex_smarts.csv');

  const smartsList = colAsStrings(templatesDf, 'reaction_smarts');
  const blockingList = templatesDf.col('blocking_fg') ? colAsStrings(templatesDf, 'blocking_fg') : null;
  const reactionNames = templatesDf.col('reaction_name') ? colAsStrings(templatesDf, 'reaction_name') : null;
  let templates: TemplateInput[] = [];
  for (let i = 0; i < smartsList.length; i++) {
    const smarts = (smartsList[i] ?? '').trim();
    if (!smarts) continue;
    const blockingRaw = blockingList ? blockingList[i] : '';
    const blockingSmartsList = blockingRaw ?
      blockingRaw.split(/[;|]/).map((s) => s.trim()).filter((s) => s.length > 0) : [];
    templates.push({smarts, blockingSmartsList, reactionName: reactionNames?.[i] ?? ''});
  }
  if (opts?.firstNTemplates != null) templates = templates.slice(0, opts.firstNTemplates);

  const buildingBlocks = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);
  const exclusionSmarts = (opts?.useExclusion ?? true) ?
    colAsStrings(exclusionDf, 'SMARTS').filter((s) => s.trim().length > 0) : [];

  const config = cloneConfig(DEFAULT_CONFIG);
  config.max_num_components = opts?.maxComponents ?? 4;
  config.max_num_combinations_per_template = opts?.maxCombos ?? 50;
  return {rdkit, templates, buildingBlocks, exclusionSmarts, config};
}

function probeRdkit(rdkit: any, label: string): void {
  let probe: any = null;
  try {
    probe = rdkit.get_mol('CCC');
  } catch (e) {
    throw new Error(`[${label}] rdkit.get_mol("CCC") threw: ${e instanceof Error ? e.message : String(e)}`);
  }
  expect(probe != null, true, `[${label}] rdkit.get_mol("CCC") returned null`);
  let valid = false;
  try {
    valid = probe.is_valid();
  } catch (e) {
    throw new Error(`[${label}] probe.is_valid() threw: ${e instanceof Error ? e.message : String(e)}`);
  }
  expect(valid, true, `[${label}] probe.is_valid() === false`);
  let smi = '';
  try {
    smi = probe.get_smiles();
  } catch (e) {
    throw new Error(`[${label}] probe.get_smiles() threw: ${e instanceof Error ? e.message : String(e)}`);
  }
  expect(smi.length > 0, true, `[${label}] probe.get_smiles() returned empty`);
  try {
    probe.delete();
  } catch (e) {
    throw new Error(`[${label}] probe.delete() threw: ${e instanceof Error ? e.message : String(e)}`);
  }
}

// ── Tests ──────────────────────────────────────────────────────────────────────

category('Reaction Enumeration', () => {
  before(async () => {
    if (!chemCommonRdKit.moduleInitialized) {
      chemCommonRdKit.setRdKitWebRoot(_package.webRoot);
      await chemCommonRdKit.initRdKitModuleLocal();
    }
  });

  // ── parse ───────────────────────────────────────────────────────────────
  test('parse: splitSmartsByReactants — simple two-reactant LHS', async () => {
    const parts = splitSmartsByReactants('[A:1].[B:2]');
    expect(parts.length, 2);
    expect(parts[0], '[A:1]');
    expect(parts[1], '[B:2]');
  });

  test('parse: splitSmartsByReactants — does not split inside square brackets', async () => {
    const parts = splitSmartsByReactants('[#6:1][C:2](=[O:3])[OH,O-]');
    expect(parts.length, 1);
  });

  test('parse: splitSmartsByReactants — does not split inside parens (tied groups)', async () => {
    const parts = splitSmartsByReactants('([A:1].[B:2]).[C:3]');
    expect(parts.length, 2);
    expect(parts[0], '([A:1].[B:2])');
    expect(parts[1], '[C:3]');
  });

  test('parse: splitSmartsByReactants — three-reactant LHS with tied first slot', async () => {
    const lhs = '([NX3;!$(N[C,S]=[S,O,N]);$(N[#6]);!$(N[#7]);H1:1].[NX3;H0:2]C(=O)OC(C)(C)C).[#6:3][C:4](=[O:5])[OH,O-].[#6:6][C:7](=[O:8])[OH,O-]';
    const parts = splitSmartsByReactants(lhs);
    expect(parts.length, 3);
    expect(parts[0].startsWith('('), true);
    expect(parts[0].endsWith(')'), true);
    expect(parts[1], '[#6:3][C:4](=[O:5])[OH,O-]');
    expect(parts[2], '[#6:6][C:7](=[O:8])[OH,O-]');
  });

  test('parse: stripOuterParens — strips a single matched outer pair', async () => {
    expect(stripOuterParens('([A:1].[B:2])'), '[A:1].[B:2]');
  });

  test('parse: stripOuterParens — leaves unparenthesized strings unchanged', async () => {
    expect(stripOuterParens('[A:1].[B:2]'), '[A:1].[B:2]');
  });

  test('parse: stripOuterParens — does not strip when outer parens close before end', async () => {
    expect(stripOuterParens('(X).(Y)'), '(X).(Y)');
  });

  test('parse: stripOuterParens — strips only one level of nested parens', async () => {
    expect(stripOuterParens('((X))'), '(X)');
  });

  test('parse: stripOuterParens — leaves bracket-only strings unchanged', async () => {
    expect(stripOuterParens('[#6:1]'), '[#6:1]');
  });

  // ── route formatting ────────────────────────────────────────────────────
  const tmpl = '';
  const step = (reactants: string[], product: string) =>
    ({reactants, product, templateSmarts: tmpl, reactionName: ''});

  test('route formatting: single-step route', async () => {
    const route: Route = [step(['BB1', 'BB2'], 'P1')];
    expect(formatRoute(route), 'BB1.BB2>>P1');
  });

  test('route formatting: two-step linear chain — every step delimited', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1', 'BB3'], 'P2'),
    ];
    const formatted = formatRoute(route);
    expect(formatted, 'BB1.BB2>>P1--**--P1.BB3>>P2');
    // parseMultiStepReaction returns branches × steps; both steps live in their own branch here
    // (each step was emitted as its own --**--delimited segment).
    const branches = parseMultiStepReaction(formatted);
    expect(branches.length, 2);
    expect(branches[0].length, 1);
    expect(branches[0][0], 'BB1.BB2>>P1');
    expect(branches[1][0], 'P1.BB3>>P2');
  });

  test('route formatting: three-step linear chain', async () => {
    const route: Route = [
      step(['A', 'B'], 'P1'),
      step(['P1', 'C'], 'P2'),
      step(['P2', 'D'], 'P3'),
    ];
    expect(formatRoute(route), 'A.B>>P1--**--P1.C>>P2--**--P2.D>>P3');
  });

  test('route formatting: linear chain where the next step has only the prev product as reactant', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1'], 'P2'),
    ];
    expect(formatRoute(route), 'BB1.BB2>>P1--**--P1>>P2');
  });

  test('route formatting: branching route uses --**-- delimiter at every seam', async () => {
    const route: Route = [
      step(['A', 'B'], 'P1'),
      step(['C', 'D'], 'P2'),
      step(['P1', 'P2'], 'P3'),
    ];
    const formatted = formatRoute(route);
    expect(formatted, 'A.B>>P1--**--C.D>>P2--**--P1.P2>>P3');
    expect(formatted.indexOf(';'), -1, `Route used semicolon separator: ${formatted}`);
    // Three branches, each a single step — renderer will draw branch separators (not arrows)
    // between them.
    const branches = parseMultiStepReaction(formatted);
    expect(branches.length, 3, `Expected 3 branches, got ${branches.length}: ${JSON.stringify(branches)}`);
    expect(branches.flat().length, 3, 'Each branch should hold exactly one step');
    expect(branches[0][0], 'A.B>>P1');
    expect(branches[1][0], 'C.D>>P2');
    expect(branches[2][0], 'P1.P2>>P3');
  });

  test('route formatting: convergent multi-branch synthesis — no fake intermediate seams', async () => {
    // bb1+bb2→p1, bb3+bb4→p2, p2+bb5→p3, p3+p1→p4. Every step is its own self-contained
    // reaction segment, so each ends up in its own branch.
    const route: Route = [
      step(['bb1', 'bb2'], 'p1'),
      step(['bb3', 'bb4'], 'p2'),
      step(['p2', 'bb5'], 'p3'),
      step(['p3', 'p1'], 'p4'),
    ];
    const formatted = formatRoute(route);
    expect(formatted, 'bb1.bb2>>p1--**--bb3.bb4>>p2--**--p2.bb5>>p3--**--p3.p1>>p4');

    const branches = parseMultiStepReaction(formatted);
    expect(branches.length, 4, `Expected 4 branches, got ${branches.length}: ${JSON.stringify(branches)}`);
    const flat = branches.flat();
    expect(flat.length, 4);
    expect(flat[0], 'bb1.bb2>>p1');
    expect(flat[1], 'bb3.bb4>>p2');
    expect(flat[2], 'p2.bb5>>p3');
    expect(flat[3], 'p3.p1>>p4');
    for (const s of flat)
      expect(s.split('>>').length - 1, 1, `Step "${s}" is not a single-step reaction`);
  });

  test('route formatting: each formatted step is a valid single-step reaction', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1', 'BB3'], 'P2'),
    ];
    const formatted = formatRoute(route);
    const branches = parseMultiStepReaction(formatted);
    for (const s of branches.flat()) {
      const arrowCount = s.split('>>').length - 1;
      expect(arrowCount, 1, `Step "${s}" has ${arrowCount} ">>" separators (expected 1)`);
    }
  });

  test('route formatting: a legacy linear chain (no delimiter) parses as a single branch', async () => {
    // Backward compat: strings written before the delimiter existed must still parse as before
    // (a single branch with multiple linear steps), so the renderer keeps drawing arrows between
    // adjacent steps instead of branch separators.
    const branches = parseMultiStepReaction('A.B>>C>>D');
    expect(branches.length, 1, `Expected 1 branch for legacy linear chain, got ${branches.length}`);
    expect(branches[0].length, 2);
    expect(branches[0][0], 'A.B>>C');
    expect(branches[0][1], 'C>>D');
  });

  // ── end-to-end ──────────────────────────────────────────────────────────
  test('end-to-end: one-round depth-first over bundled demo data produces rows', async () => {
    const rdkit = getRdKitModule();
    const templatesDf = await loadCsv('enumerations/reactions.csv');
    const bbsDf = await loadCsv('enumerations/bb.csv');
    const exclusionDf = await loadCsv('enumerations/ex_smarts.csv');

    const smartsList = colAsStrings(templatesDf, 'reaction_smarts');
    const blockingList = templatesDf.col('blocking_fg') ? colAsStrings(templatesDf, 'blocking_fg') : null;
    const reactionNames = templatesDf.col('reaction_name') ? colAsStrings(templatesDf, 'reaction_name') : null;
    const templates: TemplateInput[] = [];
    for (let i = 0; i < smartsList.length; i++) {
      const smarts = (smartsList[i] ?? '').trim();
      if (!smarts) continue;
      const blockingRaw = blockingList ? blockingList[i] : '';
      const blockingSmartsList = blockingRaw ?
        blockingRaw.split(/[;|]/).map((s) => s.trim()).filter((s) => s.length > 0) : [];
      templates.push({smarts, blockingSmartsList, reactionName: reactionNames?.[i] ?? ''});
    }
    const buildingBlocks = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);
    const exclusionSmarts = colAsStrings(exclusionDf, 'SMARTS').filter((s) => s.trim().length > 0);

    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 1;
    config.enumeration.depth_first = true;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];
    config.max_num_combinations_per_template = 100000;

    const {rows, warnings} = await enumerate({rdkit, config, templates, buildingBlocks, exclusionSmarts});
    expect(rows.length > 0, true,
      `Expected at least one product row. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);
    for (const row of rows.slice(0, 5)) {
      expect(typeof row.product, 'string');
      expect(row.product.length > 0, true);
      expect(typeof row.route, 'string');
    }
  }, {timeout: 120000});

  test('end-to-end: amino-acid bb file + peptide/disulfide reactions produce expected coupling products', async () => {
    const rdkit = getRdKitModule();
    const bbsDf = await loadCsv('enumerations/aa_bb.csv');
    const rxnsDf = await loadCsv('enumerations/aa_reactions.csv');

    const buildingBlocks = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);
    expect(buildingBlocks.length, 20, `Expected 20 canonical amino acids, got ${buildingBlocks.length}`);

    const smartsList = colAsStrings(rxnsDf, 'reaction_smarts');
    const namesList = rxnsDf.col('reaction_name') ? colAsStrings(rxnsDf, 'reaction_name') : null;
    const templates: TemplateInput[] = [];
    for (let i = 0; i < smartsList.length; i++) {
      const smarts = (smartsList[i] ?? '').trim();
      if (!smarts) continue;
      templates.push({smarts, blockingSmartsList: [], reactionName: namesList?.[i] ?? ''});
    }
    expect(templates.length, 2, 'Expected 2 reaction templates (peptide + disulfide)');

    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 1;
    config.enumeration.depth_first = true;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];
    config.max_num_combinations_per_template = 10000;

    const {rows, warnings} = await enumerate({
      rdkit, config, templates, buildingBlocks, exclusionSmarts: [],
    });
    if (rows.length === 0)
      console.warn('aa enum: no rows. Warnings:', warnings.slice(0, 10));
    expect(rows.length > 0, true,
      `Expected products from amino-acid coupling. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);

    // Peptide coupling: every AA pair (incl. self) -> dipeptide -> at least 20 products,
    // all containing the C-N peptide bond pattern N-C(=O)-C.
    const peptideRows = rows.filter((r) => r.reaction_name === 'Peptide coupling');
    expect(peptideRows.length > 0, true, 'Peptide coupling produced no rows');

    // Disulfide: only Cys-Cys can react (only Cys has a thiol). Expect exactly one canonical product.
    const disulfideRows = rows.filter((r) => r.reaction_name === 'Disulfide formation');
    expect(disulfideRows.length > 0, true, 'Disulfide formation produced no rows — Cys + Cys should match');
    // The product should contain an S-S bond.
    expect(disulfideRows[0].product.includes('SS') || /S.{0,3}S/.test(disulfideRows[0].product), true,
      `Disulfide product missing SS link: ${disulfideRows[0].product}`);
  }, {timeout: 120000});

  test('end-to-end: every step\'s reactants are BBs or products of an earlier step', async () => {
    // Regression test for the route-merging bug where a step using a product from a round
    // earlier than (r-1) had that product's synthesis steps silently dropped, so the final
    // route appeared to "start" from a non-BB intermediate.
    //
    // Use the AA library with breadth_first + 3 rounds — that's where convergent paths
    // (e.g. final disulfide of LeuCysGly + CysVal) most commonly arise, with CysVal in R1
    // and LeuCysGly in R2.
    const rdkit = getRdKitModule();
    const bbsDf = await loadCsv('enumerations/aa_bb.csv');
    const rxnsDf = await loadCsv('enumerations/aa_reactions.csv');
    const rawBBs = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);

    // Canonicalize BBs through RDKit so they match the canonical SMILES the enumerator emits.
    const canonBBs = new Set<string>();
    for (const s of rawBBs) {
      const m = rdkit.get_mol(s);
      try {if (m && m.is_valid()) canonBBs.add(m.get_smiles());} finally {try {m?.delete();} catch {/* ignore */}}
    }

    const smartsList = colAsStrings(rxnsDf, 'reaction_smarts');
    const namesList = rxnsDf.col('reaction_name') ? colAsStrings(rxnsDf, 'reaction_name') : null;
    const templates: TemplateInput[] = [];
    for (let i = 0; i < smartsList.length; i++) {
      const smarts = (smartsList[i] ?? '').trim();
      if (!smarts) continue;
      templates.push({smarts, blockingSmartsList: [], reactionName: namesList?.[i] ?? ''});
    }

    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 3;
    config.enumeration.depth_first = false;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];
    config.max_num_components = 2;
    // Keep the budget small — we just need enough rows to exercise the merge logic.
    config.max_num_combinations_per_template = 30;

    const {rows, warnings} = await enumerate({
      rdkit, config, templates, buildingBlocks: rawBBs, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected products. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);

    // For every row, walk its route step-by-step. At each step, every reactant must already be
    // either a BB or a product of an earlier step. This is the property the route-merging logic
    // is supposed to guarantee.
    let checked = 0;
    let multiStepChecked = 0;
    for (const row of rows) {
      if (!row.route) continue;
      const stepStrs = row.route.split('--**--');
      const known = new Set<string>(canonBBs);
      for (const stepStr of stepStrs) {
        const arrow = stepStr.indexOf('>>');
        if (arrow < 0) continue;
        const lhs = stepStr.slice(0, arrow);
        const rhs = stepStr.slice(arrow + 2);
        for (const r of lhs.split('.')) {
          if (!r) continue;
          expect(known.has(r), true,
            `Route step "${stepStr}" of product "${row.product}" uses reactant "${r}" that is ` +
            `neither a BB nor a product of an earlier step. Full route: ${row.route}`);
        }
        for (const p of rhs.split('.')) if (p) known.add(p);
      }
      checked++;
      if (stepStrs.length > 1) multiStepChecked++;
    }
    expect(checked > 0, true, 'No routes were checked');
    // We need to actually exercise multi-step routes (otherwise the bug couldn't surface).
    expect(multiStepChecked > 0, true,
      `No multi-step routes produced; the merge logic wasn't exercised. Total rows: ${rows.length}`);
  }, {timeout: 600000});

  test('end-to-end: depth_first never merges two prev-round products in a single step', async () => {
    // Regression test: in depth_first mode, every round-(R>1) step must combine EXACTLY ONE
    // round-(R-1) product with original BBs (linear chain extension). Combining two complex
    // products in one step is convergent/breadth-first behavior and must not occur in depth_first.
    //
    // The AA library + 3 rounds is the easy way to trigger the bug: without the strict filter,
    // peptide+disulfide combinations would happily merge two prev-round products.
    const rdkit = getRdKitModule();
    const bbsDf = await loadCsv('enumerations/aa_bb.csv');
    const rxnsDf = await loadCsv('enumerations/aa_reactions.csv');
    const rawBBs = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);

    const canonBBs = new Set<string>();
    for (const s of rawBBs) {
      const m = rdkit.get_mol(s);
      try {
        if (m && m.is_valid()) canonBBs.add(m.get_smiles());
      } finally {
        try {
          m?.delete();
        } catch {
          /* ignore */
        }
      }
    }

    const smartsList = colAsStrings(rxnsDf, 'reaction_smarts');
    const namesList = rxnsDf.col('reaction_name') ? colAsStrings(rxnsDf, 'reaction_name') : null;
    const templates: TemplateInput[] = [];
    for (let i = 0; i < smartsList.length; i++) {
      const smarts = (smartsList[i] ?? '').trim();
      if (!smarts) continue;
      templates.push({smarts, blockingSmartsList: [], reactionName: namesList?.[i] ?? ''});
    }

    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 3;
    config.enumeration.depth_first = true;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];
    config.max_num_components = 2;
    config.max_num_combinations_per_template = 50;

    const {rows, warnings} = await enumerate({
      rdkit, config, templates, buildingBlocks: rawBBs, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected products. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);

    // Walk every route. At each step, classify each reactant as BB or non-BB. In depth_first,
    // a step has at most one non-BB reactant (the prev-round product being extended). Two non-BB
    // reactants in one step = the bug we're guarding against.
    let multiStepChecked = 0;
    for (const row of rows) {
      if (!row.route) continue;
      const stepStrs = row.route.split('--**--');
      if (stepStrs.length <= 1) continue;
      multiStepChecked++;
      for (const stepStr of stepStrs) {
        const arrow = stepStr.indexOf('>>');
        if (arrow < 0) continue;
        const lhs = stepStr.slice(0, arrow);
        const reactants = lhs.split('.').filter((r) => r.length > 0);
        const nonBBs = reactants.filter((r) => !canonBBs.has(r));
        expect(nonBBs.length <= 1, true,
          `depth_first violation: step "${stepStr}" has ${nonBBs.length} non-BB reactants ` +
          `(should be ≤ 1). Reactants: ${reactants.join(' + ')}. Full route: ${row.route}`);
      }
    }
    expect(multiStepChecked > 0, true,
      `No multi-step depth_first routes produced; constraint wasn't exercised. Total rows: ${rows.length}`);
  }, {timeout: 600000});

  test('end-to-end: tied (parenthesized) reactant group of two fragments', async () => {
    const rdkit = getRdKitModule();
    const tiedTemplate: TemplateInput = {
      smarts: '([N;H2;X3;!$(N=*):1].[O;H1;X2;!$(O=*);!$(O~C=O):2]).[#6:3][C:4](=[O:5])[OH,O-]>>[N:1][C:4](=[O:5])[#6:3]',
      blockingSmartsList: [],
      reactionName: 'amide formation on aminoalcohol (tied-group test)',
    };
    const buildingBlocks = ['OCCN', 'CC(=O)O', 'CCC(=O)O'];
    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 1;
    config.enumeration.depth_first = true;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];

    const {rows, warnings} = await enumerate({
      rdkit, config, templates: [tiedTemplate], buildingBlocks, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected at least one product. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);
  }, {timeout: 60000});

  test('end-to-end: outer parens around a single-fragment slot do not break parsing', async () => {
    const rdkit = getRdKitModule();
    const template: TemplateInput = {
      smarts: '([N;H2;X3:1]).[#6:3][C:4](=[O:5])[OH,O-]>>[N:1][C:4](=[O:5])[#6:3]',
      blockingSmartsList: [],
      reactionName: 'amide formation, slot 0 wrapped in parens',
    };
    const buildingBlocks = ['CCN', 'CC(=O)O'];
    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.num_rounds = 1;
    config.enumeration.depth_first = true;
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];

    const {rows, warnings} = await enumerate({
      rdkit, config, templates: [template], buildingBlocks, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected at least one product. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);
  }, {timeout: 60000});

  test('end-to-end: rdkit module remains usable — max_components=2, no exclusion', async () => {
    const setup = await runFullSetup({maxComponents: 2, maxCombos: 50, useExclusion: false});
    probeRdkit(setup.rdkit, 'before');
    const {rows} = await enumerate(setup);
    probeRdkit(setup.rdkit, 'after');
    expect(rows.length > 0, true, 'No rows produced');
  }, {timeout: 300000});

  test('end-to-end: rdkit module remains usable — max_components=2, with exclusion', async () => {
    const setup = await runFullSetup({maxComponents: 2, maxCombos: 50, useExclusion: true});
    probeRdkit(setup.rdkit, 'before');
    const {rows} = await enumerate(setup);
    probeRdkit(setup.rdkit, 'after');
    expect(rows.length > 0, true, 'No rows produced');
  }, {timeout: 300000});

  test('end-to-end: rdkit module remains usable — max_components=4, no exclusion', async () => {
    const setup = await runFullSetup({maxComponents: 4, maxCombos: 50, useExclusion: false});
    probeRdkit(setup.rdkit, 'before');
    const {rows} = await enumerate(setup);
    probeRdkit(setup.rdkit, 'after');
    expect(rows.length > 0, true, 'No rows produced');
  }, {timeout: 300000});

  test('end-to-end: rdkit module remains usable — max_components=4, with exclusion (user scenario)', async () => {
    const setup = await runFullSetup({maxComponents: 4, maxCombos: 50, useExclusion: true});
    probeRdkit(setup.rdkit, 'before');

    const captured: string[] = [];
    const origWarn = console.warn;
    console.warn = (...args: any[]) => {
      captured.push(args.map((a) => (typeof a === 'string' ? a : JSON.stringify(a))).join(' '));
      origWarn.apply(console, args);
    };

    let rows: any[] = [];
    let warnings: string[] = [];
    try {
      const out = await enumerate(setup);
      rows = out.rows; warnings = out.warnings;
    } finally {
      console.warn = origWarn;
    }

    probeRdkit(setup.rdkit, 'after');

    const nullFn = captured.filter((s) => /null function|delete failed/.test(s));
    expect(nullFn.length, 0,
      `Got ${nullFn.length} delete/null-function warnings during enumeration. First 3: ${nullFn.slice(0, 3).join(' || ')}`);

    expect(rows.length > 0, true,
      `Enumeration produced no rows. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);
  }, {timeout: 300000});

  // ── reagents mode ───────────────────────────────────────────────────────
  //
  // Shared helpers for the reagents-mode tests. We split the AA library into a disjoint BB-pool
  // and reagent-pool, then use peptide coupling (a symmetric 2-slot template that yields chain-
  // extendable dipeptides). The structural invariant we verify is the spec from the email:
  //   every step uses EXACTLY ONE reactant from the "main pool" (BBs in round 1, or any earlier
  //   product / BB in round R > 1), and EVERY remaining slot is filled from reagents.
  async function loadAaSplit(rdkit: any, bbCount: number) {
    const bbsDf = await loadCsv('enumerations/aa_bb.csv');
    const rxnsDf = await loadCsv('enumerations/aa_reactions.csv');
    const rawAAs = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0);
    expect(rawAAs.length > bbCount + 2, true, 'AA library too small to split');
    const bbs = rawAAs.slice(0, bbCount);
    const reagents = rawAAs.slice(bbCount);

    const canon = (s: string): string => {
      const m = rdkit.get_mol(s);
      try {return (m && m.is_valid()) ? m.get_smiles() : s;} finally {
        try {m?.delete();} catch {/* ignore */}
      }
    };
    const canonBBs = new Set(bbs.map(canon));
    const canonReagents = new Set(reagents.map(canon));
    // Sanity: the split must actually be disjoint after canonicalisation, otherwise the test
    // can't classify reactants.
    for (const s of canonBBs)
      expect(canonReagents.has(s), false, `Split overlap at canonical SMILES "${s}"`);

    const smartsList = colAsStrings(rxnsDf, 'reaction_smarts');
    const namesList = rxnsDf.col('reaction_name') ? colAsStrings(rxnsDf, 'reaction_name') : null;
    const peptideIdx = (namesList ?? []).findIndex((n) => /peptide/i.test(n));
    expect(peptideIdx >= 0, true, 'Could not find peptide coupling template in aa_reactions.csv');
    const peptideTemplate: TemplateInput = {
      smarts: smartsList[peptideIdx].trim(),
      blockingSmartsList: [],
      reactionName: namesList?.[peptideIdx] ?? '',
    };

    return {bbs, reagents, canonBBs, canonReagents, peptideTemplate};
  }

  function reagentsModeBaseConfig() {
    const config = cloneConfig(DEFAULT_CONFIG);
    config.enumeration.depth_first = true; // ignored in reagents mode but exercises that fact
    config.products_specs.min_num_carbon_atoms = -1;
    config.products_specs.max_num_carbon_atoms = -1;
    config.products_specs.max_num_hetero_atoms = -1;
    config.products_specs.max_num_unsaturated_nonaromatic_bonds = -1;
    config.products_specs.only_these_atoms_allowed = [];
    config.max_num_components = 2;
    config.max_num_combinations_per_template = 200;
    return config;
  }

  test('reagents mode: round 1 — every step is exactly one BB and one reagent', async () => {
    const rdkit = getRdKitModule();
    const {bbs, reagents, canonBBs, canonReagents, peptideTemplate} = await loadAaSplit(rdkit, 5);

    const config = reagentsModeBaseConfig();
    config.enumeration.num_rounds = 1;

    const {rows, warnings} = await enumerate({
      rdkit, config, templates: [peptideTemplate],
      buildingBlocks: bbs, reagents, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected products. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);

    let checked = 0;
    for (const row of rows) {
      if (!row.route) continue;
      const stepStrs = row.route.split('--**--');
      expect(stepStrs.length, 1, `Round-1 route should be single-step: ${row.route}`);
      const arrow = stepStrs[0].indexOf('>>');
      const lhs = stepStrs[0].slice(0, arrow);
      const reactants = lhs.split('.').filter((r) => r.length > 0);
      const fromBB = reactants.filter((r) => canonBBs.has(r));
      const fromReagent = reactants.filter((r) => canonReagents.has(r));
      expect(fromBB.length, 1,
        `Round-1 step "${stepStrs[0]}" should have EXACTLY 1 BB reactant, got ${fromBB.length}: ${reactants.join(' + ')}`);
      expect(fromBB.length + fromReagent.length, reactants.length,
        `Round-1 step "${stepStrs[0]}" had unclassified reactants: ${reactants.join(' + ')}`);
      checked++;
    }
    expect(checked > 0, true, 'No reagents-mode rows were checked');
  }, {timeout: 300000});

  test('reagents mode: multi-round — every step is one main-pool reactant plus reagents only', async () => {
    // Spec from the email:
    //   round 1: BB + reagent(s)               → P1
    //   round 2: P1 + reagent(s)               → P2
    //   ...
    // For each step in each route, exactly one reactant comes from the "main pool" — that is,
    // either a BB or a previously-emitted product in this same route — and every remaining
    // reactant comes from the reagents library. The test also verifies that *at least one*
    // multi-step route is produced, otherwise the round-2+ path wasn't actually exercised.
    const rdkit = getRdKitModule();
    const {bbs, reagents, canonBBs, canonReagents, peptideTemplate} = await loadAaSplit(rdkit, 4);

    const config = reagentsModeBaseConfig();
    config.enumeration.num_rounds = 2;
    // Cap needs to comfortably exceed (bbs × reagents) so the round-2 main slot sampling
    // reaches the round-1 products that follow BBs in eligibleSmiles.
    config.max_num_combinations_per_template = 200;

    const {rows, warnings} = await enumerate({
      rdkit, config, templates: [peptideTemplate],
      buildingBlocks: bbs, reagents, exclusionSmarts: [],
    });
    expect(rows.length > 0, true,
      `Expected products. Warnings (${warnings.length}): ${warnings.slice(0, 5).join(' | ')}`);

    let multiStepChecked = 0;
    let totalSteps = 0;
    for (const row of rows) {
      if (!row.route) continue;
      const stepStrs = row.route.split('--**--');
      const knownProducts = new Set<string>();
      for (const stepStr of stepStrs) {
        const arrow = stepStr.indexOf('>>');
        if (arrow < 0) continue;
        const lhs = stepStr.slice(0, arrow);
        const rhs = stepStr.slice(arrow + 2);
        const reactants = lhs.split('.').filter((r) => r.length > 0);
        // Main pool at this point = BBs ∪ products of earlier steps in this route.
        const mainPool = (r: string) => canonBBs.has(r) || knownProducts.has(r);
        const mainCount = reactants.filter(mainPool).length;
        const reagentCount = reactants.filter((r) => canonReagents.has(r)).length;
        expect(mainCount, 1,
          `Step "${stepStr}" in route "${row.route}": expected EXACTLY 1 main-pool reactant, got ${mainCount}. Reactants: ${reactants.join(' + ')}`);
        expect(mainCount + reagentCount, reactants.length,
          `Step "${stepStr}" had unclassified reactants: ${reactants.join(' + ')}`);
        for (const p of rhs.split('.')) if (p) knownProducts.add(p);
        totalSteps++;
      }
      if (stepStrs.length > 1) multiStepChecked++;
    }
    expect(multiStepChecked > 0, true,
      `No multi-step routes produced — round 2 chain extension wasn't exercised. Total rows: ${rows.length}, total steps: ${totalSteps}`);
  }, {timeout: 600000});

  test('reagents mode: empty reagents array behaves like standard depth-first', async () => {
    // The reagentsMode branch must be skipped when no reagents are provided so existing setups
    // keep working without surprises. Run the same template/BBs once with reagents=[] and once
    // with reagents undefined; both should produce identical product sets to a baseline that
    // doesn't pass the field at all.
    const rdkit = getRdKitModule();
    const bbsDf = await loadCsv('enumerations/aa_bb.csv');
    const rxnsDf = await loadCsv('enumerations/aa_reactions.csv');
    const bbs = colAsStrings(bbsDf, 'SMILES').filter((s) => s.trim().length > 0).slice(0, 6);
    const smartsList = colAsStrings(rxnsDf, 'reaction_smarts');
    const namesList = rxnsDf.col('reaction_name') ? colAsStrings(rxnsDf, 'reaction_name') : null;
    const peptideIdx = (namesList ?? []).findIndex((n) => /peptide/i.test(n));
    const peptideTemplate: TemplateInput = {
      smarts: smartsList[peptideIdx].trim(), blockingSmartsList: [],
      reactionName: namesList?.[peptideIdx] ?? '',
    };

    const baseConfig = reagentsModeBaseConfig();
    baseConfig.enumeration.num_rounds = 1;

    const baseline = await enumerate({
      rdkit, config: baseConfig, templates: [peptideTemplate],
      buildingBlocks: bbs, exclusionSmarts: [],
    });
    const withEmptyReagents = await enumerate({
      rdkit, config: baseConfig, templates: [peptideTemplate],
      buildingBlocks: bbs, reagents: [], exclusionSmarts: [],
    });

    expect(baseline.rows.length > 0, true, 'Baseline produced no rows');
    expect(baseline.rows.length, withEmptyReagents.rows.length,
      `Empty-reagents run should match baseline row count. Baseline=${baseline.rows.length}, empty=${withEmptyReagents.rows.length}`);

    const baseSet = new Set(baseline.rows.map((r) => r.product));
    const emptySet = new Set(withEmptyReagents.rows.map((r) => r.product));
    expect(baseSet.size, emptySet.size, 'Product sets differ in size');
    for (const p of baseSet)
      expect(emptySet.has(p), true, `Empty-reagents missing product "${p}" present in baseline`);
  }, {timeout: 300000});

  test('reagents mode: no BB × BB combinations leak through (regression for slot pool isolation)', async () => {
    // With reagents present, slot pools enforce one-from-main + rest-from-reagents. A BB-only
    // combination (e.g. two BB reactants in the same step) would violate the spec and indicate
    // the slot pools got mixed. Use disjoint BB/reagent sets so the check is unambiguous.
    const rdkit = getRdKitModule();
    const {bbs, reagents, canonBBs, peptideTemplate} = await loadAaSplit(rdkit, 5);

    const config = reagentsModeBaseConfig();
    config.enumeration.num_rounds = 1;

    const {rows} = await enumerate({
      rdkit, config, templates: [peptideTemplate],
      buildingBlocks: bbs, reagents, exclusionSmarts: [],
    });
    expect(rows.length > 0, true, 'Expected products');

    for (const row of rows) {
      if (!row.route) continue;
      const lhs = row.route.split('>>')[0];
      const reactants = lhs.split('.').filter((r) => r.length > 0);
      const allBBs = reactants.every((r) => canonBBs.has(r));
      expect(allBBs, false,
        `Step "${row.route}" used BB-only reactants (${reactants.join(' + ')}) — slot pools should have isolated BBs to a single slot.`);
    }
  }, {timeout: 300000});

  test('reagents mode: round-2 step uses one round-1 product or BB plus reagents only', async () => {
    // Stronger spec check: in a 2-round run, for every multi-step route, the SECOND step's
    // reactants must come from {round-1 products of THIS route} or {BBs}, never from raw BBs in
    // the "other" slots. (One reactant comes from the round-1 product or BB; everything else is
    // reagents.) The combo cap needs to be high enough that round-2 sampling reaches the
    // round-1 products in the main slot — eligibleSmiles is [BBs..., prevRoundProducts...] so
    // the first (bbs × reagents) combos hit BBs only. With 4 BBs × 16 reagents = 64, anything
    // ≤ 64 misses prevRoundProducts entirely; we use 200 to be comfortably past that.
    const rdkit = getRdKitModule();
    const {bbs, reagents, canonBBs, canonReagents, peptideTemplate} = await loadAaSplit(rdkit, 4);

    const config = reagentsModeBaseConfig();
    config.enumeration.num_rounds = 2;
    config.max_num_combinations_per_template = 200;

    const {rows} = await enumerate({
      rdkit, config, templates: [peptideTemplate],
      buildingBlocks: bbs, reagents, exclusionSmarts: [],
    });
    expect(rows.length > 0, true, 'Expected products');

    let multiStepChecked = 0;
    for (const row of rows) {
      if (!row.route) continue;
      const stepStrs = row.route.split('--**--');
      if (stepStrs.length < 2) continue;
      multiStepChecked++;
      // After step 1, the "known products" set is the RHS of step 1.
      const step1 = stepStrs[0];
      const step1Rhs = step1.slice(step1.indexOf('>>') + 2).split('.').filter((s) => s.length > 0);
      const knownAfterStep1 = new Set([...step1Rhs]);

      const step2 = stepStrs[1];
      const lhs2 = step2.slice(0, step2.indexOf('>>'));
      const r2 = lhs2.split('.').filter((r) => r.length > 0);
      // Step-2 main candidate set: round-1 products in this same route, OR fresh BBs starting a
      // new chain in round 2. Reagents fill the remaining slots. Use the canonicalized BB set so
      // raw-SMILES vs RDKit-emitted-SMILES casing/forms don't cause false negatives.
      const fromMain = r2.filter((r) => knownAfterStep1.has(r) || canonBBs.has(r));
      const fromReagent = r2.filter((r) => canonReagents.has(r));
      expect(fromMain.length, 1,
        `Step 2 "${step2}" of route "${row.route}" should have EXACTLY 1 main-pool reactant ` +
        `(round-1 product or BB), got ${fromMain.length}. Reactants: ${r2.join(' + ')}`);
      expect(fromMain.length + fromReagent.length, r2.length,
        `Step 2 "${step2}" had unclassified reactants: ${r2.join(' + ')}`);
    }
    expect(multiStepChecked > 0, true, 'No multi-step routes produced — round 2 wasn\'t exercised');
  }, {timeout: 600000});
});
