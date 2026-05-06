import * as DG from 'datagrok-api/dg';
import {before, category, test, expect} from '@datagrok-libraries/test/src/test';
import {_package} from '../package-test';
import {cloneConfig, DEFAULT_CONFIG} from '../utils/reaction-enumeration/config';
import {
  enumerate, formatRoute, Route,
  splitSmartsByReactants, stripOuterParens, TemplateInput,
} from '../utils/reaction-enumeration/enumerate';
import * as chemCommonRdKit from '../utils/chem-common-rdkit';
import {getRdKitModule} from '../utils/chem-common-rdkit';

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

// Mirror of packages/Chem/src/rendering/rdkit-reaction-renderer.ts:parseMultiStepReaction.
// Used to verify our route format renders correctly.
function parseMultiStepReaction(reactionString: string): string[] {
  const parts = reactionString.split('>>');
  if (parts.length <= 2) return [reactionString];
  const steps: string[] = [];
  for (let i = 0; i < parts.length - 1; i++) {
    const left = parts[i].trim();
    const right = parts[i + 1].trim();
    if (left && right) steps.push(`${left}>>${right}`);
  }
  return steps.length > 0 ? steps : [reactionString];
}

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

  test('route formatting: two-step linear chain — Chem-renderer-compatible format', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1', 'BB3'], 'P2'),
    ];
    const formatted = formatRoute(route);
    expect(formatted, 'BB1.BB2>>P1.BB3>>P2');
    const steps = parseMultiStepReaction(formatted);
    expect(steps.length, 2);
    expect(steps[0], 'BB1.BB2>>P1.BB3');
    expect(steps[1], 'P1.BB3>>P2');
  });

  test('route formatting: three-step linear chain', async () => {
    const route: Route = [
      step(['A', 'B'], 'P1'),
      step(['P1', 'C'], 'P2'),
      step(['P2', 'D'], 'P3'),
    ];
    expect(formatRoute(route), 'A.B>>P1.C>>P2.D>>P3');
  });

  test('route formatting: linear chain where the next step has only the prev product as reactant', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1'], 'P2'),
    ];
    expect(formatRoute(route), 'BB1.BB2>>P1>>P2');
  });

  test('route formatting: branching route concatenates with >> only — never uses semicolons', async () => {
    const route: Route = [
      step(['A', 'B'], 'P1'),
      step(['C', 'D'], 'P2'),
      step(['P1', 'P2'], 'P3'),
    ];
    const formatted = formatRoute(route);
    expect(formatted.indexOf('; '), -1, `Route used "; " separator: ${formatted}`);
    expect(formatted.indexOf(';'), -1, `Route used ";" separator: ${formatted}`);
    const steps = parseMultiStepReaction(formatted);
    expect(steps.length >= 1, true, `parseMultiStepReaction returned no steps for ${formatted}`);
  });

  test('route formatting: each step is a valid single-step reaction after parseMultiStepReaction', async () => {
    const route: Route = [
      step(['BB1', 'BB2'], 'P1'),
      step(['P1', 'BB3'], 'P2'),
    ];
    const formatted = formatRoute(route);
    const steps = parseMultiStepReaction(formatted);
    for (const s of steps) {
      const arrowCount = s.split('>>').length - 1;
      expect(arrowCount, 1, `Step "${s}" has ${arrowCount} ">>" separators (expected 1)`);
    }
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
});
