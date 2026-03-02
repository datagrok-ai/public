// Tests covering IVP parser correctness, negative cases, and script generation

import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';
import {getIVP, getScriptLines, getScriptParams} from '../scripting-tools';
import {ModelError} from '../error-utils';
import {TEMPLATES, ENERGY_N_CONTROL} from '../templates';
import {USE_CASES} from '../use-cases';
import {MAX_LINE_CHART} from '../constants';

// Helpers

/** Minimal valid IVP text used as a base for negative-test mutations */
const MINIMAL_VALID = `#name: MinimalModel
#equations:
  dy/dt = -y

#inits:
  y = 1

#argument: t
  initial = 0
  final = 10
  step = 0.1`;

/** Replace a section in the base model text */
function withoutBlock(text: string, tag: string): string {
  const lines = text.split('\n');
  let inside = false;
  return lines.filter((line) => {
    if (line.trimStart().startsWith(tag)) {inside = true; return false;}
    if (inside && line.trimStart().startsWith('#')) inside = false;
    return !inside;
  }).join('\n');
}

/** Return true if fn() throws a ModelError */
function throwsModelError(fn: () => void): boolean {
  try {fn(); return false;} catch (e) {return e instanceof ModelError;}
}

// Negative tests — parser must throw ModelError

category('Parser: Negative', () => {
  // Tests that getIVP() throws ModelError when the #name block is absent.
  // checkCorrectness() in scripting-tools.ts validates all required block presence before
  // building the IVP object; #name is required to identify the model in the UI and script output.
  test('TC-N-01: Missing #name block', async () => {
    const text = withoutBlock(MINIMAL_VALID, '#name');
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for missing #name');
  });

  // Tests that getIVP() throws ModelError when the #equations block is missing.
  // Without #equations there are no ODEs to solve — the model is structurally incomplete.
  // checkCorrectness() enforces that at least one differential equation must be present.
  test('TC-N-02: Missing #equations block', async () => {
    const text = withoutBlock(MINIMAL_VALID, '#equations');
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for missing #equations');
  });

  // Tests that getIVP() throws ModelError when the #inits block is absent.
  // #inits provides the initial conditions (y(t0) values) for every ODE state variable;
  // without them the solver cannot be initialised and the IVP is undefined.
  test('TC-N-03: Missing #inits block', async () => {
    const text = withoutBlock(MINIMAL_VALID, '#inits');
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for missing #inits');
  });

  // Tests that getIVP() throws ModelError when the #argument block (integration interval) is absent.
  // The argument block defines t0 (initial), t_final (final), and step — all required to
  // construct the time grid and call the solver.
  test('TC-N-04: Missing #argument block', async () => {
    const text = withoutBlock(MINIMAL_VALID, '#argument');
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for missing #argument');
  });

  // Tests that getIVP() throws ModelError when the argument interval is reversed (initial >= final).
  // checkCorrectness() validates that initial < final; an empty or reversed interval makes
  // the integration meaningless and would cause the solver to run zero or negative steps.
  test('TC-N-05: initial >= final in #argument', async () => {
    const text = `#name: BadRange
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 10
  final = 0
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for initial >= final');
  });

  // Tests that getIVP() throws ModelError when the step value is <= 0.
  // A non-positive step would cause the solver to run indefinitely or produce an empty result;
  // checkCorrectness() rejects non-positive step values as structurally invalid.
  test('TC-N-06: Non-positive step in #argument', async () => {
    const text = `#name: BadStep
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = -0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for step <= 0');
  });

  // Tests that getIVP() throws ModelError when the step is larger than (final - initial).
  // A step exceeding the interval width would produce at most one solution point, which is
  // not a valid numerical integration. checkCorrectness() rejects this configuration.
  test('TC-N-07: Step exceeds interval length', async () => {
    const text = `#name: BigStep
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 0
  final = 1
  step = 5`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for step > interval');
  });

  // Tests that getIVP() throws ModelError for variable names beginning with '_'.
  // The '_' prefix is reserved for internally generated output columns (e.g., '_Stage');
  // checkCorrectness() prevents user-defined variables from colliding with this namespace.
  test('TC-N-08: Variable name starts with _ (reserved prefix)', async () => {
    const text = `#name: BadVarName
#equations:
  d(_y)/dt = -_y
#inits:
  _y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for variable starting with _');
  });

  // Tests that getIVP() throws ModelError when two variable names differ only by case (x vs X).
  // The parser normalises names for lookup; duplicate names (case-insensitive) create ambiguity
  // in the equation system and are rejected by checkCorrectness().
  test('TC-N-09: Duplicate variable name (case-insensitive) in #inits', async () => {
    const text = `#name: DuplicateVar
#equations:
  dx/dt = -x
  dX/dt = -X
#inits:
  x = 1
  X = 2
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for duplicate variable name');
  });

  // Tests that getIVP() throws ModelError when both #loop and #update blocks are present.
  // These two blocks represent mutually exclusive process types:
  //   #loop  = a single-stage cyclic process repeated N times
  //   #update = a multi-stage process with different conditions per stage
  // checkCorrectness() enforces that only one can appear in any given model.
  test('TC-N-10: Simultaneous #loop and #update blocks', async () => {
    const text = `#name: LoopAndUpdate
#equations:
  dy/dt = -y
#inits:
  y = 1
#loop:
  count = 2
  y += 0.1
#update: stage2
  duration = 5
  y += 0.1
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for simultaneous #loop and #update');
  });

  // Tests that getIVP() throws ModelError when a variable defined in #equations has no
  // corresponding entry in #inits. Every ODE state variable requires an initial condition;
  // a missing one prevents the solver from constructing its initial state vector.
  test('TC-N-11: Missing initial value for equation variable', async () => {
    const text = `#name: MissingInit
#equations:
  dy/dt = -y
  dz/dt = -z
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for missing initial value for z');
  });

  // Tests that getIVP() throws ModelError when a #loop block specifies count = 0.
  // Zero iterations is semantically meaningless for a cyclic process;
  // checkCorrectness() requires count >= 1 for any loop block.
  test('TC-N-12: Loop count less than 1', async () => {
    const text = `#name: BadCount
#equations:
  dy/dt = -y
#inits:
  y = 1
#loop:
  count = 0
  y += 0.1
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for loop count < 1');
  });

  // Tests that getIVP() throws ModelError when a #meta.solver option exceeds its allowed range.
  // SOLVER_OPTIONS_RANGES in constants.ts defines bounds (e.g., maxTime: 1–10000);
  // checkCorrectness() validates all parsed solver option values against these bounds and throws
  // ModelError for any out-of-range value.
  test('TC-N-13: Solver settings — maxTime out of allowed range', async () => {
    const text = `#name: BadSolverOpts
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1
#meta.solver: {maxTime: 99999}`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for maxTime out of range');
  });

  // Tests that getIVP() throws ModelError when #meta.solver content is not enclosed in curly braces.
  // The parser requires solver settings to be a JSON-like object: #meta.solver: { key: value };
  // plain key-value text without braces is not a valid format and must be rejected.
  test('TC-N-14: Solver settings without curly braces', async () => {
    const text = `#name: NoBraces
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1
#meta.solver: maxTime: 500`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for solver settings without braces');
  });

  // Tests that getIVP() throws ModelError when a value in #inits is a non-numeric string.
  // Initial conditions must be parseable as floating-point numbers; string values like 'abc'
  // cannot be used as ODE starting points and are rejected during checkCorrectness() validation.
  test('TC-N-15: Non-numeric value in #inits', async () => {
    const text = `#name: NonNumericInit
#equations:
  dy/dt = -y
#inits:
  y = abc
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    expect(throwsModelError(() => getIVP(text)), true, 'Expected ModelError for non-numeric value in #inits');
  });
}); // Parser: Negative

// Structural tests — IVP object properties after successful parsing

category('Parser: Structure', () => {
  // Tests that a model containing only the four required blocks (#name, #equations, #inits, #argument)
  // parses without throwing and returns an IVP with a non-empty name, at least one equation,
  // and null loop/updates (no cyclic or multi-stage blocks). This is the minimal passing baseline;
  // if it fails, all other positive tests are suspect.
  test('TC-S-01: Minimal valid model parses without error', async () => {
    let ivp = null;
    try {ivp = getIVP(MINIMAL_VALID);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse of minimal model');
    expect(ivp!.name.length > 0, true, 'Expected non-empty model name');
    expect(ivp!.deqs.equations.size >= 1, true, 'Expected at least one equation');
    expect(ivp!.loop === null, true, 'Expected no loop block');
    expect(ivp!.updates === null, true, 'Expected no update blocks');
  });

  // Tests that the PK-PD model (pharmacokinetics + pharmacodynamics) parses to an IVP with
  // a non-null ivp.loop object. getIVP() reads the #loop tag and fills Loop { count, updates };
  // the UI uses ivp.loop to display the Count clicker input for repeated-dosing simulations.
  // Also verifies mutual exclusion: ivp.updates must be null when a loop block is present.
  test('TC-S-02: PK-PD model has a valid #loop block', async () => {
    let ivp = null;
    try {ivp = getIVP(USE_CASES.PK_PD);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse of PK-PD model');
    expect(ivp!.loop !== null, true, 'Expected loop block to be present');
    expect(ivp!.loop!.count.value >= 1, true, 'Expected loop count >= 1');
    expect(ivp!.updates === null, true, 'Expected no update blocks alongside #loop');
  });

  // Tests that the ACID_PROD (gluconic acid production) model parses to an IVP with a non-null
  // ivp.updates array. Each #update entry defines a stage transition: duration + state mutations.
  // The multi-stage solver iterates over this array, applying transitions between stages and
  // appending a '_Stage' column to the output. Verifies mutual exclusion with #loop.
  test('TC-S-03: GA-production model has valid #update blocks', async () => {
    let ivp = null;
    try {ivp = getIVP(USE_CASES.ACID_PROD);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse of GA-production model');
    expect(ivp!.updates !== null, true, 'Expected update blocks to be present');
    expect(ivp!.updates!.length >= 1, true, 'Expected at least one update stage');
    expect(ivp!.loop === null, true, 'Expected no loop block alongside #update');
  });

  // Tests that the Bioreactor model's output column count exceeds MAX_LINE_CHART (= 4).
  // getIVP() populates ivp.outputs (or falls back to ivp.inits.size when #output is absent);
  // this count is the gating condition in getScriptLines() for injecting multiAxis: "true"
  // into the viewer annotation — which is what enables the Multiaxis plot tab in the Diff Studio UI.
  test('TC-S-04: Bioreactor output column count exceeds MAX_LINE_CHART', async () => {
    let ivp = null;
    try {ivp = getIVP(USE_CASES.BIOREACTOR);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse of Bioreactor model');
    const outputCount = ivp!.outputs ? ivp!.outputs.size : ivp!.inits.size;
    expect(outputCount > MAX_LINE_CHART, true,
      `Expected output columns (${outputCount}) > MAX_LINE_CHART (${MAX_LINE_CHART}) to trigger multiaxis mode`);
  });

  // Tests that // comments (both standalone lines and inline after expressions) are stripped
  // by the pre-processing stage of getIVP() before blocks are parsed.
  // Without comment stripping, comment text could corrupt equation RHS strings or block detection.
  test('TC-S-05: Comments are stripped from the model', async () => {
    const text = `#name: CommentModel
// This entire line is a comment and must be ignored
#equations:
  dy/dt = -y // inline comment should be ignored
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1`;
    let ivp = null;
    try {ivp = getIVP(text);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse with inline comments');
    expect(ivp!.deqs.equations.size, 1, 'Expected exactly 1 equation after comment stripping');
  });

  // Tests that multi-line equation definitions in the Bioreactor model are correctly concatenated
  // into single right-hand side strings. The parser merges continuation lines until the next
  // block-level or top-level line begins; verifies no RHS is empty or null after merging.
  test('TC-S-06: Multi-line formula is concatenated correctly', async () => {
    const ivp = getIVP(USE_CASES.BIOREACTOR);
    // Bioreactor has multi-line equations; verify they are merged (no empty right-hand sides)
    let allValid = true;
    ivp.deqs.equations.forEach((rhs) => {if (!rhs || rhs.trim() === '') allValid = false;});
    expect(allValid, true, 'Expected all equation right-hand sides to be non-empty after multi-line concat');
  });

  // Tests that a #tolerance block is parsed and its value stored in ivp.tolerance as a string.
  // The tolerance value is passed directly to the adaptive solver (ROS3PRw / ROS34PRw) at runtime;
  // an exact string match is required so the solver receives the correct precision setting.
  test('TC-S-07: #tolerance block sets the correct value', async () => {
    const text = `#name: TolModel
#equations:
  dy/dt = -y
#inits:
  y = 1
#argument: t
  initial = 0
  final = 10
  step = 0.1
#tolerance: 1e-7`;
    let ivp = null;
    try {ivp = getIVP(text);} catch (e) {}
    expect(ivp !== null, true, 'Expected successful parse');
    expect(ivp!.tolerance, '1e-7', 'Expected tolerance to be "1e-7"');
  });

  // Tests that a #constants block is parsed and stored in ivp.consts as a Map with at least
  // one entry. Constants differ from parameters: they are fixed values not exposed as UI sliders,
  // but are substituted into equation expressions before the solver runs.
  test('TC-S-08: #constants block is parsed into ivp.consts', async () => {
    const ivp = getIVP(TEMPLATES.ADVANCED);
    expect(ivp.consts !== null, true, 'Expected consts block to be parsed');
    expect(ivp.consts!.size >= 1, true, 'Expected at least one constant');
  });

  // Tests that the Nimotuzumab model's #output block is parsed with formula === null for all entries.
  // When formula is null, getScriptLines() uses the output column name as a direct ODE variable
  // reference, rather than emitting a computed expression. This controls the exact set of columns
  // placed in the output table and shown in the line chart viewers.
  test('TC-S-09: Nimotuzumab #output block has no formula entries', async () => {
    const ivp = getIVP(USE_CASES.NIMOTUZUMAB);
    expect(ivp.outputs !== null, true, 'Expected outputs block to be present');
    let allNoFormula = true;
    ivp.outputs!.forEach((out) => {if (out.formula !== null) allNoFormula = false;});
    expect(allNoFormula, true, 'Expected all #output entries to have formula === null');
  });

  // Tests that the ENERGY_N_CONTROL model's #output block produces exactly 3 entries: t, energy, control.
  // These entries reference expression-derived columns (not raw ODE variables);
  // the exact count is used by getScriptLines() to build the output column list in the script header.
  test('TC-S-10: ENERGY_N_CONTROL #output entries match block definition', async () => {
    const ivp = getIVP(ENERGY_N_CONTROL);
    expect(ivp.outputs !== null, true, 'Expected outputs block to be present');
    // #output block defines: t, energy, control  => 3 entries
    expect(ivp.outputs!.size, 3, 'Expected 3 output entries for ENERGY_N_CONTROL');
  });
}); // Parser: Structure

// Script generation content tests

category('Script: Generation', () => {
  // Tests that getScriptLines() always emits //language: javascript as the first annotation.
  // This annotation is required by the Datagrok script runner to identify the script language;
  // without it, the runner would reject the generated script entirely.
  test('TC-G-01: Generated script contains language annotation', async () => {
    const ivp = getIVP(TEMPLATES.BASIC);
    const lines = getScriptLines(ivp);
    const hasLang = lines.some((l) => l.includes('language: javascript'));
    expect(hasLang, true, 'Expected //language: javascript in generated script');
  });

  // Tests that getScriptLines() injects multiAxis: "true" into the line-chart viewer annotation
  // for the Bioreactor model, whose output column count exceeds MAX_LINE_CHART (= 4).
  // This flag causes the Datagrok line-chart viewer to use separate Y-axes for each output column,
  // which is the expected visual layout for models with many output curves.
  test('TC-G-02: multiAxis "true" for model with > MAX_LINE_CHART output columns', async () => {
    const ivp = getIVP(USE_CASES.BIOREACTOR);
    const script = getScriptLines(ivp).join('\n');
    const hasMultiAxis = script.includes('multiAxis: "true"');
    expect(hasMultiAxis, true,
      `Expected multiAxis: "true" in viewer annotation for Bioreactor (${ivp.outputs!.size} columns)`);
  });

  // Tests that getScriptLines() injects multiAxis: "false" for a model with few output columns.
  // The Basic template produces only 2 columns (t and y), which is <= MAX_LINE_CHART;
  // the generator must explicitly set multiAxis: "false" to keep a single shared Y-axis.
  test('TC-G-03: multiAxis "false" for model with <= MAX_LINE_CHART output columns', async () => {
    const ivp = getIVP(TEMPLATES.BASIC); // outputs: t, y => 2 columns
    const script = getScriptLines(ivp).join('\n');
    const hasMultiAxisFalse = script.includes('multiAxis: "false"');
    expect(hasMultiAxisFalse, true,
      'Expected multiAxis: "false" in viewer annotation for Basic template (2 columns)');
  });

  // Tests that getScriptParams() includes a '_count' key for the PK model, which has a #loop block.
  // '_count' is the generated parameter name for the cyclic loop iteration count; it is emitted
  // as a script parameter so the Diff Studio UI renders it as a Count input/slider for the user.
  test('TC-G-04: Script params contain _count for a cyclic model', async () => {
    const ivp = getIVP(USE_CASES.PK);
    const params = getScriptParams(ivp);
    expect('_count' in params, true, 'Expected _count key in script params for PK model with #loop');
  });

  // Tests that getScriptParams() does NOT emit a '_count' key for a non-cyclic model.
  // The Basic template has no #loop block; '_count' would be meaningless and must not appear
  // among the script parameters, as it would create an unexpected slider in the UI.
  test('TC-G-05: Script params do not contain _count for a non-cyclic model', async () => {
    const ivp = getIVP(TEMPLATES.BASIC);
    const params = getScriptParams(ivp);
    expect('_count' in params, false, 'Expected no _count key in script params for Basic template');
  });

  // Tests that getScriptLines() emits //meta.runOnOpen: true for the Extended template.
  // This Datagrok annotation causes the script runner to execute the script automatically
  // when the user opens it, showing pre-computed results without requiring a manual run.
  test('TC-G-06: Generated script contains meta.runOnOpen annotation', async () => {
    const ivp = getIVP(TEMPLATES.EXTENDED);
    const script = getScriptLines(ivp).join('\n');
    const hasMeta = script.includes('meta.runOnOpen: true');
    expect(hasMeta, true, 'Expected meta.runOnOpen: true in generated script for Extended template');
  });

  // Tests that getScriptParams() returns all expected parameter keys for the Extended template.
  // The Extended template defines:
  //   '_t0', '_t1', '_h'  ← from #argument (integration interval controls)
  //   'x', 'y'            ← from #inits (initial conditions exposed as inputs)
  //   'P1', 'P2'          ← from #parameters (user-tunable model parameters)
  // This verifies the generator correctly aggregates all user-visible inputs into the params map.
  test('TC-G-07: All script params keys are present for Extended template', async () => {
    const ivp = getIVP(TEMPLATES.EXTENDED);
    const params = getScriptParams(ivp);
    // Extended has: _t0, _t1, _h (arg), x, y (inits), P1, P2 (params)
    const requiredKeys = ['_t0', '_t1', '_h', 'x', 'y', 'P1', 'P2'];
    const missingKeys = requiredKeys.filter((k) => !(k in params));
    expect(missingKeys.length, 0,
      `Missing script param keys for Extended template: ${missingKeys.join(', ')}`);
  });
}); // Script: Generation
