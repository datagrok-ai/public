// Tests of equations parsing

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from '../package-test';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {TEMPLATES, ENERGY_N_CONTROL} from '../templates';
import {USE_CASES} from '../use-cases';
import {testTemplate} from './test-utils';
import {getIVP} from '../scripting-tools';
import {MAX_LINE_CHART} from '../constants';

// ---------------------------------------------------------------------------
// Full solve-pipeline tests (parse → generate script → execute → verify output)
//
// Each testTemplate() call:
//   1. Runs getIVP(text) to parse the DSL model into a typed IVP object.
//   2. Calls getScriptLines(ivp) to emit a Datagrok JS script string.
//   3. Creates a DG.Script from the emitted text, invokes it via call.call().
//   4. Asserts that the resulting solution table has at least one row.
// ---------------------------------------------------------------------------
category('Features', () => {
  // --- Templates ---

  // Tests the simplest template: one ODE, four required blocks (#name, #equations, #inits, #argument).
  // Baseline check — if this test fails, a core parser or script-runner regression has occurred.
  testTemplate('Basic project', TEMPLATES.BASIC);

  // Tests the Advanced template, which adds #constants and #expressions blocks.
  // Verifies that fixed constants and expression-derived intermediates are correctly substituted in equations.
  testTemplate('Project structs', TEMPLATES.ADVANCED);

  // Tests the Extended template, which adds #parameters with UI annotations (units, min, max, step).
  // Verifies that annotated parameters are parsed and passed as script inputs to the solver.
  testTemplate('Annotating params', TEMPLATES.EXTENDED);

  // --- Library use-cases ---

  // Tests the Chemical Reaction model (coupled ODEs for reactant/product concentrations).
  // Covers the common case of multiple linked ODEs without any special (#loop/#update) blocks.
  testTemplate('Chem react model', USE_CASES.CHEM_REACT);

  // Tests the Robertson stiff ODE benchmark (3 equations with vastly different time scales).
  // Exercises the stiff solver path; the model sets #tolerance to 1e-7 for sufficient accuracy.
  testTemplate('Robertson stiff ODEs', USE_CASES.ROBERTSON);

  // Tests the Fermentation model (substrate, biomass, and product concentration dynamics).
  // Typical biochemical engineering use-case with kinetic constants defined in #parameters.
  testTemplate('Fermentation model', USE_CASES.FERMENTATION);

  // Tests the PK (pharmacokinetics) model, which includes a #loop block for repeated drug dosing.
  // Verifies that cyclic process logic (loop count, dose-reset expressions) runs correctly end-to-end.
  testTemplate('PK model', USE_CASES.PK);

  // Tests the PK-PD model (pharmacokinetics + pharmacodynamics with effect compartment).
  // Also uses a #loop block; confirms the combined PK and PD equations solve without errors.
  testTemplate('Cyclic process', USE_CASES.PK_PD);

  // Tests the Bioreactor model: multi-line equations, #constants, #parameters, #expressions,
  // and more output columns than MAX_LINE_CHART — exercises the multi-axis output path in the generator.
  testTemplate('Bioreactor model', USE_CASES.BIOREACTOR);

  // Tests the Acid Production (ACID_PROD) model with #update blocks (multi-stage process).
  // Verifies the multi-stage solve pipeline: each #update stage modifies state, and a '_Stage'
  // column is appended to the output table.
  testTemplate('Multistage model', USE_CASES.ACID_PROD);

  // Tests the Nimotuzumab model with a custom #output block where no entry has a formula.
  // Verifies that output column selection from raw ODE variable names (not computed expressions) works.
  testTemplate('Output control', USE_CASES.NIMOTUZUMAB);

  // Reuses the Nimotuzumab model to exercise the #meta.inputs value-lookup table path.
  // Verifies that categorical input lookup tables (mapping string keys to numeric values) are handled.
  testTemplate('Value lookups', USE_CASES.NIMOTUZUMAB);

  // Tests the Pollution model (20 coupled ODEs), verifying large-system parsing and solving.
  // Edge case for parser correctness with many state variables and many initial conditions.
  testTemplate('Pollution model', USE_CASES.POLLUTION);

  // Tests ENERGY_N_CONTROL, which uses JavaScript expressions in #expressions and a custom #output.
  // Verifies that JS code embedded in model expressions is correctly emitted and evaluated at runtime.
  testTemplate('Output expressions & use of JS code in model', ENERGY_N_CONTROL);
}); // Features
// ---------------------------------------------------------------------------
// Structural checks derived from manual test observations:
// verify IVP object properties that drive UI/viewer behaviour.
// Each test calls getIVP() directly and inspects fields of the returned IVP object
// without running the full solver — faster checks for structural invariants.
// ---------------------------------------------------------------------------
category('Features: Structure', () => {
  // Tests that Bioreactor's output column count exceeds MAX_LINE_CHART (= 4).
  // getIVP() populates ivp.outputs (or falls back to ivp.inits.size when no #output block exists);
  // getScriptLines() reads this count and injects multiAxis: "true" into the viewer annotation
  // when count > MAX_LINE_CHART — this is what enables the Multiaxis chart tab in the Diff Studio UI.
  test('Bioreactor output columns exceed MAX_LINE_CHART (multiaxis mode expected)', async () => {
    const ivp = getIVP(USE_CASES.BIOREACTOR);
    const outputCount = ivp.outputs ? ivp.outputs.size : ivp.inits.size;
    expect(
      outputCount > MAX_LINE_CHART,
      true,
      `Bioreactor output columns (${outputCount}) must exceed MAX_LINE_CHART (${MAX_LINE_CHART}) to trigger multiaxis`,
    );
  });

  // Tests that the PK-PD model's #loop block is parsed into a non-null ivp.loop object.
  // getIVP() reads the #loop tag and fills Loop { count, updates }; the Diff Studio UI uses
  // ivp.loop to show the Count clicker input. Also verifies mutual exclusion: ivp.updates === null
  // when a #loop block is present (the two block types cannot coexist in the same model).
  test('PK-PD model carries a #loop block (cyclic process)', async () => {
    const ivp = getIVP(USE_CASES.PK_PD);
    expect(ivp.loop !== null, true, 'PK-PD model must have a #loop block');
    expect(ivp.updates === null, true, '#loop and #update must not coexist');
    expect(ivp.loop!.count.value >= 1, true, 'Loop count must be >= 1');
  });

  // Tests that the PK model's #loop block contains at least one state-variable update expression.
  // Unlike PK-PD, PK's loop includes explicit dose-reset formulas (e.g., Ac += dose);
  // ivp.loop.updates.length >= 1 confirms those formulas were parsed and stored correctly.
  test('PK model carries a #loop block (cyclic process)', async () => {
    const ivp = getIVP(USE_CASES.PK);
    expect(ivp.loop !== null, true, 'PK model must have a #loop block');
    expect(ivp.loop!.updates.length >= 1, true, 'Loop must have at least one update expression');
  });

  // Tests that the ACID_PROD (gluconic acid) model's #update blocks are parsed into ivp.updates.
  // Each #update entry defines a stage transition: duration + state mutations applied between stages.
  // The multi-stage solver uses this array to drive sequential stage execution, appending a '_Stage'
  // column to the output. Also verifies mutual exclusion: ivp.loop === null when #update is present.
  test('GA-production (ACID_PROD) model carries #update stages (multistage process)', async () => {
    const ivp = getIVP(USE_CASES.ACID_PROD);
    expect(ivp.updates !== null, true, 'ACID_PROD model must have #update blocks');
    expect(ivp.updates!.length >= 1, true, 'Must have at least one update stage');
    expect(ivp.loop === null, true, '#loop and #update must not coexist');
  });

  // Tests that all entries in the Nimotuzumab #output block have formula === null.
  // When formula is null, getScriptLines() maps the output column name directly to an ODE variable,
  // rather than emitting a derived expression — this controls exactly which columns appear in the
  // output table and the line chart viewers.
  test('Nimotuzumab model has custom #output with no formula entries', async () => {
    const ivp = getIVP(USE_CASES.NIMOTUZUMAB);
    expect(ivp.outputs !== null, true, 'Nimotuzumab model must have an #output block');
    let allNoFormula = true;
    ivp.outputs!.forEach((out) => {if (out.formula !== null) allNoFormula = false;});
    expect(allNoFormula, true, 'All Nimotuzumab #output entries must have formula === null');
  });

  // Tests that ENERGY_N_CONTROL's #output block parses to exactly 3 entries: t, energy, control.
  // These entries reference expression-derived values (not raw ODE variables); the count is used
  // by getScriptLines() to build the output column declaration in the generated script header.
  test('ENERGY_N_CONTROL model has custom #output referencing expression-based columns', async () => {
    const ivp = getIVP(ENERGY_N_CONTROL);
    expect(ivp.outputs !== null, true, 'ENERGY_N_CONTROL must have an #output block');
    expect(ivp.outputs!.size, 3, 'ENERGY_N_CONTROL #output must define exactly 3 columns (t, energy, control)');
  });

  // Tests that the Pollution model contains exactly 20 differential equations in ivp.deqs.equations.
  // This large-system model stresses the parser's equation-collection loop; verifies all 20 ODEs
  // are stored in the Map without truncation, merging errors, or off-by-one mistakes.
  test('Pollution model has 20 differential equations', async () => {
    const ivp = getIVP(USE_CASES.POLLUTION);
    expect(ivp.deqs.equations.size, 20, 'Pollution model must have exactly 20 differential equations');
  });

  // Tests that the Bioreactor model's three supporting block types are all parsed correctly:
  //   ivp.consts  ← #constants  (fixed numeric values, not exposed as UI inputs or sliders)
  //   ivp.params  ← #parameters (mutable values shown as sliders / numeric inputs in the UI)
  //   ivp.exprs   ← #expressions (computed intermediate variables referenced by the ODEs)
  // Their co-presence exercises the full block-parsing path in getIVP(), including
  // the order-of-processing rules between constants, parameters, and expressions.
  test('Bioreactor model carries #constants, #parameters, and #expressions blocks', async () => {
    const ivp = getIVP(USE_CASES.BIOREACTOR);
    expect(ivp.consts !== null, true, 'Bioreactor must have a #constants block');
    expect(ivp.params !== null, true, 'Bioreactor must have a #parameters block');
    expect(ivp.exprs !== null, true, 'Bioreactor must have an #expressions block');
  });

  // Tests that the Robertson model's #tolerance block is stored as ivp.tolerance = '1e-7'.
  // Robertson is a benchmark stiff system where tight tolerance is essential for correct solution;
  // the string value is passed as-is to the adaptive solver (ROS3PRw / ROS34PRw) at runtime,
  // so an exact string match is required — not just numeric equivalence.
  test('Robertson model carries #tolerance set to 1e-7', async () => {
    const ivp = getIVP(USE_CASES.ROBERTSON);
    expect(ivp.tolerance, '1e-7', 'Robertson model tolerance must be 1e-7');
  });
}); // Features: Structure

