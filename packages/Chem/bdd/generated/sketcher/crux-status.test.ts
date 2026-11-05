/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/crux-status.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/crux.js';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {pressKeyIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {areaAtLeastTall, clickArea, doubleClickArea, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, areaLiesBeside} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Crux in the platform's automation", () => {
  const session = feature(test, "features/sketcher/crux-status.feature", import.meta.url);
  test("Crux's canvas, toolbars and label editor are parts of its widget and areas of its status", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"],["c1ccccc1"]]), [["molecule"],["CCO"],["c1ccccc1"]]);
    await session.step(20, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(21, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(22, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(23, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"CCO\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCO"));
    await session.step(27, "Then canvas of crux sketcher widget should be visible", () => shouldBe(page, el("canvas of crux sketcher widget"), "visible"));
    await session.step(28, "And actions of crux sketcher widget should be visible", () => shouldBe(page, el("actions of crux sketcher widget"), "visible"));
    await session.step(29, "And tools of crux sketcher widget should be visible", () => shouldBe(page, el("tools of crux sketcher widget"), "visible"));
    await session.step(30, "And elements of crux sketcher widget should be visible", () => shouldBe(page, el("elements of crux sketcher widget"), "visible"));
    await session.step(31, "And templates of crux sketcher widget should be visible", () => shouldBe(page, el("templates of crux sketcher widget"), "visible"));
    await session.step(32, "And the \"canvas\" area of crux sketcher widget should lie below the \"actions\" area", () => areaLies(page, "canvas", el("crux sketcher widget"), "below", "actions"));
    await session.step(33, "And the \"canvas\" area of crux sketcher widget should lie to the right of the \"tools\" area", () => areaLiesBeside(page, "canvas", el("crux sketcher widget"), "right", "tools"));
    await session.step(34, "And the \"elements\" area of crux sketcher widget should lie to the right of the \"canvas\" area", () => areaLiesBeside(page, "elements", el("crux sketcher widget"), "right", "canvas"));
    await session.step(35, "And the \"templates\" area of crux sketcher widget should lie below the \"canvas\" area", () => areaLies(page, "templates", el("crux sketcher widget"), "below", "canvas"));
    await session.step(36, "And label editor of crux sketcher widget should be absent", () => shouldBe(page, el("label editor of crux sketcher widget"), "absent"));
    await session.step(37, "When user double-clicks on the \"atom 1\" area of crux sketcher widget", () => doubleClickArea(page, "atom 1", el("crux sketcher widget")));
    await session.step(38, "Then label editor of crux sketcher widget should be visible", () => shouldBe(page, el("label editor of crux sketcher widget"), "visible"));
    await session.step(39, "And the \"label editor\" area of crux sketcher widget should be at least 10 pixels tall", () => areaAtLeastTall(page, "label editor", el("crux sketcher widget"), 10));
  });
  test("A click on a tool's area chooses that tool, and its edit on an atom's area is read back", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"],["c1ccccc1"]]), [["molecule"],["CCO"],["c1ccccc1"]]);
    await session.step(20, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(21, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(22, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(23, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"CCO\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCO"));
    await session.step(43, "Then the \"tool\" reading of crux sketcher widget should be \"bond.single\"", () => readingReads(page, "tool", el("crux sketcher widget"), "bond.single"));
    await session.step(44, "And the \"empty\" reading of crux sketcher widget should be \"false\"", () => readingReads(page, "empty", el("crux sketcher widget"), "false"));
    await session.step(45, "And the \"changes\" reading of crux sketcher widget should be 1", () => readingIs(page, "changes", el("crux sketcher widget"), 1));
    await session.step(46, "When user clicks on the \"tool element.n\" area of crux sketcher widget", () => clickArea(page, "tool element.n", el("crux sketcher widget")));
    await session.step(47, "Then the \"tool\" reading of crux sketcher widget should be \"element.n\"", () => readingReads(page, "tool", el("crux sketcher widget"), "element.n"));
    await session.step(48, "When user clicks on the \"atom 2\" area of crux sketcher widget", () => clickArea(page, "atom 2", el("crux sketcher widget")));
    await session.step(49, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"CCN\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCN"));
    await session.step(50, "And the \"changes\" reading of crux sketcher widget should be 2", () => readingIs(page, "changes", el("crux sketcher widget"), 2));
    await session.step(51, "And the \"query\" reading of crux sketcher widget should be \"false\"", () => readingReads(page, "query", el("crux sketcher widget"), "false"));
  });
  test("The selected atoms are read in molfile order, and a cleared canvas reads as empty", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"Crux\"", () => sketcherIs(page, "Crux"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["CCO"],["c1ccccc1"]]), [["molecule"],["CCO"],["c1ccccc1"]]);
    await session.step(20, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(21, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(22, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(23, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"CCO\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "CCO"));
    await session.step(55, "Then the \"selected atoms\" reading of crux sketcher widget should be \"\"", () => readingReads(page, "selected atoms", el("crux sketcher widget"), ""));
    await session.step(56, "When user presses Control+A in crux canvas", () => pressKeyIn(page, "Control+A", el("crux canvas")));
    await session.step(57, "Then the \"selected atoms\" reading of crux sketcher widget should be \"0, 1, 2\"", () => readingReads(page, "selected atoms", el("crux sketcher widget"), "0, 1, 2"));
    await session.step(58, "And the \"selected bonds\" reading of crux sketcher widget should be \"0, 1\"", () => readingReads(page, "selected bonds", el("crux sketcher widget"), "0, 1"));
    await session.step(59, "And the \"tool\" reading of crux sketcher widget should be \"select.rect\"", () => readingReads(page, "tool", el("crux sketcher widget"), "select.rect"));
    await session.step(60, "When user clicks on the \"tool clear\" area of crux sketcher widget", () => clickArea(page, "tool clear", el("crux sketcher widget")));
    await session.step(61, "Then the \"empty\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "empty", el("crux sketcher widget"), "true"));
    await session.step(62, "And the \"atoms\" reading of crux sketcher widget should be 0", () => readingIs(page, "atoms", el("crux sketcher widget"), 0));
    await session.step(63, "And the \"selected atoms\" reading of crux sketcher widget should be \"\"", () => readingReads(page, "selected atoms", el("crux sketcher widget"), ""));
  });
});
