/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/crux-choice.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {sketcherChosen} from '../../bindings/crux.js';
import {clipboardMolecule} from '../../bindings/molecules.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, openTableOf, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, readingIsMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {doubleClickArea, menuLists, noErrors, pickFromOpenMenu, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Crux as a sketcher choice", () => {
  const session = feature(test, "features/sketcher/crux-choice.feature", import.meta.url);
  test("The sketcher's options menu offers Crux among the sketchers, OpenChemLib checked", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(22, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(26, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(27, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(28, "And crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(29, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(30, "Then the open menu should list \"Crux\"", () => menuLists(page, "Crux"));
    await session.step(31, "And \"OpenChemLib\" menu item should be selected", () => shouldBe(page, el("\"OpenChemLib\" menu item"), "selected"));
    await session.step(32, "And \"Crux\" menu item should not be selected", () => shouldNotBe(page, el("\"Crux\" menu item"), "selected"));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("Picking Crux in the options menu makes it the account's sketcher, also after a reload", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(22, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(37, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(38, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(39, "And crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(40, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(41, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(42, "Then the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(43, "And the session's sketcher and the account's choice should be \"Crux\"", () => sketcherChosen(page, "Crux"));
    await session.step(44, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(45, "And user reloads the page", () => reloadPage(page));
    await session.step(46, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"]]);
    await session.step(50, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(51, "And user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(52, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(53, "And the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(54, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
  });
  test("Switching sketchers in place keeps the molecule, its stereo included", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(15, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(16, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(22, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(58, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(59, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(60, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(61, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(62, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
    await session.step(63, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(64, "And user picks \"Ketcher\" from the open menu", () => pickFromOpenMenu(page, "Ketcher"));
    await session.step(65, "Then crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(66, "And Ketcher canvas should be visible", () => shouldBe(page, el("Ketcher canvas"), "visible"));
    await session.step(67, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(68, "And user picks \"Copy as MOLBLOCK\" from the open menu", () => pickFromOpenMenu(page, "Copy as MOLBLOCK"));
    await session.step(69, "Then the clipboard should hold the molecule of row 1 of \"molecule\" column", () => clipboardMolecule(page, 1, "molecule"));
    await session.step(70, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(71, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(72, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
  });
});
