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
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(29, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(30, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(31, "And crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(32, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(33, "Then the open menu should list \"Crux\"", () => menuLists(page, "Crux"));
    await session.step(34, "And \"OpenChemLib\" menu item should be selected", () => shouldBe(page, el("\"OpenChemLib\" menu item"), "selected"));
    await session.step(35, "And \"Crux\" menu item should not be selected", () => shouldNotBe(page, el("\"Crux\" menu item"), "selected"));
    await session.step(36, "And no errors should have been logged", () => noErrors(page));
  });
  test("Picking Crux in the options menu makes it the account's sketcher, also after a reload", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(40, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(41, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(42, "And crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(43, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(44, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(45, "Then the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(46, "And the session's sketcher and the account's choice should be \"Crux\"", () => sketcherChosen(page, "Crux"));
    await session.step(47, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(48, "And user reloads the page", () => reloadPage(page));
    await session.step(49, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"]]);
    await session.step(53, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(54, "And user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(55, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(56, "And the \"ready\" reading of crux sketcher widget should be \"true\"", () => readingReads(page, "ready", el("crux sketcher widget"), "true"));
    await session.step(57, "And the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
  });
  test("Switching sketchers in place keeps the molecule, its stereo included", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(61, "When user double-clicks on the \"cell 1 of molecule\" area of grid", () => doubleClickArea(page, "cell 1 of molecule", el("grid")));
    await session.step(62, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(63, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(64, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(65, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
    await session.step(66, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(67, "And user picks \"Ketcher\" from the open menu", () => pickFromOpenMenu(page, "Ketcher"));
    await session.step(68, "Then crux sketcher widget should be absent", () => shouldBe(page, el("crux sketcher widget"), "absent"));
    await session.step(69, "And Ketcher canvas should be visible", () => shouldBe(page, el("Ketcher canvas"), "visible"));
    await session.step(70, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(71, "And user picks \"Copy as MOLBLOCK\" from the open menu", () => pickFromOpenMenu(page, "Copy as MOLBLOCK"));
    await session.step(72, "Then the clipboard should hold the molecule of row 1 of \"molecule\" column", () => clipboardMolecule(page, 1, "molecule"));
    await session.step(73, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(74, "And user picks \"Crux\" from the open menu", () => pickFromOpenMenu(page, "Crux"));
    await session.step(75, "Then the \"smiles\" reading of crux sketcher widget should be the molecule \"C[C@H](N)C(=O)O\"", () => readingIsMolecule(page, "smiles", el("crux sketcher widget"), "C[C@H](N)C(=O)O"));
    await session.step(76, "And no errors should have been logged", () => noErrors(page));
  });
  test("Ketcher picked in the filter card's sketcher is the session's and the account's, and the card's sketcher reopens after the switch back", {tag: ["@sketcher-controls"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(19, "And user opens a table \"molecules\" with:", () => openTableOf(page, "molecules", [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]), [["molecule"],["C[C@H](N)C(=O)O"],["c1ccccc1"],["CCO"],["CC(=O)O"]]);
    await session.step(25, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(80, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
    await session.step(81, "Then \"molecule\" filter card should be visible", () => shouldBe(page, el("\"molecule\" filter card"), "visible"));
    await session.step(82, "When user clicks on \"Sketch\" text in \"molecule\" filter card", () => clickOn(page, el("\"Sketch\" text in \"molecule\" filter card")));
    await session.step(83, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(84, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(85, "And user picks \"Ketcher\" from the open menu", () => pickFromOpenMenu(page, "Ketcher"));
    await session.step(86, "Then Ketcher canvas should be visible", () => shouldBe(page, el("Ketcher canvas"), "visible"));
    await session.step(87, "And the session's sketcher and the account's choice should be \"Ketcher\"", () => sketcherChosen(page, "Ketcher"));
    await session.step(88, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(89, "And user picks \"OpenChemLib\" from the open menu", () => pickFromOpenMenu(page, "OpenChemLib"));
    await session.step(90, "Then the session's sketcher and the account's choice should be \"OpenChemLib\"", () => sketcherChosen(page, "OpenChemLib"));
    await session.step(91, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
    await session.step(92, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
    await session.step(93, "When user clicks on \"Sketch\" text in \"molecule\" filter card", () => clickOn(page, el("\"Sketch\" text in \"molecule\" filter card")));
    await session.step(94, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
    await session.step(95, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
    await session.step(96, "Then \"OpenChemLib\" menu item should be selected", () => shouldBe(page, el("\"OpenChemLib\" menu item"), "selected"));
    await session.step(97, "And \"Ketcher\" menu item should not be selected", () => shouldNotBe(page, el("\"Ketcher\" menu item"), "selected"));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
});
