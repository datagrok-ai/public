/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sequence-translator/markush-enumerator-templates.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.menu.transform.markush-enumeration]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnSemType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, dialogCloses, openTableOf, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {detectTypes, rowMolecule} from '@datagrok-libraries/bdd/bindings/tiers/molecules/molecules';
import {noBalloons, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Markush Enumerator: R-group templates on every position of a stereo alkene core", () => {
  const session = feature(test, "features/sequence-translator/markush-enumerator-templates.feature", import.meta.url);
  test("Templates on all four R-groups of a stereo alkene enumerate every molecule with its geometry", {tag: ["@realizes:chem.menu.transform.markush-enumeration"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(18, "And user opens a table \"alkene\" with:", () => openTableOf(page, "alkene", [["Core"],["[*:1]/C([*:2])=C(\\[*:3])[*:4]"]]), [["Core"],["[*:1]/C([*:2])=C(\\[*:3])[*:4]"]]);
    await session.step(21, "And the semantic types of the current table are detected", () => detectTypes(page));
    await session.step(24, "Then \"Core\" column should have semantic type \"Molecule\"", () => columnSemType(page, "Core", "Molecule"));
    await session.step(25, "When user picks \"Enumerate Markush Structure...\" from the context menu of the \"cell 1 of Core\" area of grid", () => pickFromAreaContextMenu(page, "Enumerate Markush Structure...", "cell 1 of Core", el("grid")));
    await session.step(26, "Then \"Markush Enumerator\" dialog should be visible", () => shouldBe(page, el("\"Markush Enumerator\" dialog"), "visible"));
    await session.step(27, "And \"Markush Enumerator\" dialog should contain the text \"R4\"", () => shouldContainText(page, el("\"Markush Enumerator\" dialog"), "R4"));
    await session.step(28, "When user clicks on \"Employ templates\" button in \"Markush Enumerator\" dialog", () => clickOn(page, el("\"Employ templates\" button in \"Markush Enumerator\" dialog")));
    await session.step(29, "Then \"Add R-group template\" dialog should be visible", () => shouldBe(page, el("\"Add R-group template\" dialog"), "visible"));
    await session.step(30, "And \"Template\" input in \"Add R-group template\" dialog should have the value \"Hydrogen and Alkyl (C1-C8)\"", () => shouldHaveValue(page, el("\"Template\" input in \"Add R-group template\" dialog"), "Hydrogen and Alkyl (C1-C8)"));
    await session.step(31, "And \"Target R#\" input in \"Add R-group template\" dialog should have the value \"1\"", () => shouldHaveValue(page, el("\"Target R#\" input in \"Add R-group template\" dialog"), "1"));
    await session.step(32, "And \"Add R-group template\" dialog should contain the text \"Add 9 substituents to R1 (new slot).\"", () => shouldContainText(page, el("\"Add R-group template\" dialog"), "Add 9 substituents to R1 (new slot)."));
    await session.step(33, "When user clicks on OK button in \"Add R-group template\" dialog", () => clickOn(page, el("OK button in \"Add R-group template\" dialog")));
    await session.step(34, "Then the \"Add R-group template\" dialog should close", () => dialogCloses(page, "Add R-group template"));
    await session.step(35, "When user clicks on \"Employ templates\" button in \"Markush Enumerator\" dialog", () => clickOn(page, el("\"Employ templates\" button in \"Markush Enumerator\" dialog")));
    await session.step(36, "Then \"Target R#\" input in \"Add R-group template\" dialog should have the value \"2\"", () => shouldHaveValue(page, el("\"Target R#\" input in \"Add R-group template\" dialog"), "2"));
    await session.step(37, "When user clicks on OK button in \"Add R-group template\" dialog", () => clickOn(page, el("OK button in \"Add R-group template\" dialog")));
    await session.step(38, "Then the \"Add R-group template\" dialog should close", () => dialogCloses(page, "Add R-group template"));
    await session.step(39, "When user clicks on \"Employ templates\" button in \"Markush Enumerator\" dialog", () => clickOn(page, el("\"Employ templates\" button in \"Markush Enumerator\" dialog")));
    await session.step(40, "Then \"Target R#\" input in \"Add R-group template\" dialog should have the value \"3\"", () => shouldHaveValue(page, el("\"Target R#\" input in \"Add R-group template\" dialog"), "3"));
    await session.step(41, "When user clicks on OK button in \"Add R-group template\" dialog", () => clickOn(page, el("OK button in \"Add R-group template\" dialog")));
    await session.step(42, "Then the \"Add R-group template\" dialog should close", () => dialogCloses(page, "Add R-group template"));
    await session.step(43, "When user clicks on \"Employ templates\" button in \"Markush Enumerator\" dialog", () => clickOn(page, el("\"Employ templates\" button in \"Markush Enumerator\" dialog")));
    await session.step(44, "Then \"Target R#\" input in \"Add R-group template\" dialog should have the value \"4\"", () => shouldHaveValue(page, el("\"Target R#\" input in \"Add R-group template\" dialog"), "4"));
    await session.step(45, "When user clicks on OK button in \"Add R-group template\" dialog", () => clickOn(page, el("OK button in \"Add R-group template\" dialog")));
    await session.step(46, "Then the \"Add R-group template\" dialog should close", () => dialogCloses(page, "Add R-group template"));
    await session.step(47, "And \"Markush Enumerator\" dialog should contain the text \"r group 9 · R4\"", () => shouldContainText(page, el("\"Markush Enumerator\" dialog"), "r group 9 · R4"));
    await session.step(48, "And \"Markush Enumerator\" dialog should contain the text \"6,561 molecules will be generated\"", () => shouldContainText(page, el("\"Markush Enumerator\" dialog"), "6,561 molecules will be generated"));
    await session.step(49, "And OK button in \"Markush Enumerator\" dialog should be enabled", () => shouldBe(page, el("OK button in \"Markush Enumerator\" dialog"), "enabled"));
    await session.step(50, "When user enters \"BddMarkushAlkene\" into \"Table name\" input in \"Markush Enumerator\" dialog", () => enterInto(page, "BddMarkushAlkene", el("\"Table name\" input in \"Markush Enumerator\" dialog")));
    await session.step(51, "And user clicks on OK button in \"Markush Enumerator\" dialog", () => clickOn(page, el("OK button in \"Markush Enumerator\" dialog")));
    await session.step(52, "Then the \"BddMarkushAlkene\" view should be current", () => viewIsCurrent(page, "BddMarkushAlkene"));
    await session.step(53, "And the table should have 1701 rows", () => rowCount(page, 1701));
    await session.step(54, "And \"Enumerated\" column should have no missing values", () => columnComplete(page, "Enumerated"));
    await session.step(55, "And the molecule in row 1 of \"Enumerated\" column should be \"C=C\"", () => rowMolecule(page, 1, "Enumerated", "C=C"));
    await session.step(56, "And the molecule in row 47 of \"Enumerated\" column should be \"C(/C)=C/CC\"", () => rowMolecule(page, 47, "Enumerated", "C(/C)=C/CC"));
    await session.step(57, "And the molecule in row 63 of \"Enumerated\" column should be \"C/C=C/CC\"", () => rowMolecule(page, 63, "Enumerated", "C/C=C/CC"));
    await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
});
