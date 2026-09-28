/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/sketcher/cell-editor.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.sketcher-cell-editor]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clipboardMolecule} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, noErrors, pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The sketcher opened from a molecule cell", () => {
  const session = feature(test, "features/sketcher/cell-editor.feature", import.meta.url);
  test("The sketcher opened from a molecule cell", {tag: ["@journey", "@realizes:chem.cp.sketcher-cell-editor"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens smiles-50 dataset", () => openDataset(page, ds("smiles-50")));
    await run.scenario("A double-click opens the sketcher on the cell's molecule", async () => {
      await session.step(16, "When user double-clicks on the \"cell 1 of canonical_smiles\" area of grid", () => doubleClickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(17, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(18, "And molecule input of sketcher dialog should be visible", () => shouldBe(page, el("molecule input of sketcher dialog"), "visible"));
      await session.step(19, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
      await session.step(20, "And user picks \"Copy as SMILES\" from the open menu", () => pickFromOpenMenu(page, "Copy as SMILES"));
      await session.step(21, "Then the clipboard should hold the molecule of row 1 of \"canonical_smiles\" column", () => clipboardMolecule(page, 1, "canonical_smiles"));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A molecule typed into the sketcher replaces the one in the cell", async () => {
      await session.step(25, "When user types \"C1CCCCC1\" into molecule input of sketcher dialog", () => typeInto(page, "C1CCCCC1", el("molecule input of sketcher dialog")));
      await session.step(26, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(27, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(28, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
      await session.step(29, "And the value of \"canonical_smiles\" column in row 1 should be \"C1CCCCC1\"", () => valueInRow(page, "canonical_smiles", 1, "C1CCCCC1"));
      await session.step(30, "And the table should have 50 rows", () => rowCount(page, 50));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Reopening the cell brings the sketcher up on the new molecule", async () => {
      await session.step(34, "When user double-clicks on the \"cell 1 of canonical_smiles\" area of grid", () => doubleClickArea(page, "cell 1 of canonical_smiles", el("grid")));
      await session.step(35, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(36, "When user clicks on \"Options\" icon in sketcher dialog", () => clickOn(page, el("\"Options\" icon in sketcher dialog")));
      await session.step(37, "And user picks \"Copy as SMILES\" from the open menu", () => pickFromOpenMenu(page, "Copy as SMILES"));
      await session.step(38, "Then the clipboard should hold the molecule of row 1 of \"canonical_smiles\" column", () => clipboardMolecule(page, 1, "canonical_smiles"));
      await session.step(39, "When user clicks on CANCEL button in sketcher dialog", () => clickOn(page, el("CANCEL button in sketcher dialog")));
      await session.step(40, "Then sketcher dialog should be absent", () => shouldBe(page, el("sketcher dialog"), "absent"));
      await session.step(41, "And the value of \"canonical_smiles\" column in row 1 should be \"C1CCCCC1\"", () => valueInRow(page, "canonical_smiles", 1, "C1CCCCC1"));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
