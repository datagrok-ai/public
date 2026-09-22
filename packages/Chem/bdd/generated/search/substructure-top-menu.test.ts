/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/search/substructure-top-menu.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [chem.cp.substructure-search-top-menu]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {filterPassesMatching} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll, rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Substructure search from the Search menu", () => {
  const session = feature(test, "features/search/substructure-top-menu.feature", import.meta.url);
  test("Substructure search from the Search menu", {tag: ["@journey", "@realizes:chem.cp.substructure-search-top-menu"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("The command prepares an empty card and opens its sketcher", async () => {
      await session.step(17, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(18, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(19, "And the \"cards\" reading of filter panel should contain \"canonical_smiles\"", () => readingContains(page, "cards", el("filter panel"), "canonical_smiles"));
      await session.step(20, "And the \"structure of canonical_smiles\" reading of filter panel should be \"\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), ""));
      await session.step(21, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(22, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Benzene keeps the molecules that contain it", async () => {
      await session.step(25, "When user types \"c1ccccc1\" into molecule input of sketcher dialog", () => typeInto(page, "c1ccccc1", el("molecule input of sketcher dialog")));
      await session.step(26, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(27, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(28, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(29, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(30, "And the filter should pass exactly the molecules of \"canonical_smiles\" column containing \"c1ccccc1\"", () => filterPassesMatching(page, "canonical_smiles", "c1ccccc1"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A gold atom keeps no molecule and deletes no row", async () => {
      await session.step(34, "When user clicks on the \"card canonical_smiles\" area of filter panel", () => clickArea(page, "card canonical_smiles", el("filter panel")));
      await session.step(35, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(36, "When user types \"[Au]\" into molecule input of sketcher dialog", () => typeInto(page, "[Au]", el("molecule input of sketcher dialog")));
      await session.step(37, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(38, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(39, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(40, "And the filter should pass exactly the molecules of \"canonical_smiles\" column containing \"[Au]\"", () => filterPassesMatching(page, "canonical_smiles", "[Au]"));
      await session.step(41, "And the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(42, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second invocation asks for the column and starts from an empty card", async () => {
      await session.step(45, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(46, "Then \"Substructure search\" dialog should be visible", () => shouldBe(page, el("\"Substructure search\" dialog"), "visible"));
      await session.step(47, "When user clicks on OK button in \"Substructure search\" dialog", () => clickOn(page, el("OK button in \"Substructure search\" dialog")));
      await session.step(48, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(49, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(50, "When user types \"C(=O)O\" into molecule input of sketcher dialog", () => typeInto(page, "C(=O)O", el("molecule input of sketcher dialog")));
      await session.step(51, "And user presses Enter in molecule input of sketcher dialog", () => pressKeyIn(page, "Enter", el("molecule input of sketcher dialog")));
      await session.step(52, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(53, "Then 314 rows should pass the filter", () => filterPasses(page, 314));
      await session.step(54, "And the filter should pass exactly the molecules of \"canonical_smiles\" column containing \"C(=O)O\"", () => filterPassesMatching(page, "canonical_smiles", "C(=O)O"));
      await session.step(55, "And the \"filters\" reading of filter panel should be 1", () => readingIs(page, "filters", el("filter panel"), 1));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
