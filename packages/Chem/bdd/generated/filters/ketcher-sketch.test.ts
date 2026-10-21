/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/ketcher-sketch.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters]
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
import {checkFilterAsYouDraw, uncheckFilterAsYouDraw} from '../../bindings/dialogs.js';
import {placeBenzene} from '../../bindings/ketcher.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs, userSettingsPutBack} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A structure drawn in Ketcher reaches the substructure card as drawn", () => {
  const session = feature(test, "features/filters/ketcher-sketch.feature", import.meta.url);
  test("A structure drawn in Ketcher reaches the substructure card as drawn", {tag: ["@journey", "@sketcher-controls", "@realizes:filters.cp.chem-and-bio-filters", "@serial"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the \"chem-filter\" user settings are put back at feature end", () => userSettingsPutBack(page, "chem-filter"));
    await session.step(19, "And the molecule sketcher is \"Ketcher\"", () => sketcherIs(page, "Ketcher"));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("With Filter as you draw cleared, OK right after the stroke filters by it", async () => {
      await session.step(24, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(25, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(26, "When user unchecks \"Filter as you draw\" input in sketcher dialog", () => uncheckFilterAsYouDraw(page));
      await session.step(27, "And user places the benzene template on the Ketcher canvas", () => placeBenzene(page));
      await session.step(28, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(29, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(30, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(31, "And \"Sketch\" text in \"canonical_smiles\" filter card should not be visible", () => shouldNotBe(page, el("\"Sketch\" text in \"canonical_smiles\" filter card"), "visible"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw checked, the rows follow the stroke while the pointer rests on the canvas", async () => {
      await session.step(35, "When user hovers over \"canonical_smiles\" filter card", () => hoverOver(page, el("\"canonical_smiles\" filter card")));
      await session.step(36, "And user clicks on close of \"canonical_smiles\" filter card", () => clickOn(page, el("close of \"canonical_smiles\" filter card")));
      await session.step(37, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(38, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(39, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(40, "When user checks \"Filter as you draw\" input in sketcher dialog", () => checkFilterAsYouDraw(page));
      await session.step(41, "And user places the benzene template on the Ketcher canvas", () => placeBenzene(page));
      await session.step(42, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(43, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(44, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(45, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(46, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
