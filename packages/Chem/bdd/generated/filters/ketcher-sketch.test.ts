/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/ketcher-sketch.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [filters.cp.chem-and-bio-filters]
--- */
import {test} from '@playwright/test';
import '../../bindings/datasets.js';
import '../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {placeBenzene} from '../../bindings/ketcher.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, hoverOver, shouldBe, shouldNotBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPasses, filterPassesAll} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A structure drawn in Ketcher reaches the substructure card as drawn", () => {
  const session = feature(test, "features/filters/ketcher-sketch.feature", import.meta.url);
  test("A structure drawn in Ketcher reaches the substructure card as drawn", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the molecule sketcher is \"Ketcher\"", () => sketcherIs(page, "Ketcher"));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(15, "And user opens smiles dataset", () => openDataset(page, ds("smiles")));
    await run.scenario("With Filter as you draw cleared, OK right after the stroke filters by it", async () => {
      await session.step(18, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(19, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(20, "When user unchecks \"Filter as you draw\" input in sketcher dialog", () => uncheck(page, el("\"Filter as you draw\" input in sketcher dialog")));
      await session.step(21, "And user places the benzene template on the Ketcher canvas", () => placeBenzene(page));
      await session.step(22, "And user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(23, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(24, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(25, "And \"Sketch\" text in \"canonical_smiles\" filter card should not be visible", () => shouldNotBe(page, el("\"Sketch\" text in \"canonical_smiles\" filter card"), "visible"));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("With Filter as you draw checked, the rows follow the stroke while the pointer rests on the canvas", async () => {
      await session.step(29, "When user hovers over \"canonical_smiles\" filter card", () => hoverOver(page, el("\"canonical_smiles\" filter card")));
      await session.step(30, "And user clicks on close of \"canonical_smiles\" filter card", () => clickOn(page, el("close of \"canonical_smiles\" filter card")));
      await session.step(31, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(32, "When user picks \"Chem > Search > Substructure Search...\" from the top menu", () => pickFromTopMenu(page, "Chem > Search > Substructure Search..."));
      await session.step(33, "Then sketcher dialog should be visible", () => shouldBe(page, el("sketcher dialog"), "visible"));
      await session.step(34, "When user checks \"Filter as you draw\" input in sketcher dialog", () => check(page, el("\"Filter as you draw\" input in sketcher dialog")));
      await session.step(35, "And user places the benzene template on the Ketcher canvas", () => placeBenzene(page));
      await session.step(36, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(37, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(38, "When user clicks on OK button in sketcher dialog", () => clickOn(page, el("OK button in sketcher dialog")));
      await session.step(39, "Then 924 rows should pass the filter", () => filterPasses(page, 924));
      await session.step(40, "And the \"structure of canonical_smiles\" reading of filter panel should be \"c1ccccc1\"", () => readingReads(page, "structure of canonical_smiles", el("filter panel"), "c1ccccc1"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
