/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/filters/column-popup.feature
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
import {filterPassesMatching} from '../../bindings/molecules.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, pressKey, pressKeyIn, shouldBe, shouldHaveText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, filterPassesAll, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, openDataset, sketcherIs} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A column popup's substructure filter moved to the filter panel", () => {
  const session = feature(test, "features/filters/column-popup.feature", import.meta.url);
  test("A column popup's substructure filter moved to the filter panel", {tag: ["@journey", "@realizes:filters.cp.chem-and-bio-filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the molecule sketcher is \"OpenChemLib\"", () => sketcherIs(page, "OpenChemLib"));
    await session.step(13, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(14, "And user opens spgi-100 dataset", () => openDataset(page, ds("spgi-100")));
    await run.scenario("The popup's filter narrows the rows on its own", async () => {
      await session.step(17, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(18, "And user hovers over the \"header Structure\" area of grid", () => hoverArea(page, "header Structure", el("grid")));
      await session.step(19, "And user clicks on \"Column options\" icon in grid", () => clickOn(page, el("\"Column options\" icon in grid")));
      await session.step(20, "Then column popup should be visible", () => shouldBe(page, el("column popup"), "visible"));
      await session.step(21, "And title of column popup should have text \"Structure\"", () => shouldHaveText(page, el("title of column popup"), "Structure"));
      await session.step(22, "And \"Structure\" filter card should be absent", () => shouldBe(page, el("\"Structure\" filter card"), "absent"));
      await session.step(23, "When user types \"c1ccncc1\" into molecule input of column popup", () => typeInto(page, "c1ccncc1", el("molecule input of column popup")));
      await session.step(24, "And user presses Enter in molecule input of column popup", () => pressKeyIn(page, "Enter", el("molecule input of column popup")));
      await session.step(25, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(26, "And the filter should pass exactly the molecules of \"Structure\" column containing \"c1ccncc1\"", () => filterPassesMatching(page, "Structure", "c1ccncc1"));
      await session.step(27, "And \"Structure\" filter card should be absent", () => shouldBe(page, el("\"Structure\" filter card"), "absent"));
      await session.step(28, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Add filter moves the criterion to the panel", async () => {
      await session.step(31, "When user clicks on \"Add filter\" action in column popup", () => clickOn(page, el("\"Add filter\" action in column popup")));
      await session.step(32, "Then \"Structure\" filter card should be visible", () => shouldBe(page, el("\"Structure\" filter card"), "visible"));
      await session.step(33, "And molecule input of column popup should be absent", () => shouldBe(page, el("molecule input of column popup"), "absent"));
      await session.step(34, "And the \"structure of Structure\" reading of filter panel should be \"c1ccncc1\"", () => readingReads(page, "structure of Structure", el("filter panel"), "c1ccncc1"));
      await session.step(35, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(36, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(37, "Then column popup should be absent", () => shouldBe(page, el("column popup"), "absent"));
      await session.step(38, "And 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(39, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The panel card switched off lets every row through", async () => {
      await session.step(42, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(43, "And user unchecks checkbox of \"Structure\" filter card", () => uncheck(page, el("checkbox of \"Structure\" filter card")));
      await session.step(44, "Then \"Structure\" filter card should be disabled", () => shouldBe(page, el("\"Structure\" filter card"), "disabled"));
      await session.step(45, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(46, "When user hovers over \"Structure\" filter card", () => hoverOver(page, el("\"Structure\" filter card")));
      await session.step(47, "And user clicks on checkbox of \"Structure\" filter card", () => clickOn(page, el("checkbox of \"Structure\" filter card")));
      await session.step(48, "Then 17 rows should pass the filter", () => filterPasses(page, 17));
      await session.step(49, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
