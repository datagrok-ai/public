/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/clustering/sort-and-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [GROK-13041]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noErrors, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, readingsEqual, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The tree next to the grid under a filter and under a sort", () => {
  const session = feature(test, "features/clustering/sort-and-filter.feature", import.meta.url);
  test("The tree next to the grid under a filter and under a sort", {tag: ["@journey", "@realizes:GROK-13041"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens mol1K dataset", () => openDataset(page, ds("mol1K")));
    await session.step(11, "And user watches the task bar", () => watchTaskBar(page));
    await session.step(12, "When user picks \"Chem > Analyze > Hierarchical Clustering...\" from the top menu", () => pickFromTopMenu(page, "Chem > Analyze > Hierarchical Clustering..."));
    await session.step(13, "And user clicks on OK button in \"Hierarchical Clustering\" dialog", () => clickOn(page, el("OK button in \"Hierarchical Clustering\" dialog")));
    await session.step(14, "Then the task bar should have finished \"Creating dendrogram\"", () => taskBarFinished(page, "Creating dendrogram"));
    await session.step(15, "And the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
    await run.scenario("A filter keeps the tree in step with the grid and shows no overlay", async () => {
      await session.step(18, "When user remembers the \"tree leaves\" reading of grid", () => rememberReading(page, "tree leaves", el("grid")));
      await session.step(19, "And user clicks on first \"Toggle filters\" icon", () => clickOn(page, el("first \"Toggle filters\" icon")));
      await session.step(20, "And user clicks on the \"category Active_Integrase of Activity_Integrase\" area of filter panel", () => clickArea(page, "category Active_Integrase of Activity_Integrase", el("filter panel")));
      await session.step(21, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(22, "And the \"tree leaves\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "tree leaves", el("grid")));
      await session.step(23, "And the \"tree leaves\" and \"rows shown\" readings of grid should be the same", () => readingsEqual(page, "tree leaves", "rows shown", el("grid")));
      await session.step(24, "And grid should not contain text \"Revert columns sort order to see Dendrogram Tree\"", () => shouldNotContainText(page, el("grid"), "Revert columns sort order to see Dendrogram Tree"));
      await session.step(25, "And \"Revert sort\" button should be absent", () => shouldBe(page, el("\"Revert sort\" button"), "absent"));
      await session.step(26, "And \"Assign Clusters\" icon should be visible", () => shouldBe(page, el("\"Assign Clusters\" icon"), "visible"));
      await session.step(27, "When user clicks on the \"checkbox Inactive_Integrase of Activity_Integrase\" area of filter panel", () => clickArea(page, "checkbox Inactive_Integrase of Activity_Integrase", el("filter panel")));
      await session.step(28, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(29, "And the \"tree leaves\" reading of grid should be 1000", () => readingIs(page, "tree leaves", el("grid"), 1000));
      await session.step(30, "And grid should not contain text \"Revert columns sort order to see Dendrogram Tree\"", () => shouldNotContainText(page, el("grid"), "Revert columns sort order to see Dendrogram Tree"));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A sort shows the overlay, and Revert sort puts the rows back in the tree's order", async () => {
      await session.step(34, "When user remembers the \"row order\" reading of grid", () => rememberReading(page, "row order", el("grid")));
      await session.step(35, "And user double-clicks on the \"header pIC50_HIV_Integrase\" area of grid", () => doubleClickArea(page, "header pIC50_HIV_Integrase", el("grid")));
      await session.step(36, "Then the \"sort column\" reading of grid should be \"pIC50_HIV_Integrase\"", () => readingReads(page, "sort column", el("grid"), "pIC50_HIV_Integrase"));
      await session.step(37, "And the \"row order\" reading of grid should not be as remembered", () => readingNotAsRemembered(page, "row order", el("grid")));
      await session.step(38, "And \"Revert sort\" button should be visible", () => shouldBe(page, el("\"Revert sort\" button"), "visible"));
      await session.step(39, "And grid should contain text \"Revert columns sort order to see Dendrogram Tree\"", () => shouldContainText(page, el("grid"), "Revert columns sort order to see Dendrogram Tree"));
      await session.step(40, "When user clicks on \"Revert sort\" button", () => clickOn(page, el("\"Revert sort\" button")));
      await session.step(41, "Then the \"row order\" reading of grid should be as remembered", () => readingAsRemembered(page, "row order", el("grid")));
      await session.step(42, "And \"Revert sort\" button should be absent", () => shouldBe(page, el("\"Revert sort\" button"), "absent"));
      await session.step(43, "And grid should not contain text \"Revert columns sort order to see Dendrogram Tree\"", () => shouldNotContainText(page, el("grid"), "Revert columns sort order to see Dendrogram Tree"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
