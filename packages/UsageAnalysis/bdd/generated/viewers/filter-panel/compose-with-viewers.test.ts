/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/compose-with-viewers.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, filterPassesFewer, noneOfFiltered, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, enterIntoArea, noErrors, pickFromContextMenu, propertyShouldBe, readingAtLeast, readingIs, readingReads, resizeTo, setProperty, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The panel's criterion composes with the viewers", () => {
  const session = feature(test, "features/viewers/filter-panel/compose-with-viewers.feature", import.meta.url);
  test("The panel's criterion composes with the viewers", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(14, "And user adds a categorical filter on \"RACE\" keeping \"Caucasian\"", () => addCategoricalFilter(page, "RACE", "Caucasian"));
    await session.step(15, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
    await session.step(16, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
    await session.step(17, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
    await run.scenario("A scatter plot zoom narrows the rows the card left and Reset View gives them back", async () => {
      await session.step(20, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","HEIGHT"]]));
      await session.step(23, "Then \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(24, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(25, "Then fewer than 896 rows should pass the filter", () => filterPassesFewer(page, 896));
      await session.step(26, "And no rows where \"RACE\" is \"Black\" should pass the filter", () => noneOfFiltered(page, "RACE", "Black"));
      await session.step(27, "And the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(28, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
      await session.step(29, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(30, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(31, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(32, "When user clicks on close icon of scatter plot viewer", () => clickOn(page, el("close icon of scatter plot viewer")));
      await session.step(33, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A bar click keeps the bar's category inside the card's", async () => {
      await session.step(37, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","SEX"]]));
      await session.step(39, "And user sets \"On Click\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("bar chart viewer"), "Filter"));
      await session.step(40, "And user clicks on the \"bar M\" area of bar chart viewer", () => clickArea(page, "bar M", el("bar chart viewer")));
      await session.step(41, "Then 416 rows should pass the filter", () => filterPasses(page, 416));
      await session.step(42, "And the \"rows shown\" reading of grid should be 416", () => readingIs(page, "rows shown", el("grid"), 416));
      await session.step(43, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
      await session.step(44, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(45, "When user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(46, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A histogram range narrows the card's rows further", async () => {
      await session.step(50, "When user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]));
      await session.step(54, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
      await session.step(55, "And user enters \"30\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "30", "range min input", el("histogram viewer")));
      await session.step(56, "Then 771 rows should pass the filter", () => filterPasses(page, 771));
      await session.step(57, "When user enters \"60\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "60", "range max input", el("histogram viewer")));
      await session.step(58, "Then 633 rows should pass the filter", () => filterPasses(page, 633));
      await session.step(59, "And the \"selected categories of RACE\" reading of filter panel should be \"Caucasian\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Caucasian"));
      await session.step(60, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(61, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(62, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The grid shows exactly the rows the card keeps", async () => {
      await session.step(66, "Then the \"rows shown\" reading of grid should be 896", () => readingIs(page, "rows shown", el("grid"), 896));
      await session.step(67, "When user clicks on the \"category Asian of RACE\" area of filter panel", () => clickArea(page, "category Asian of RACE", el("filter panel")));
      await session.step(68, "Then 15 rows should pass the filter", () => filterPasses(page, 15));
      await session.step(69, "And the \"rows shown\" reading of grid should be 15", () => readingIs(page, "rows shown", el("grid"), 15));
      await session.step(70, "When user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
      await session.step(71, "Then 896 rows should pass the filter", () => filterPasses(page, 896));
      await session.step(72, "And the \"rows shown\" reading of grid should be 896", () => readingIs(page, "rows shown", el("grid"), 896));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
