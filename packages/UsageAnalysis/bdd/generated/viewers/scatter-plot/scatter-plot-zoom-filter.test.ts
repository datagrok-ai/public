/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-zoom-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.scatter-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, filterBetween, filterIsExactly, filterPasses, filterPassesAll, filterPassesFewer, openFilterPanel, removeColumn, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, pickFromContextMenu, propertyShouldBe, readingAsRemembered, readingLower, readingSame, rememberRange, rememberReading, rememberedRange, setProperties, setProperty, showsFewerRows, showsRows, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot zoom and filter synchronization", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-zoom-filter.feature", import.meta.url);
  test("Scatter plot zoom and filter synchronization", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]));
    await session.step(17, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(18, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A zoom filters the table and Reset View gives it back", async () => {
      await session.step(21, "Then \"Zoom and Filter\" property of scatter plot viewer should be \"filter by zoom\"", () => propertyShouldBe(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(22, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(23, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(24, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(25, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(26, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(27, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(28, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Pack and zoom by filter does not leave the table filtered", async () => {
      await session.step(32, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(33, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(34, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"pack and zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "pack and zoom by filter"));
      await session.step(35, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(36, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(37, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(38, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(39, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(40, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Under zoom by filter an external filter narrows the viewport", async () => {
      await session.step(44, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"zoom by filter\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "zoom by filter"));
      await session.step(45, "And user remembers the value range of scatter plot viewer", () => rememberRange(page, el("scatter plot viewer")));
      await session.step(46, "And user remembers the \"x axis span\" reading of scatter plot viewer", () => rememberReading(page, "x axis span", el("scatter plot viewer")));
      await session.step(47, "And user filters rows where \"WEIGHT\" is between 90.96 and 115.64", () => filterBetween(page, "WEIGHT", 90.96, 115.64));
      await session.step(48, "Then 212 rows should pass the filter", () => filterPasses(page, 212));
      await session.step(49, "And the \"x axis span\" reading of scatter plot viewer should be lower than before", () => readingLower(page, "x axis span", el("scatter plot viewer")));
      await session.step(50, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(51, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(52, "When user resets the filter", () => resetFilter(page));
      await session.step(53, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(54, "And the \"x axis span\" reading of scatter plot viewer should be as remembered", () => readingAsRemembered(page, "x axis span", el("scatter plot viewer")));
      await session.step(55, "And scatter plot viewer should show the remembered value range", () => rememberedRange(page, el("scatter plot viewer")));
      await session.step(56, "When user sets \"Zoom and Filter\" property of scatter plot viewer to \"filter by zoom\"", () => setProperty(page, "Zoom and Filter", el("scatter plot viewer"), "filter by zoom"));
      await session.step(57, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The filter panel's reset undoes what the zoom filtered", async () => {
      await session.step(60, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(61, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(62, "When user scrolls the mouse wheel up over the \"view\" area of scatter plot viewer", () => wheelOverArea(page, "up", "view", el("scatter plot viewer")));
      await session.step(63, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
      await session.step(64, "And scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(65, "And the \"rows selected\" reading of scatter plot viewer should be the same as before", () => readingSame(page, "rows selected", el("scatter plot viewer")));
      await session.step(66, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(67, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(68, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(69, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(70, "When user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(71, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A logarithmic axis alone filters nothing; Filter Out Invalid drops the non-positive rows", async () => {
      await session.step(75, "When user adds a calculated column \"AGE_SHIFT\" with formula \"${AGE} - 40\"", () => addCalculated(page, "AGE_SHIFT", "${AGE} - 40"));
      await session.step(76, "And user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["X","WEIGHT"],["Y","AGE_SHIFT"]]));
      await session.step(79, "Then scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(80, "And \"Filter Out Invalid\" property of scatter plot viewer should be \"false\"", () => propertyShouldBe(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(81, "When user sets \"Y Axis Type\" property of scatter plot viewer to \"logarithmic\"", () => setProperty(page, "Y Axis Type", el("scatter plot viewer"), "logarithmic"));
      await session.step(82, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(83, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"true\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "true"));
      await session.step(84, "Then 635 rows should pass the filter", () => filterPasses(page, 635));
      await session.step(85, "And the filter should pass exactly the rows where \"AGE\" is between 40.5 and 100", () => filterIsExactly(page, "AGE", 40.5, 100));
      await session.step(86, "When user sets \"Filter Out Invalid\" property of scatter plot viewer to \"false\"", () => setProperty(page, "Filter Out Invalid", el("scatter plot viewer"), "false"));
      await session.step(87, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(88, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Y Axis Type","linear"],["Filter Out Invalid","true"]]));
      await session.step(91, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(92, "And scatter plot viewer should show 1000 rows", () => showsRows(page, el("scatter plot viewer"), 1000));
      await session.step(93, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Filter Out Invalid","false"],["X","WEIGHT"],["Y","HEIGHT"]]));
      await session.step(97, "And user removes \"AGE_SHIFT\" column", () => removeColumn(page, "AGE_SHIFT"));
      await session.step(98, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")));
      await session.step(99, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(100, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
