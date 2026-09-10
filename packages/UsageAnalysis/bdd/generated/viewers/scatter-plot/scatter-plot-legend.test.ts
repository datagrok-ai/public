/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/scatter-plot/scatter-plot-legend.feature
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
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, filterPassesAll, openFilterPanel, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickLegendItem, legendLists, legendSide, lessInk, noErrors, propertyShouldBe, setProperties, setProperty, showsFewerRows, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Scatter plot legend", () => {
  const session = feature(test, "features/viewers/scatter-plot/scatter-plot-legend.feature", import.meta.url);
  test("Scatter plot legend", {tag: ["@journey", "@viewers", "@realizes:viewers.scatter-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]));
    await session.step(19, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(20, "And legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
    await run.scenario("Color and Markers build one legend of two sections", async () => {
      await session.step(23, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(24, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
      await session.step(25, "And the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(26, "And legend of scatter plot viewer should contain text \"Asian\"", () => shouldContainText(page, el("legend of scatter plot viewer"), "Asian"));
      await session.step(27, "And legend of scatter plot viewer should contain text \"Caucasian\"", () => shouldContainText(page, el("legend of scatter plot viewer"), "Caucasian"));
      await session.step(28, "When user sets \"Markers\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "SEX"));
      await session.step(29, "Then the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(30, "When user sets \"Color\" property of scatter plot viewer to \"AGE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "AGE"));
      await session.step(31, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(32, "When user sets \"Color\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Color", el("scatter plot viewer"), "RACE"));
      await session.step(33, "Then the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(34, "And legend of scatter plot viewer should contain text \"Asian\"", () => shouldContainText(page, el("legend of scatter plot viewer"), "Asian"));
      await session.step(35, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Clearing Markers leaves the color entries alone", async () => {
      await session.step(38, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","SEX"],["Markers","SEX"]]));
      await session.step(41, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(42, "When user sets \"Markers\" property of scatter plot viewer to \"\"", () => setProperty(page, "Markers", el("scatter plot viewer"), ""));
      await session.step(43, "Then \"Markers\" property of scatter plot viewer should be \"\"", () => propertyShouldBe(page, "Markers", el("scatter plot viewer"), ""));
      await session.step(44, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(45, "And legend of scatter plot viewer should contain text \"F\"", () => shouldContainText(page, el("legend of scatter plot viewer"), "F"));
      await session.step(46, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Color","RACE"],["Markers","SEX"]]));
      await session.step(49, "Then the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table filter drops the filtered-out categories and Markers on the color column adds no section", async () => {
      await session.step(53, "When user sets \"Markers\" property of scatter plot viewer to \"\"", () => setProperty(page, "Markers", el("scatter plot viewer"), ""));
      await session.step(54, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(55, "When user opens the filter panel", () => openFilterPanel(page));
      await session.step(56, "And user adds a categorical filter on \"RACE\" keeping \"Asian, Caucasian\"", () => addCategoricalFilter(page, "RACE", "Asian, Caucasian"));
      await session.step(57, "Then 911 rows should pass the filter", () => filterPasses(page, 911));
      await session.step(58, "And the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(59, "When user sets \"Markers\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "RACE"));
      await session.step(60, "Then the legend of scatter plot viewer should list 2 items", () => legendLists(page, el("scatter plot viewer"), 2));
      await session.step(61, "When user sets \"Markers\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "SEX"));
      await session.step(62, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(63, "When user resets the filter", () => resetFilter(page));
      await session.step(64, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(65, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(66, "When user sets \"Markers\" property of scatter plot viewer to \"RACE\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "RACE"));
      await session.step(67, "Then the legend of scatter plot viewer should list 4 items", () => legendLists(page, el("scatter plot viewer"), 4));
      await session.step(68, "When user sets \"Markers\" property of scatter plot viewer to \"SEX\"", () => setProperty(page, "Markers", el("scatter plot viewer"), "SEX"));
      await session.step(69, "Then the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(70, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A click on a legend entry hides that category on the canvas, not in the table", async () => {
      await session.step(73, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(74, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(75, "When user clicks on \"Asian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Asian", el("scatter plot viewer")));
      await session.step(76, "Then scatter plot viewer should show fewer rows than before", () => showsFewerRows(page, el("scatter plot viewer")));
      await session.step(77, "And scatter plot viewer should have less ink than before", () => lessInk(page, el("scatter plot viewer")));
      await session.step(78, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(79, "When user clicks on \"Asian\" item in the legend of scatter plot viewer", () => clickLegendItem(page, "Asian", el("scatter plot viewer")));
      await session.step(80, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(81, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(82, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Legend Visibility and Legend Position decide whether and where it shows", async () => {
      await session.step(85, "When user sets \"Legend Visibility\" property of scatter plot viewer to \"Never\"", () => setProperty(page, "Legend Visibility", el("scatter plot viewer"), "Never"));
      await session.step(86, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(87, "When user sets \"Legend Visibility\" property of scatter plot viewer to \"Always\"", () => setProperty(page, "Legend Visibility", el("scatter plot viewer"), "Always"));
      await session.step(88, "Then legend of scatter plot viewer should be visible", () => shouldBe(page, el("legend of scatter plot viewer"), "visible"));
      await session.step(89, "And the legend of scatter plot viewer should list 6 items", () => legendLists(page, el("scatter plot viewer"), 6));
      await session.step(90, "When user sets \"Legend Position\" property of scatter plot viewer to \"Left\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Left"));
      await session.step(91, "Then the legend of scatter plot viewer should be on the left", () => legendSide(page, el("scatter plot viewer"), "left"));
      await session.step(92, "When user sets \"Legend Position\" property of scatter plot viewer to \"Right\"", () => setProperty(page, "Legend Position", el("scatter plot viewer"), "Right"));
      await session.step(93, "Then the legend of scatter plot viewer should be on the right", () => legendSide(page, el("scatter plot viewer"), "right"));
      await session.step(94, "When user sets properties of scatter plot viewer:", () => setProperties(page, el("scatter plot viewer"), [["Legend Position","Auto"],["Legend Visibility","Auto"],["Color",""],["Markers",""]]));
      await session.step(99, "Then legend of scatter plot viewer should be hidden", () => shouldBe(page, el("legend of scatter plot viewer"), "hidden"));
      await session.step(100, "And scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
