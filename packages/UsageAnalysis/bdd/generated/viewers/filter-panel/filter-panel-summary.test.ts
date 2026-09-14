/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-summary.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {addCardFor} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, shouldBe, shouldContainText, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addRangeFilter, filterPasses, filterPassesAll, filterPassesFewer, noneOfFiltered, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragZoomOverArea, noErrors, readingAtLeast, readingIs, readingLowerThanRemembered, readingReads, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel indicator with viewers filtering", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-summary.feature", import.meta.url);
  test("Filter panel indicator with viewers filtering", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(20, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(21, "And user clicks on the \"category Caucasian of RACE\" area of filter panel", () => clickArea(page, "category Caucasian of RACE", el("filter panel")));
    await session.step(22, "And user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
    await session.step(23, "Then 633 rows should pass the filter", () => filterPasses(page, 633));
    await session.step(24, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
    await run.scenario("Viewers filtering on top of the cards leave the counter at the cards", async () => {
      await session.step(27, "When user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","HEIGHT"]]));
      await session.step(30, "And user drags a zoom box over the \"view\" area of scatter plot viewer", () => dragZoomOverArea(page, "view", el("scatter plot viewer")));
      await session.step(31, "Then fewer than 633 rows should pass the filter", () => filterPassesFewer(page, 633));
      await session.step(32, "And counter of filter panel should be visible", () => shouldBe(page, el("counter of filter panel"), "visible"));
      await session.step(33, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(34, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","SEX"]]));
      await session.step(36, "And user sets \"On Click\" property of bar chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("bar chart viewer"), "Filter"));
      await session.step(37, "And user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(38, "And user clicks on the \"bar M\" area of bar chart viewer", () => clickArea(page, "bar M", el("bar chart viewer")));
      await session.step(39, "Then no rows where \"SEX\" is \"F\" should pass the filter", () => noneOfFiltered(page, "SEX", "F"));
      await session.step(40, "And the \"rows shown\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(41, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(42, "When user adds a pie chart viewer with:", () => addViewerWith(page, "pie chart", [["Category","DIS_POP"]]));
      await session.step(44, "And user sets \"On Click\" property of pie chart viewer to \"Filter\"", () => setProperty(page, "On Click", el("pie chart viewer"), "Filter"));
      await session.step(45, "And user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(46, "And user clicks on the \"slice RA\" area of pie chart viewer", () => clickArea(page, "slice RA", el("pie chart viewer")));
      await session.step(47, "Then no rows where \"DIS_POP\" is \"UC\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "UC"));
      await session.step(48, "And the \"rows shown\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(49, "And the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(50, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(51, "When user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]));
      await session.step(53, "And user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(54, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 450 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 450));
      await session.step(55, "Then the \"filtering\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "filtering", el("pc plot viewer"), "true"));
      await session.step(56, "And the \"rows shown\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(57, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(58, "When user hovers over counter of filter panel", () => hoverOver(page, el("counter of filter panel")));
      await session.step(59, "Then tooltip should contain the text \"[30,60]\"", () => shouldContainText(page, el("tooltip"), "[30,60]"));
      await session.step(60, "And tooltip should contain the text \"RACE\"", () => shouldContainText(page, el("tooltip"), "RACE"));
      await session.step(61, "And tooltip should contain the text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(62, "And tooltip should contain the text \"AGE\"", () => shouldContainText(page, el("tooltip"), "AGE"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A trellis cell filters on top of the cards and the counter stays at the cards", async () => {
      await session.step(66, "When user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEVERITY"],["Y Column Names","SEX"],["Viewer Type","Scatter plot"]]));
      await session.step(70, "And user sets \"On Click\" property of trellis plot viewer to \"Filter\"", () => setProperty(page, "On Click", el("trellis plot viewer"), "Filter"));
      await session.step(71, "And user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(72, "And user clicks on the \"cell None | M\" area of trellis plot viewer", () => clickArea(page, "cell None | M", el("trellis plot viewer")));
      await session.step(73, "Then no rows where \"SEVERITY\" is \"High\" should pass the filter", () => noneOfFiltered(page, "SEVERITY", "High"));
      await session.step(74, "And the \"rows shown\" reading of filter panel should be lower than remembered", () => readingLowerThanRemembered(page, "rows shown", el("filter panel")));
      await session.step(75, "And the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(76, "And counter of filter panel should have text \"2\"", () => shouldHaveText(page, el("counter of filter panel"), "2"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The panel's summary icon names the viewers' filtering as well as the cards", async () => {
      await session.step(81, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(82, "Then question-circle icon in filter panel should be visible", () => shouldBe(page, el("question-circle icon in filter panel"), "visible"));
      await session.step(83, "When user hovers over question-circle icon in filter panel", () => hoverOver(page, el("question-circle icon in filter panel")));
      await session.step(84, "Then tooltip should contain the text \"RACE\"", () => shouldContainText(page, el("tooltip"), "RACE"));
      await session.step(85, "And tooltip should contain the text \"Caucasian\"", () => shouldContainText(page, el("tooltip"), "Caucasian"));
      await session.step(86, "And tooltip should contain the text \"[30,60]\"", () => shouldContainText(page, el("tooltip"), "[30,60]"));
      await session.step(87, "And tooltip should contain the text \"Scatter plot\"", () => shouldContainText(page, el("tooltip"), "Scatter plot"));
      await session.step(88, "And tooltip should contain the text \"Bar chart\"", () => shouldContainText(page, el("tooltip"), "Bar chart"));
      await session.step(89, "And tooltip should contain the text \"Pie chart\"", () => shouldContainText(page, el("tooltip"), "Pie chart"));
      await session.step(90, "And tooltip should contain the text \"Trellis plot\"", () => shouldContainText(page, el("tooltip"), "Trellis plot"));
      await session.step(91, "And tooltip should contain the text \"PC Plot\"", () => shouldContainText(page, el("tooltip"), "PC Plot"));
    }, {knownFailure: true});
    await run.scenario("Reset filters empties the counter and gives back every row", async () => {
      await session.step(94, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(95, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(96, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(97, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(98, "And the \"filters\" reading of filter panel should be 0", () => readingIs(page, "filters", el("filter panel"), 0));
      await session.step(99, "And the \"cards\" reading of filter panel should be \"AGE, RACE\"", () => readingReads(page, "cards", el("filter panel"), "AGE, RACE"));
      await session.step(100, "And the \"filtering\" reading of pc plot viewer should be \"false\"", () => readingReads(page, "filtering", el("pc plot viewer"), "false"));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
