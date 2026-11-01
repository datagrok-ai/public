/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, public#790, public#791, public#3059, public#1464, public#319, public#2501, public#2943, public#1242]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, hoverOver, keepsWidth, pressKey, rememberSize, shouldBe, shouldNotBe, sizeAgainstRemembered} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {filterPanelCount, filterPasses, filterPassesAll, openEmptyFilterPanel, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeCurrentView, openDataset, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, dragCard, openCardIndicatorMenu, pickCardIndicatorMenu, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewer, clickArea, closeContextMenu, dragAreaToArea, hoverArea, loadLayout, menuLists, noErrors, readingReads, saveLayoutToServer, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel regressions", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-regressions.feature", import.meta.url);
  test("Filter panel regressions", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@realizes:public#790", "@realizes:public#791", "@realizes:public#3059", "@realizes:public#1464", "@realizes:public#319", "@realizes:public#2501", "@realizes:public#2943", "@realizes:public#1242"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 6, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(24, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(25, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("A missing-values choice that came back with a layout can be undone, and Reset filters puts it back to Keep (public#790, public#791)", async () => {
      await session.step(29, "When user adds a card for \"HEIGHT\" to the filter panel", () => addCardFor(page, "HEIGHT"));
      await session.step(30, "And user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
      await session.step(31, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(32, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(33, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(34, "And user loads the saved layout", () => loadLayout(page));
      await session.step(35, "Then \"HEIGHT\" filter card should be visible", () => shouldBe(page, el("\"HEIGHT\" filter card"), "visible"));
      await session.step(36, "And 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(37, "When user picks \"Missing values | Keep missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Keep missing values", "HEIGHT"));
      await session.step(38, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(39, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(40, "When user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
      await session.step(41, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(42, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(43, "And the \"filtering of HEIGHT\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of HEIGHT", el("filter panel"), "true"));
      await session.step(44, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(45, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(46, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(47, "And the \"filtering of HEIGHT\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of HEIGHT", el("filter panel"), "false"));
      await session.step(48, "When user opens the indicator menu of the \"HEIGHT\" filter card", () => openCardIndicatorMenu(page, "HEIGHT"));
      await session.step(49, "Then the open menu should list \"Missing values > Keep missing values\"", () => menuLists(page, "Missing values > Keep missing values"));
      await session.step(50, "And \"Keep missing values\" menu item should be selected", () => shouldBe(page, el("\"Keep missing values\" menu item"), "selected"));
      await session.step(51, "And \"Filter out missing values\" menu item should not be selected", () => shouldNotBe(page, el("\"Filter out missing values\" menu item"), "selected"));
      await session.step(52, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(53, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(54, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved after a card search was reset brings back no filtering (public#3059)", async () => {
      await session.step(59, "When user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
      await session.step(60, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(61, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(62, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(63, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(64, "And the \"filtering of DIS_POP\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of DIS_POP", el("filter panel"), "true"));
      await session.step(65, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(66, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(67, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(68, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(69, "And user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
      await session.step(70, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(71, "When user loads the saved layout", () => loadLayout(page));
      await session.step(72, "Then \"DIS_POP\" filter card should be visible", () => shouldBe(page, el("\"DIS_POP\" filter card"), "visible"));
      await session.step(73, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(74, "And the \"filtering of DIS_POP\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of DIS_POP", el("filter panel"), "false"));
      await session.step(75, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(76, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(77, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(78, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The order the cards were dragged into comes back with the layout (public#1464)", async () => {
      await session.step(82, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(83, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
      await session.step(84, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(85, "Then the \"cards\" reading of filter panel should be \"AGE, SEX, RACE\"", () => readingReads(page, "cards", el("filter panel"), "AGE, SEX, RACE"));
      await session.step(86, "When user drags the \"RACE\" filter card above the \"SEX\" filter card", () => dragCard(page, "RACE", "SEX"));
      await session.step(87, "Then the \"cards\" reading of filter panel should be \"AGE, RACE, SEX\"", () => readingReads(page, "cards", el("filter panel"), "AGE, RACE, SEX"));
      await session.step(88, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(89, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(90, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(91, "When user loads the saved layout", () => loadLayout(page));
      await session.step(92, "Then the \"cards\" reading of filter panel should be \"AGE, RACE, SEX\"", () => readingReads(page, "cards", el("filter panel"), "AGE, RACE, SEX"));
      await session.step(93, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(94, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(95, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A view cloned after Reset filters comes up unfiltered (public#319)", async () => {
      await session.step(99, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(100, "And user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(101, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(102, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(103, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(104, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(105, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(106, "Then the \"demog-1000 copy\" view should be current", () => viewIsCurrent(page, "demog-1000 copy"));
      await session.step(107, "And filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(108, "And \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(109, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(110, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(111, "When user closes the current view", () => closeCurrentView(page));
      await session.step(112, "And user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(113, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(114, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Adding viewers does not shrink the filter panel (public#2501)", async () => {
      await session.step(119, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(120, "And user remembers the size of filter panel", () => rememberSize(page, el("filter panel")));
      await session.step(121, "And user remembers the size of grid", () => rememberSize(page, el("grid")));
      await session.step(122, "And user adds histogram viewer", () => addViewer(page, "histogram"));
      await session.step(123, "And user adds bar chart viewer", () => addViewer(page, "bar chart"));
      await session.step(124, "And user adds line chart viewer", () => addViewer(page, "line chart"));
      await session.step(125, "Then the open tableview should have 1 line chart viewer", () => viewerCount(page, 1, "line chart"));
      await session.step(126, "And grid should be narrower than remembered", () => sizeAgainstRemembered(page, el("grid"), "narrower"));
      await session.step(127, "And filter panel should keep its remembered width", () => keepsWidth(page, el("filter panel")));
      await session.step(128, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header popup still filters after the column's card left the panel (public#2943, public#1242)", async () => {
      await session.step(132, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(133, "And user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(134, "And user clicks on close of \"AGE\" filter card", () => clickOn(page, el("close of \"AGE\" filter card")));
      await session.step(135, "Then \"AGE\" filter card should be absent", () => shouldBe(page, el("\"AGE\" filter card"), "absent"));
      await session.step(136, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(137, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(138, "And user clicks on \"Column options\" icon in grid", () => clickOn(page, el("\"Column options\" icon in grid")));
      await session.step(139, "Then column popup should be visible", () => shouldBe(page, el("column popup"), "visible"));
      await session.step(140, "When user hovers over the \"view\" area of histogram viewer in column popup", () => hoverArea(page, "view", el("histogram viewer in column popup")));
      await session.step(141, "And user drags the \"range min handle\" area of histogram viewer in column popup to the \"bin 10\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in column popup"), "bin 10"));
      await session.step(142, "Then 326 rows should pass the filter", () => filterPasses(page, 326));
      await session.step(143, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(144, "And user resets the filter", () => resetFilter(page));
      await session.step(145, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(146, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
