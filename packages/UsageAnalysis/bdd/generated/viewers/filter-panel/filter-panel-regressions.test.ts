/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters, public#790, public#2364, public#1343, public#2359, public#1986, public#3027, core#479, public#2501, public#2943, public#1242]
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
import {clickOn, clipboardHas, hoverOver, keepsWidth, pressKey, pressKeyIn, rememberSize, shouldBe, shouldHaveText, shouldNotBe, typeAtCaret} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCalculated, filterIsExactlyCategory, filterPanelCount, filterPanelHas, filterPasses, filterPassesAll, openEmptyFilterPanel, resetFilter} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, enterCardBound, openCardIndicatorMenu, pickCardIndicatorMenu, pickPanelMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {addViewer, closeContextMenu, dragAreaBy, dragAreaToArea, hasArea, hoverArea, menuLists, noErrors, readingIs, readingReads, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLies, dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel regressions", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-regressions.feature", import.meta.url);
  test("Filter panel regressions", {tag: ["@journey", "@viewers", "@realizes:viewers.filters", "@realizes:public#790", "@realizes:public#2364", "@realizes:public#1343", "@realizes:public#2359", "@realizes:public#1986", "@realizes:public#3027", "@realizes:core#479", "@realizes:public#2501", "@realizes:public#2943", "@realizes:public#1242"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(28, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(29, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
    await session.step(30, "And all rows should pass the filter", () => filterPassesAll(page));
    await run.scenario("Reset filters gives back the rows a missing-values choice dropped (public#790)", async () => {
      await session.step(34, "When user adds a card for \"HEIGHT\" to the filter panel", () => addCardFor(page, "HEIGHT"));
      await session.step(35, "And user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
      await session.step(36, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(37, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(38, "And the \"filtering of HEIGHT\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of HEIGHT", el("filter panel"), "true"));
      await session.step(39, "When user hovers over filter panel", () => hoverOver(page, el("filter panel")));
      await session.step(40, "And user clicks on reset icon of filter panel", () => clickOn(page, el("reset icon of filter panel")));
      await session.step(41, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(42, "And the \"filtering of HEIGHT\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of HEIGHT", el("filter panel"), "false"));
      await session.step(43, "When user opens the indicator menu of the \"HEIGHT\" filter card", () => openCardIndicatorMenu(page, "HEIGHT"));
      await session.step(44, "Then the open menu should list \"Missing values > Keep missing values\"", () => menuLists(page, "Missing values > Keep missing values"));
      await session.step(45, "And \"Keep missing values\" menu item should be selected", () => shouldBe(page, el("\"Keep missing values\" menu item"), "selected"));
      await session.step(46, "And \"Filter out missing values\" menu item should not be selected", () => shouldNotBe(page, el("\"Filter out missing values\" menu item"), "selected"));
      await session.step(47, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(48, "And user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(49, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(50, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A min/max range wholly outside the column's values keeps no row (public#2364)", async () => {
      await session.step(54, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(55, "And user picks \"Min / max\" from the indicator menu of the \"AGE\" filter card", () => pickCardIndicatorMenu(page, "Min / max", "AGE"));
      await session.step(56, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(57, "And user enters \"88.2\" into the min field of the \"AGE\" filter card", () => enterCardBound(page, "88.2", "min", "AGE"));
      await session.step(58, "And user enters \"88.8\" into the max field of the \"AGE\" filter card", () => enterCardBound(page, "88.8", "max", "AGE"));
      await session.step(59, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(60, "And the \"filtering of AGE\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of AGE", el("filter panel"), "true"));
      await session.step(61, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(62, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(63, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(64, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Range handles brought back to the ends keep every row (public#1343)", async () => {
      await session.step(68, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(69, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
      await session.step(70, "And user drags the \"range min handle\" area of histogram viewer in filter panel to the \"bin 10\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in filter panel"), "bin 10"));
      await session.step(71, "Then 306 rows should pass the filter", () => filterPasses(page, 306));
      await session.step(72, "When user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
      await session.step(73, "And user drags the \"range max handle\" area of histogram viewer in filter panel to the \"bin 15\" area", () => dragAreaToArea(page, "range max handle", el("histogram viewer in filter panel"), "bin 15"));
      await session.step(74, "Then 260 rows should pass the filter", () => filterPasses(page, 260));
      await session.step(75, "When user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
      await session.step(76, "And user drags the \"range min handle\" area of histogram viewer in filter panel by 400 pixels to the left", () => dragAreaBy(page, "range min handle", el("histogram viewer in filter panel"), 400, "left"));
      await session.step(77, "And user hovers over the \"view\" area of histogram viewer in filter panel", () => hoverArea(page, "view", el("histogram viewer in filter panel")));
      await session.step(78, "And user drags the \"range max handle\" area of histogram viewer in filter panel by 400 pixels to the right", () => dragAreaBy(page, "range max handle", el("histogram viewer in filter panel"), 400, "right"));
      await session.step(79, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(80, "And the \"min of AGE\" reading of filter panel should be 18", () => readingIs(page, "min of AGE", el("filter panel"), 18));
      await session.step(81, "And the \"max of AGE\" reading of filter panel should be 89", () => readingIs(page, "max of AGE", el("filter panel"), 89));
      await session.step(82, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(83, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(84, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A second card's search icon takes the typing while the first search is open (public#2359)", async () => {
      await session.step(88, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(89, "And user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
      await session.step(90, "And user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(91, "And user clicks on search icon of \"RACE\" filter card", () => clickOn(page, el("search icon of \"RACE\" filter card")));
      await session.step(92, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(93, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(94, "And user types \"RA\" at the caret", () => typeAtCaret(page, "RA"));
      await session.step(95, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(96, "And the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
      await session.step(97, "And the \"filtering of DIS_POP\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of DIS_POP", el("filter panel"), "true"));
      await session.step(98, "And the \"filtering of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "false"));
      await session.step(99, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(100, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(101, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A header dragged onto the panel stays where it was in the grid (public#1986)", async () => {
      await session.step(105, "Then the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(106, "When user drags the \"header RACE\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header RACE", el("grid"), "view", el("filter panel")));
      await session.step(107, "Then \"RACE\" filter card should be visible", () => shouldBe(page, el("\"RACE\" filter card"), "visible"));
      await session.step(108, "And the \"cards\" reading of filter panel should be \"RACE\"", () => readingReads(page, "cards", el("filter panel"), "RACE"));
      await session.step(109, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(110, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(111, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(112, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Text selected in a card's search is copied as text (public#3027)", async () => {
      await session.step(116, "When user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
      await session.step(117, "And user hovers over \"DIS_POP\" filter card", () => hoverOver(page, el("\"DIS_POP\" filter card")));
      await session.step(118, "And user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(119, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(120, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(121, "When user presses Control+A in search of \"DIS_POP\" filter card", () => pressKeyIn(page, "Control+A", el("search of \"DIS_POP\" filter card")));
      await session.step(122, "And user presses Control+C", () => pressKey(page, "Control+C"));
      await session.step(123, "Then the clipboard should have text \"Ps\"", () => clipboardHas(page, "Ps"));
      await session.step(124, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(125, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(126, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column holding one value is listed as a row after the cards (core#479)", async () => {
      await session.step(130, "When user adds a calculated column \"probe_class\" with formula '\"low\"'", () => addCalculated(page, "probe_class", "\"low\""));
      await session.step(131, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(132, "And user drags the \"header probe_class\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header probe_class", el("grid"), "view", el("filter panel")));
      await session.step(133, "Then the filter panel should have a filter on \"probe_class\" column", () => filterPanelHas(page, "probe_class"));
      await session.step(134, "And filter panel should have a \"card probe_class\" area", () => hasArea(page, el("filter panel"), "card probe_class"));
      await session.step(135, "And \"probe_class\" filter card should be hidden", () => shouldBe(page, el("\"probe_class\" filter card"), "hidden"));
      await session.step(136, "And the \"card probe_class\" area of filter panel should lie below the \"card RACE\" area", () => areaLies(page, "card probe_class", el("filter panel"), "below", "card RACE"));
      await session.step(137, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(138, "When user picks \"Remove All\" from the filter panel menu", () => pickPanelMenu(page, "Remove All"));
      await session.step(139, "Then the filter panel should have 0 filters", () => filterPanelCount(page, 0));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Adding viewers does not shrink the filter panel (public#2501)", async () => {
      await session.step(144, "When user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
      await session.step(145, "And user remembers the size of filter panel", () => rememberSize(page, el("filter panel")));
      await session.step(146, "And user adds histogram viewer", () => addViewer(page, "histogram"));
      await session.step(147, "And user adds bar chart viewer", () => addViewer(page, "bar chart"));
      await session.step(148, "And user adds line chart viewer", () => addViewer(page, "line chart"));
      await session.step(149, "Then the open tableview should have 1 line chart viewer", () => viewerCount(page, 1, "line chart"));
      await session.step(150, "And filter panel should keep its remembered width", () => keepsWidth(page, el("filter panel")));
      await session.step(151, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The header popup still filters after the column's card left the panel (public#2943, public#1242)", async () => {
      await session.step(155, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(156, "And user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(157, "And user clicks on close of \"AGE\" filter card", () => clickOn(page, el("close of \"AGE\" filter card")));
      await session.step(158, "Then \"AGE\" filter card should be absent", () => shouldBe(page, el("\"AGE\" filter card"), "absent"));
      await session.step(159, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(160, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(161, "And user clicks on \"Column options\" icon in grid", () => clickOn(page, el("\"Column options\" icon in grid")));
      await session.step(162, "Then column popup should be visible", () => shouldBe(page, el("column popup"), "visible"));
      await session.step(163, "When user hovers over the \"view\" area of histogram viewer in column popup", () => hoverArea(page, "view", el("histogram viewer in column popup")));
      await session.step(164, "And user drags the \"range min handle\" area of histogram viewer in column popup to the \"bin 10\" area", () => dragAreaToArea(page, "range min handle", el("histogram viewer in column popup"), "bin 10"));
      await session.step(165, "Then 326 rows should pass the filter", () => filterPasses(page, 326));
      await session.step(166, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(167, "And user resets the filter", () => resetFilter(page));
      await session.step(168, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(169, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
