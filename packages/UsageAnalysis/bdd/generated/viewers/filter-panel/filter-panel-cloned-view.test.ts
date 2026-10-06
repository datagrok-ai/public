/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-cloned-view.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openSharedTable} from '../../../bindings/filter-panel.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, hoverOver, shouldBe, shouldHaveText, shouldNotBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCategoricalFilter, addRangeFilter, filterPasses, filterPassesAll, filterPassesFewer, noneOfFiltered, openEmptyFilterPanel, tableFilterCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openDatasetRowsAs, switchView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, openCardIndicatorMenu, pickCardIndicatorMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, loadLayout, menuLists, noErrors, readingAsRemembered, readingAtLeast, readingNotAsRemembered, readingReads, rememberReading, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragAreaOntoWidget} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Filter panel of a cloned view", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-cloned-view.feature", import.meta.url);
  test("Filter panel of a cloned view", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 8, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(29, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(30, "And user adds a card for \"RACE\" to the filter panel", () => addCardFor(page, "RACE"));
    await session.step(31, "And user adds a card for \"SEX\" to the filter panel", () => addCardFor(page, "SEX"));
    await session.step(32, "And user adds a categorical filter on \"RACE\" keeping \"Caucasian, Black, Other\"", () => addCategoricalFilter(page, "RACE", "Caucasian, Black, Other"));
    await session.step(33, "Then 985 rows should pass the filter", () => filterPasses(page, 985));
    await session.step(34, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
    await run.scenario("A cloned view comes up with the same cards, criteria and switches", async () => {
      await session.step(37, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(38, "Then the \"demog-1000 copy\" view should be current", () => viewIsCurrent(page, "demog-1000 copy"));
      await session.step(39, "And filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(40, "And the \"cards\" reading of filter panel should be \"SEX, RACE\"", () => readingReads(page, "cards", el("filter panel"), "SEX, RACE"));
      await session.step(41, "And the \"selected categories of RACE\" reading of filter panel should be \"Black, Caucasian, Other\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black, Caucasian, Other"));
      await session.step(42, "And the \"filtering of RACE\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "true"));
      await session.step(43, "And the \"filtering of SEX\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of SEX", el("filter panel"), "false"));
      await session.step(44, "And \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(45, "And \"SEX\" filter card should be enabled", () => shouldBe(page, el("\"SEX\" filter card"), "enabled"));
      await session.step(46, "And 985 rows should pass the filter", () => filterPasses(page, 985));
      await session.step(47, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A category clicked on the original's card is followed by the clone's card", async () => {
      await session.step(50, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(51, "And user clicks on the \"category Black of RACE\" area of filter panel", () => clickArea(page, "category Black of RACE", el("filter panel")));
      await session.step(52, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(53, "And the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(54, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(55, "Then the \"selected categories of RACE\" reading of filter panel should be \"Black\"", () => readingReads(page, "selected categories of RACE", el("filter panel"), "Black"));
      await session.step(56, "And the \"filtering of RACE\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of RACE", el("filter panel"), "true"));
      await session.step(57, "And the \"filtering of SEX\" reading of filter panel should be \"false\"", () => readingReads(page, "filtering of SEX", el("filter panel"), "false"));
      await session.step(58, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A card switched off in one view is switched off in the other", async () => {
      await session.step(61, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(62, "And user clicks on the \"category F of SEX\" area of filter panel", () => clickArea(page, "category F of SEX", el("filter panel")));
      await session.step(63, "Then 18 rows should pass the filter", () => filterPasses(page, 18));
      await session.step(64, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(65, "And user unchecks checkbox of \"SEX\" filter card", () => uncheck(page, el("checkbox of \"SEX\" filter card")));
      await session.step(66, "Then 27 rows should pass the filter", () => filterPasses(page, 27));
      await session.step(67, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(68, "Then \"SEX\" filter card should be disabled", () => shouldBe(page, el("\"SEX\" filter card"), "disabled"));
      await session.step(69, "And the \"enabled of SEX\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of SEX", el("filter panel"), "false"));
      await session.step(70, "And \"RACE\" filter card should be enabled", () => shouldBe(page, el("\"RACE\" filter card"), "enabled"));
      await session.step(71, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(72, "And user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(73, "And user checks checkbox of \"SEX\" filter card", () => check(page, el("checkbox of \"SEX\" filter card")));
      await session.step(74, "Then 18 rows should pass the filter", () => filterPasses(page, 18));
      await session.step(75, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(76, "Then \"SEX\" filter card should be enabled", () => shouldBe(page, el("\"SEX\" filter card"), "enabled"));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Switching off a card of another type in the clone leaves the original's card alone", async () => {
      await session.step(80, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(81, "And user adds a range filter on \"AGE\" from 30 to 60", () => addRangeFilter(page, "AGE", 30, 60));
      await session.step(82, "Then fewer than 18 rows should pass the filter", () => filterPassesFewer(page, 18));
      await session.step(83, "And the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(84, "When user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(85, "And user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(86, "And user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(87, "And user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(88, "And user clicks on \"Switch to categorical filter\" icon in \"AGE\" filter card", () => clickOn(page, el("\"Switch to categorical filter\" icon in \"AGE\" filter card")));
      await session.step(89, "Then the \"type of AGE\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of AGE", el("filter panel"), "categorical"));
      await session.step(90, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(91, "And user unchecks checkbox of \"AGE\" filter card", () => uncheck(page, el("checkbox of \"AGE\" filter card")));
      await session.step(92, "Then \"AGE\" filter card should be disabled", () => shouldBe(page, el("\"AGE\" filter card"), "disabled"));
      await session.step(93, "And the \"rows shown\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "rows shown", el("filter panel")));
      await session.step(94, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(95, "Then \"AGE\" filter card should be enabled", () => shouldBe(page, el("\"AGE\" filter card"), "enabled"));
      await session.step(96, "And the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
      await session.step(97, "And the \"rows shown\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "rows shown", el("filter panel")));
      await session.step(98, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Removing a card in the clone leaves the original's card as it was", async () => {
      await session.step(101, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(102, "And user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(103, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(104, "Then \"SEX\" filter card should be absent", () => shouldBe(page, el("\"SEX\" filter card"), "absent"));
      await session.step(105, "And the \"rows shown\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "rows shown", el("filter panel")));
      await session.step(106, "When user switches to the \"demog-1000\" view", () => switchView(page, "demog-1000"));
      await session.step(107, "Then \"SEX\" filter card should be visible", () => shouldBe(page, el("\"SEX\" filter card"), "visible"));
      await session.step(108, "And \"SEX\" filter card should be enabled", () => shouldBe(page, el("\"SEX\" filter card"), "enabled"));
      await session.step(109, "And the \"selected categories of SEX\" reading of filter panel should be \"F\"", () => readingReads(page, "selected categories of SEX", el("filter panel"), "F"));
      await session.step(110, "And the \"rows shown\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "rows shown", el("filter panel")));
      await session.step(111, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A layout saved from the clone brings its cards and switches back", async () => {
      await session.step(114, "When user switches to the \"demog-1000 copy\" view", () => switchView(page, "demog-1000 copy"));
      await session.step(115, "And user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(116, "And user unchecks checkbox of \"RACE\" filter card", () => uncheck(page, el("checkbox of \"RACE\" filter card")));
      await session.step(117, "Then \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(118, "When user remembers the \"rows shown\" reading of filter panel", () => rememberReading(page, "rows shown", el("filter panel")));
      await session.step(119, "And user remembers the \"cards\" reading of filter panel", () => rememberReading(page, "cards", el("filter panel")));
      await session.step(120, "And user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
      await session.step(121, "And user hovers over \"RACE\" filter card", () => hoverOver(page, el("\"RACE\" filter card")));
      await session.step(122, "And user checks checkbox of \"RACE\" filter card", () => check(page, el("checkbox of \"RACE\" filter card")));
      await session.step(123, "Then the \"rows shown\" reading of filter panel should not be as remembered", () => readingNotAsRemembered(page, "rows shown", el("filter panel")));
      await session.step(124, "When user clicks on close icon of filters viewer", () => clickOn(page, el("close icon of filters viewer")));
      await session.step(125, "Then filter panel should be hidden", () => shouldBe(page, el("filter panel"), "hidden"));
      await session.step(126, "When user loads the saved layout", () => loadLayout(page));
      await session.step(127, "Then filter panel should be visible", () => shouldBe(page, el("filter panel"), "visible"));
      await session.step(128, "And \"RACE\" filter card should be disabled", () => shouldBe(page, el("\"RACE\" filter card"), "disabled"));
      await session.step(129, "And the \"enabled of RACE\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of RACE", el("filter panel"), "false"));
      await session.step(130, "And the \"enabled of AGE\" reading of filter panel should be \"false\"", () => readingReads(page, "enabled of AGE", el("filter panel"), "false"));
      await session.step(131, "And the \"cards\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "cards", el("filter panel")));
      await session.step(132, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A missing-values choice made from the card's menu survives the clone", async () => {
      await session.step(135, "When user opens demog-1000 dataset keeping the first 1000 rows as \"demog-missing\"", () => openDatasetRowsAs(page, ds("demog-1000"), 1000, "demog-missing"));
      await session.step(136, "And user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")));
      await session.step(137, "Then \"HEIGHT\" filter card should be visible", () => shouldBe(page, el("\"HEIGHT\" filter card"), "visible"));
      await session.step(138, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(139, "When user opens the indicator menu of the \"HEIGHT\" filter card", () => openCardIndicatorMenu(page, "HEIGHT"));
      await session.step(140, "Then the open menu should list \"Missing values > Keep missing values\"", () => menuLists(page, "Missing values > Keep missing values"));
      await session.step(141, "And \"Keep missing values\" menu item should be selected", () => shouldBe(page, el("\"Keep missing values\" menu item"), "selected"));
      await session.step(142, "And \"Filter out missing values\" menu item should not be selected", () => shouldNotBe(page, el("\"Filter out missing values\" menu item"), "selected"));
      await session.step(143, "And \"Show only missing value\" menu item should not be selected", () => shouldNotBe(page, el("\"Show only missing value\" menu item"), "selected"));
      await session.step(144, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(145, "And user picks \"Missing values | Filter out missing values\" from the indicator menu of the \"HEIGHT\" filter card", () => pickCardIndicatorMenu(page, "Missing values | Filter out missing values", "HEIGHT"));
      await session.step(146, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(147, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(148, "And the \"filtering of HEIGHT\" reading of filter panel should be \"true\"", () => readingReads(page, "filtering of HEIGHT", el("filter panel"), "true"));
      await session.step(149, "When user picks \"View > Layout > Clone View\" from the top menu", () => pickFromTopMenu(page, "View > Layout > Clone View"));
      await session.step(150, "Then the \"demog-missing copy\" view should be current", () => viewIsCurrent(page, "demog-missing copy"));
      await session.step(151, "And 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(152, "When user opens the indicator menu of the \"HEIGHT\" filter card", () => openCardIndicatorMenu(page, "HEIGHT"));
      await session.step(153, "Then the open menu should list \"Missing values > Filter out missing values\"", () => menuLists(page, "Missing values > Filter out missing values"));
      await session.step(154, "And \"Filter out missing values\" menu item should be selected", () => shouldBe(page, el("\"Filter out missing values\" menu item"), "selected"));
      await session.step(155, "And \"Keep missing values\" menu item should not be selected", () => shouldNotBe(page, el("\"Keep missing values\" menu item"), "selected"));
      await session.step(156, "And \"Show only missing value\" menu item should not be selected", () => shouldNotBe(page, el("\"Show only missing value\" menu item"), "selected"));
      await session.step(157, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(158, "And user remembers the \"cards\" reading of filter panel", () => rememberReading(page, "cards", el("filter panel")));
      await session.step(159, "And user drags the \"header USUBJID\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header USUBJID", el("grid"), "view", el("filter panel")));
      await session.step(160, "Then the \"cards\" reading of filter panel should not be as remembered", () => readingNotAsRemembered(page, "cards", el("filter panel")));
      await session.step(161, "And \"USUBJID\" filter card should be visible", () => shouldBe(page, el("\"USUBJID\" filter card"), "visible"));
      await session.step(162, "When user remembers the \"cards\" reading of filter panel", () => rememberReading(page, "cards", el("filter panel")));
      await session.step(163, "And user drags the \"header HEIGHT\" area of grid onto the \"view\" area of filter panel", () => dragAreaOntoWidget(page, "header HEIGHT", el("grid"), "view", el("filter panel")));
      await session.step(164, "Then the \"cards\" reading of filter panel should be as remembered", () => readingAsRemembered(page, "cards", el("filter panel")));
      await session.step(165, "When user switches to the \"demog-missing\" view", () => switchView(page, "demog-missing"));
      await session.step(166, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(167, "When user opens the indicator menu of the \"HEIGHT\" filter card", () => openCardIndicatorMenu(page, "HEIGHT"));
      await session.step(168, "Then the open menu should list \"Missing values > Filter out missing values\"", () => menuLists(page, "Missing values > Filter out missing values"));
      await session.step(169, "And \"Filter out missing values\" menu item should be selected", () => shouldBe(page, el("\"Filter out missing values\" menu item"), "selected"));
      await session.step(170, "And \"Keep missing values\" menu item should not be selected", () => shouldNotBe(page, el("\"Keep missing values\" menu item"), "selected"));
      await session.step(171, "And \"Show only missing value\" menu item should not be selected", () => shouldNotBe(page, el("\"Show only missing value\" menu item"), "selected"));
      await session.step(172, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(173, "Then 872 rows should pass the filter", () => filterPasses(page, 872));
      await session.step(174, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A table built on the same column object is not filtered by the first table's card", async () => {
      await session.step(177, "When user opens a table \"shared SEX\" that shares the \"SEX\" column of the current table", () => openSharedTable(page, "shared SEX", "SEX"));
      await session.step(178, "Then 1000 rows of table \"shared SEX\" should pass the filter", () => tableFilterCount(page, 1000, "shared SEX"));
      await session.step(179, "When user switches to the \"demog-missing\" view", () => switchView(page, "demog-missing"));
      await session.step(180, "And user adds a categorical filter on \"SEX\" keeping \"F\"", () => addCategoricalFilter(page, "SEX", "F"));
      await session.step(181, "Then fewer than 553 rows should pass the filter", () => filterPassesFewer(page, 553));
      await session.step(182, "And the \"rows shown\" reading of filter panel should be at least 1", () => readingAtLeast(page, "rows shown", el("filter panel"), 1));
      await session.step(183, "And no rows where \"SEX\" is \"M\" should pass the filter", () => noneOfFiltered(page, "SEX", "M"));
      await session.step(184, "And 1000 rows of table \"shared SEX\" should pass the filter", () => tableFilterCount(page, 1000, "shared SEX"));
      await session.step(185, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
