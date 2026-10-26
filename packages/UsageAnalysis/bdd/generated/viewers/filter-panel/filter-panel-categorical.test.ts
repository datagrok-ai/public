/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/filter-panel/filter-panel-categorical.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.filters]
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
import {clickOn, hoverOver, pressKey, shouldBe, shouldHaveText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {addCalculated, allOfFiltered, clearSelection, filterIsExactlyCategory, filterIsExactlyContains, filterPasses, filterPassesAll, noneOfFiltered, openEmptyFilterPanel, removeColumn, rowCount, selectWhereIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addCardFor, clearCardSearch, openCardIndicatorMenu, pasteIntoCardSearch, pickCardIndicatorMenu, typeIntoCardSearch} from '@datagrok-libraries/bdd/bindings/tiers/viewers/filter-panel';
import {clickArea, closeContextMenu, hasArea, hasNoArea, menuDoesNotList, menuLists, noErrors, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Categorical filter card", () => {
  const session = feature(test, "features/viewers/filter-panel/filter-panel-categorical.feature", import.meta.url);
  test("Categorical filter card", {tag: ["@journey", "@viewers", "@realizes:viewers.filters"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(26, "And user opens an empty filter panel", () => openEmptyFilterPanel(page));
    await session.step(27, "And user adds a card for \"DIS_POP\" to the filter panel", () => addCardFor(page, "DIS_POP"));
    await session.step(28, "Then the \"type of DIS_POP\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of DIS_POP", el("filter panel"), "categorical"));
    await session.step(29, "And the \"categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis, RA, UC\"", () => readingReads(page, "categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis, RA, UC"));
    await session.step(30, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(31, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
    await run.scenario("A name click keeps one category, a checkbox click adds another", async () => {
      await session.step(34, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
      await session.step(35, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(36, "And the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
      await session.step(37, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
      await session.step(38, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(39, "When user clicks on the \"checkbox UC of DIS_POP\" area of filter panel", () => clickArea(page, "checkbox UC of DIS_POP", el("filter panel")));
      await session.step(40, "Then 557 rows should pass the filter", () => filterPasses(page, 557));
      await session.step(41, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA, UC"));
      await session.step(42, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The indicator menu selects, deselects and inverts every category", async () => {
      await session.step(46, "When user picks \"Deselect all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Deselect all", "DIS_POP"));
      await session.step(47, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(48, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), ""));
      await session.step(49, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(50, "When user picks \"Invert all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Invert all", "DIS_POP"));
      await session.step(51, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(52, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis, RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis, RA, UC"));
      await session.step(53, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(54, "When user picks \"Deselect all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Deselect all", "DIS_POP"));
      await session.step(55, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
      await session.step(56, "When user picks \"Select all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Select all", "DIS_POP"));
      await session.step(57, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(58, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(59, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Radio mode keeps exactly one category and offers no batch operations", async () => {
      await session.step(62, "When user clicks on the \"category RA of DIS_POP\" area of filter panel", () => clickArea(page, "category RA of DIS_POP", el("filter panel")));
      await session.step(63, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(64, "When user clicks on the \"checkbox UC of DIS_POP\" area of filter panel", () => clickArea(page, "checkbox UC of DIS_POP", el("filter panel")));
      await session.step(65, "Then 557 rows should pass the filter", () => filterPasses(page, 557));
      await session.step(66, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA, UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA, UC"));
      await session.step(67, "When user picks \"Mode | Radio\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Mode | Radio", "DIS_POP"));
      await session.step(68, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(69, "Then 434 rows should pass the filter", () => filterPasses(page, 434));
      await session.step(70, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"RA\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "RA"));
      await session.step(71, "When user opens the indicator menu of the \"DIS_POP\" filter card", () => openCardIndicatorMenu(page, "DIS_POP"));
      await session.step(72, "Then the open menu should list \"Mode\"", () => menuLists(page, "Mode"));
      await session.step(73, "And the open menu should not list \"Select all\"", () => menuDoesNotList(page, "Select all"));
      await session.step(74, "And the open menu should not list \"Deselect all\"", () => menuDoesNotList(page, "Deselect all"));
      await session.step(75, "And the open menu should not list \"Invert all\"", () => menuDoesNotList(page, "Invert all"));
      await session.step(76, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(77, "And user clicks on the \"checkbox UC of DIS_POP\" area of filter panel", () => clickArea(page, "checkbox UC of DIS_POP", el("filter panel")));
      await session.step(78, "Then 123 rows should pass the filter", () => filterPasses(page, 123));
      await session.step(79, "And the \"selected categories of DIS_POP\" reading of filter panel should be \"UC\"", () => readingReads(page, "selected categories of DIS_POP", el("filter panel"), "UC"));
      await session.step(80, "And counter of filter panel should be visible", () => shouldBe(page, el("counter of filter panel"), "visible"));
      await session.step(81, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(82, "When user picks \"Mode | Multi-Select\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Mode | Multi-Select", "DIS_POP"));
      await session.step(83, "And user closes the context menu", () => closeContextMenu(page));
      await session.step(84, "Then 123 rows should pass the filter", () => filterPasses(page, 123));
      await session.step(85, "When user opens the indicator menu of the \"DIS_POP\" filter card", () => openCardIndicatorMenu(page, "DIS_POP"));
      await session.step(86, "Then the open menu should list \"Select all\"", () => menuLists(page, "Select all"));
      await session.step(87, "And the open menu should list \"Deselect all\"", () => menuLists(page, "Deselect all"));
      await session.step(88, "And the open menu should list \"Invert all\"", () => menuLists(page, "Invert all"));
      await session.step(89, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(90, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The in-card search narrows the table to the categories that match", async () => {
      await session.step(93, "When user picks \"Select all\" from the indicator menu of the \"DIS_POP\" filter card", () => pickCardIndicatorMenu(page, "Select all", "DIS_POP"));
      await session.step(94, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(95, "When user clicks on search icon of \"DIS_POP\" filter card", () => clickOn(page, el("search icon of \"DIS_POP\" filter card")));
      await session.step(96, "And user types \"Ps\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "Ps", "DIS_POP"));
      await session.step(97, "Then 242 rows should pass the filter", () => filterPasses(page, 242));
      await session.step(98, "And the filter should pass exactly the rows where \"DIS_POP\" contains \"Ps\"", () => filterIsExactlyContains(page, "DIS_POP", "Ps"));
      await session.step(99, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(100, "When user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
      await session.step(101, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(102, "When user types \"ori\" into the search of the \"DIS_POP\" filter card", () => typeIntoCardSearch(page, "ori", "DIS_POP"));
      await session.step(103, "Then 204 rows should pass the filter", () => filterPasses(page, 204));
      await session.step(104, "And the filter should pass exactly the rows where \"DIS_POP\" contains \"ori\"", () => filterIsExactlyContains(page, "DIS_POP", "ori"));
      await session.step(105, "When user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
      await session.step(106, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(107, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(108, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A pasted list of lines keeps exactly the categories it lists", async () => {
      await session.step(111, "When user pastes \"RA\\nUC\\n\" into the search of the \"DIS_POP\" filter card", () => pasteIntoCardSearch(page, "RA\\nUC\\n", "DIS_POP"));
      await session.step(112, "Then 557 rows should pass the filter", () => filterPasses(page, 557));
      await session.step(113, "And no rows where \"DIS_POP\" is \"Psoriasis\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "Psoriasis"));
      await session.step(114, "And no rows where \"DIS_POP\" is \"PsA\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "PsA"));
      await session.step(115, "And all rows where \"DIS_POP\" is \"UC\" should pass the filter", () => allOfFiltered(page, "DIS_POP", "UC"));
      await session.step(116, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(117, "When user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
      await session.step(118, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(119, "When user pastes \"RA\\nUC\" into the search of the \"DIS_POP\" filter card", () => pasteIntoCardSearch(page, "RA\\nUC", "DIS_POP"));
      await session.step(120, "Then 557 rows should pass the filter", () => filterPasses(page, 557));
      await session.step(121, "And no rows where \"DIS_POP\" is \"Psoriasis\" should pass the filter", () => noneOfFiltered(page, "DIS_POP", "Psoriasis"));
      await session.step(122, "When user clears the search of the \"DIS_POP\" filter card", () => clearCardSearch(page, "DIS_POP"));
      await session.step(123, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(124, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(125, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Rows deleted and brought back by Ctrl+Z bring their category back to the card", async () => {
      await session.step(128, "When user selects rows where \"DIS_POP\" is \"PsA\"", () => selectWhereIs(page, "DIS_POP", "PsA"));
      await session.step(129, "And user picks \"Edit > Remove > Selected Rows\" from the top menu", () => pickFromTopMenu(page, "Edit > Remove > Selected Rows"));
      await session.step(130, "Then the table should have 962 rows", () => rowCount(page, 962));
      await session.step(131, "And the \"categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, Psoriasis, RA, UC\"", () => readingReads(page, "categories of DIS_POP", el("filter panel"), "AS, Indigestion, Psoriasis, RA, UC"));
      await session.step(132, "And filter panel should not have a \"category PsA of DIS_POP\" area", () => hasNoArea(page, el("filter panel"), "category PsA of DIS_POP"));
      await session.step(133, "When user clicks on grid", () => clickOn(page, el("grid")));
      await session.step(134, "And user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(135, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(136, "And the \"categories of DIS_POP\" reading of filter panel should be \"AS, Indigestion, PsA, Psoriasis, RA, UC\"", () => readingReads(page, "categories of DIS_POP", el("filter panel"), "AS, Indigestion, PsA, Psoriasis, RA, UC"));
      await session.step(137, "And filter panel should have a \"category PsA of DIS_POP\" area", () => hasArea(page, el("filter panel"), "category PsA of DIS_POP"));
      await session.step(138, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(139, "When user clears the row selection", () => clearSelection(page));
      await session.step(140, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A numeric card switches to categorical and back without filtering", async () => {
      await session.step(143, "When user adds a card for \"AGE\" to the filter panel", () => addCardFor(page, "AGE"));
      await session.step(144, "Then the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
      await session.step(145, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(146, "And user clicks on \"Switch to categorical filter\" icon in \"AGE\" filter card", () => clickOn(page, el("\"Switch to categorical filter\" icon in \"AGE\" filter card")));
      await session.step(147, "Then the \"type of AGE\" reading of filter panel should be \"categorical\"", () => readingReads(page, "type of AGE", el("filter panel"), "categorical"));
      await session.step(148, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(149, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(150, "When user hovers over caption of \"AGE\" filter card", () => hoverOver(page, el("caption of \"AGE\" filter card")));
      await session.step(151, "Then indicator of \"AGE\" filter card should be visible", () => shouldBe(page, el("indicator of \"AGE\" filter card"), "visible"));
      await session.step(152, "When user opens the indicator menu of the \"AGE\" filter card", () => openCardIndicatorMenu(page, "AGE"));
      await session.step(153, "Then the open menu should list \"Select all\"", () => menuLists(page, "Select all"));
      await session.step(154, "When user closes the context menu", () => closeContextMenu(page));
      await session.step(155, "And user clicks on the \"category 20 of AGE\" area of filter panel", () => clickArea(page, "category 20 of AGE", el("filter panel")));
      await session.step(156, "Then the filter should pass exactly the rows where \"AGE\" is \"20\"", () => filterIsExactlyCategory(page, "AGE", "20"));
      await session.step(157, "And the \"selected categories of AGE\" reading of filter panel should be \"20\"", () => readingReads(page, "selected categories of AGE", el("filter panel"), "20"));
      await session.step(158, "And counter of filter panel should have text \"1\"", () => shouldHaveText(page, el("counter of filter panel"), "1"));
      await session.step(159, "When user hovers over \"AGE\" filter card", () => hoverOver(page, el("\"AGE\" filter card")));
      await session.step(160, "And user clicks on \"Switch to histogram filter\" icon in \"AGE\" filter card", () => clickOn(page, el("\"Switch to histogram filter\" icon in \"AGE\" filter card")));
      await session.step(161, "Then the \"type of AGE\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of AGE", el("filter panel"), "histogram"));
      await session.step(162, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(163, "And counter of filter panel should be hidden", () => shouldBe(page, el("counter of filter panel"), "hidden"));
      await session.step(164, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column holding one value gets a card that takes a click without an error", async () => {
      await session.step(167, "When user adds a calculated column \"probe_constant\" with formula \"1\"", () => addCalculated(page, "probe_constant", "1"));
      await session.step(168, "And user adds a card for \"probe_constant\" to the filter panel", () => addCardFor(page, "probe_constant"));
      await session.step(169, "Then \"probe_constant\" filter card should be visible", () => shouldBe(page, el("\"probe_constant\" filter card"), "visible"));
      await session.step(170, "And the \"type of probe_constant\" reading of filter panel should be \"histogram\"", () => readingReads(page, "type of probe_constant", el("filter panel"), "histogram"));
      await session.step(171, "When user clicks on body of \"probe_constant\" filter card", () => clickOn(page, el("body of \"probe_constant\" filter card")));
      await session.step(172, "Then all rows should pass the filter", () => filterPassesAll(page));
      await session.step(173, "And no errors should have been logged", () => noErrors(page));
      await session.step(174, "When user removes \"probe_constant\" column", () => removeColumn(page, "probe_constant"));
      await session.step(175, "Then \"probe_constant\" filter card should be absent", () => shouldBe(page, el("\"probe_constant\" filter card"), "absent"));
      await session.step(176, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(177, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A \"(\" typed into the in-card search keeps the categories holding one", async () => {
      await session.step(180, "When user adds a calculated column \"AGE_BAND\" with formula \"${SEX} + if(${AGE} > 60, \\\" (60+)\\\", \\\" [..60]\\\")\"", () => addCalculated(page, "AGE_BAND", "${SEX} + if(${AGE} > 60, \" (60+)\", \" [..60]\")"));
      await session.step(181, "And user adds a card for \"AGE_BAND\" to the filter panel", () => addCardFor(page, "AGE_BAND"));
      await session.step(182, "Then the \"categories of AGE_BAND\" reading of filter panel should be \"F (60+), F [..60], M (60+), M [..60]\"", () => readingReads(page, "categories of AGE_BAND", el("filter panel"), "F (60+), F [..60], M (60+), M [..60]"));
      await session.step(183, "When user hovers over \"AGE_BAND\" filter card", () => hoverOver(page, el("\"AGE_BAND\" filter card")));
      await session.step(184, "And user clicks on search icon of \"AGE_BAND\" filter card", () => clickOn(page, el("search icon of \"AGE_BAND\" filter card")));
      await session.step(185, "And user types \"(\" into the search of the \"AGE_BAND\" filter card", () => typeIntoCardSearch(page, "(", "AGE_BAND"));
      await session.step(186, "Then the \"selected categories of AGE_BAND\" reading of filter panel should be \"F (60+), M (60+)\"", () => readingReads(page, "selected categories of AGE_BAND", el("filter panel"), "F (60+), M (60+)"));
      await session.step(187, "And 153 rows should pass the filter", () => filterPasses(page, 153));
      await session.step(188, "And no errors should have been logged", () => noErrors(page));
      await session.step(189, "When user removes \"AGE_BAND\" column", () => removeColumn(page, "AGE_BAND"));
      await session.step(190, "Then \"AGE_BAND\" filter card should be absent", () => shouldBe(page, el("\"AGE_BAND\" filter card"), "absent"));
      await session.step(191, "And all rows should pass the filter", () => filterPassesAll(page));
      await session.step(192, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
