/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1137, core#1193, core#1195, core#1194, core#1191, core#1196, core#1197, core#1198, core#1200]
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
import {clickOn, enterInto, isExpanded, pressKey, pressKeyIn, shouldBe, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, hasNoColumn, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, noneSelected, removeColumn, renameColumn, rowCount, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, doubleClickArea, hasArea, hoverArea, loadLayout, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingIs, readingReads, saveLayout, setProperty, showsRows, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingContains} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid edge cases", () => {
  const session = feature(test, "features/viewers/grid/grid-edge-cases.feature", import.meta.url);
  test("An HTML summary column keeps rendering after a column it shows is removed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(23, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
    await session.step(24, "And user clicks on All link in \"Select columns...\" dialog", () => clickOn(page, el("All link in \"Select columns...\" dialog")));
    await session.step(25, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(26, "Then the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(27, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
    await session.step(28, "Then the table should not have a column \"CONTROL\"", () => hasNoColumn(page, "CONTROL"));
    await session.step(29, "When user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
    await session.step(30, "And user presses Control+End", () => pressKey(page, "Control+End"));
    await session.step(31, "Then grid should have a \"cell 1000 of html\" area", () => hasArea(page, el("grid"), "cell 1000 of html"));
    await session.step(32, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("Unpin Column with Frozen Columns above the column count logs no error", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:42", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(40, "When user sets \"Frozen Columns\" property of grid to \"20\"", () => setProperty(page, "Frozen Columns", el("grid"), "20"));
      await session.step(41, "And user picks \"Pin > Unpin Column\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin Column", "header AGE", el("grid")));
      await session.step(42, "Then the \"column order\" reading of grid should contain \"AGE\"", () => readingContains(page, "column order", el("grid"), "AGE"));
      await session.step(43, "And grid should have a \"header AGE\" area", () => hasArea(page, el("grid"), "header AGE"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    }, 42);
  });
  test("Frozen Columns typed in the context panel as 20 is not taken [typed=20]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:53", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(48, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(49, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(50, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(51, "Given \"Columns\" category in context panel is expanded", () => isExpanded(page, el("\"Columns\" category in context panel")));
      await session.step(52, "When user enters \"20\" into \"Frozen Columns\" property in context panel", () => enterInto(page, "20", el("\"Frozen Columns\" property in context panel")));
      await session.step(53, "Then \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    }, 53);
  });
  test("Frozen Columns typed in the context panel as -1 is not taken [typed=-1]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:53", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(48, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(49, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(50, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(51, "Given \"Columns\" category in context panel is expanded", () => isExpanded(page, el("\"Columns\" category in context panel")));
      await session.step(52, "When user enters \"-1\" into \"Frozen Columns\" property in context panel", () => enterInto(page, "-1", el("\"Frozen Columns\" property in context panel")));
      await session.step(53, "Then \"Frozen Columns\" property of grid should be \"1\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "1"));
      await session.step(54, "And no errors should have been logged", () => noErrors(page));
    }, 53);
  });
  test("Frozen Columns is brought down to the column count when columns are removed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:74", "@realizes:core#1137"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(64, "When user sets \"Frozen Columns\" property of grid to \"6\"", () => setProperty(page, "Frozen Columns", el("grid"), "6"));
      await session.step(65, "And user removes \"DEMOG\" column", () => removeColumn(page, "DEMOG"));
      await session.step(66, "And user removes \"CONTROL\" column", () => removeColumn(page, "CONTROL"));
      await session.step(67, "And user removes \"STARTED\" column", () => removeColumn(page, "STARTED"));
      await session.step(68, "And user removes \"SEVERITY\" column", () => removeColumn(page, "SEVERITY"));
      await session.step(69, "And user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
      await session.step(70, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
      await session.step(71, "And user removes \"DIS_POP\" column", () => removeColumn(page, "DIS_POP"));
      await session.step(72, "And user removes \"RACE\" column", () => removeColumn(page, "RACE"));
      await session.step(73, "Then the \"column order\" reading of grid should be \"USUBJID, AGE, SEX\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX"));
      await session.step(74, "And \"Frozen Columns\" property of grid should be \"4\"", () => propertyShouldBe(page, "Frozen Columns", el("grid"), "4"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    }, 74);
  });
  test("A value refused in the last row keeps that cell as it was and keeps the added row", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(79, "When user sets \"Add New Row On Last Row Edit\" property of grid to \"true\"", () => setProperty(page, "Add New Row On Last Row Edit", el("grid"), "true"));
    await session.step(80, "And user clicks on the \"cell 3 of AGE\" area of grid", () => clickArea(page, "cell 3 of AGE", el("grid")));
    await session.step(81, "And user presses Control+End", () => pressKey(page, "Control+End"));
    await session.step(82, "And user double-clicks on the \"cell 1000 of AGE\" area of grid", () => doubleClickArea(page, "cell 1000 of AGE", el("grid")));
    await session.step(83, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
    await session.step(84, "And user types \"abc\" into cell editor", () => typeInto(page, "abc", el("cell editor")));
    await session.step(85, "And user presses Enter", () => pressKey(page, "Enter"));
    await session.step(86, "Then the table should have 1001 rows", () => rowCount(page, 1001));
    await session.step(87, "And the value of \"AGE\" column in row 1000 should be \"63\"", () => valueInRow(page, "AGE", 1000, "63"));
    await session.step(88, "And no errors should have been logged", () => noErrors(page));
  });
  test("A read-only grid does not take a paste", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(92, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
    await session.step(93, "And user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
    await session.step(94, "And user presses Control+C", () => pressKey(page, "Control+C"));
    await session.step(95, "And user clicks on the \"cell 11 of AGE\" area of grid", () => clickArea(page, "cell 11 of AGE", el("grid")));
    await session.step(96, "And user presses Control+V", () => pressKey(page, "Control+V"));
    await session.step(97, "Then the value of \"AGE\" column in row 11 should be \"46\"", () => valueInRow(page, "AGE", 11, "46"));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
  });
  test("A read-only grid keeps its rows on Shift+Delete", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:106", "@realizes:core#1193"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(102, "When user sets \"Allow Edit\" property of grid to \"false\"", () => setProperty(page, "Allow Edit", el("grid"), "false"));
      await session.step(103, "And user clicks on the \"row header 5\" area of grid holding Control", () => clickAreaHolding(page, "row header 5", el("grid"), "Control"));
      await session.step(104, "Then 1 row should be selected", () => selectedRowCount(page, 1));
      await session.step(105, "When user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(106, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(107, "And no errors should have been logged", () => noErrors(page));
    }, 106);
  });
  test("Empty values go last in a sort, descending and ascending", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(111, "Given user opens a table \"sort-empty\" with:", () => openTableOf(page, "sort-empty", [["id","v"],["a","2"],["b",""],["c","1"],["d","3"]]), [["id","v"],["a","2"],["b",""],["c","1"],["d","3"]]);
    await session.step(117, "When user double-clicks on the \"header v\" area of grid", () => doubleClickArea(page, "header v", el("grid")));
    await session.step(118, "Then the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
    await session.step(119, "And the \"row order\" reading of grid should be \"4, 1, 3, 2\"", () => readingReads(page, "row order", el("grid"), "4, 1, 3, 2"));
    await session.step(120, "When user double-clicks on the \"header v\" area of grid", () => doubleClickArea(page, "header v", el("grid")));
    await session.step(121, "Then the \"sort direction\" reading of grid should be \"ascending\"", () => readingReads(page, "sort direction", el("grid"), "ascending"));
    await session.step(122, "And the \"row order\" reading of grid should be \"3, 1, 4, 2\"", () => readingReads(page, "row order", el("grid"), "3, 1, 4, 2"));
    await session.step(123, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column group named with spaces only is not made", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:135", "@realizes:core#1195"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(127, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
      await session.step(128, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
      await session.step(129, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(130, "And Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
      await session.step(131, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
      await session.step(132, "Then Group input in dialog should be visible", () => shouldBe(page, el("Group input in dialog"), "visible"));
      await session.step(133, "When user types \"   \" into Group input in dialog", () => typeInto(page, "   ", el("Group input in dialog")));
      await session.step(134, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
      await session.step(135, "Then \"AGE\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "AGE", "group", ""));
      await session.step(136, "And \"SEX\" column should have tag \"group\" equal to \"\"", () => columnTag(page, "SEX", "group", ""));
      await session.step(137, "And no errors should have been logged", () => noErrors(page));
    }, 135);
  });
  test("Shift+Enter selects nothing while Allow Row Selection is off", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:144", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(141, "When user sets \"Allow Row Selection\" property of grid to \"false\"", () => setProperty(page, "Allow Row Selection", el("grid"), "false"));
      await session.step(142, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(143, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(144, "Then no rows should be selected", () => noneSelected(page));
      await session.step(145, "And no errors should have been logged", () => noErrors(page));
    }, 144);
  });
  test("Shift+Enter selects only the same rows that pass the filter", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:152", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(149, "When user filters rows where \"AGE\" is between 20 and 30", () => filterBetween(page, "AGE", 20, 30));
      await session.step(150, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(151, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(152, "Then every selected row should pass the filter", () => selectedPassFilter(page));
      await session.step(153, "And no errors should have been logged", () => noErrors(page));
    }, 152);
  });
  test("Twelve rows pinned show ten, and the grid says that pinned rows are hidden", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:162", "@realizes:core#1191"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(157, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")));
      await session.step(158, "And user clicks on the \"row header 12\" area of grid holding Shift", () => clickAreaHolding(page, "row header 12", el("grid"), "Shift"));
      await session.step(159, "Then 12 rows should be selected", () => selectedRowCount(page, 12));
      await session.step(160, "When user picks \"Pin > Pin Selected Rows\" from the context menu of the \"cell 1 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Selected Rows", "cell 1 of USUBJID", el("grid")));
      await session.step(161, "Then the \"pinned rows\" reading of grid should be 10", () => readingIs(page, "pinned rows", el("grid"), 10));
      await session.step(162, "And a warning balloon containing \"10\" should have been shown", () => warningBalloonText(page, "10"));
    }, 162);
  });
  test("Unpin All Rows unpins the pinned rows the grid does not show", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:174", "@realizes:core#1191"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(166, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")));
      await session.step(167, "And user clicks on the \"row header 12\" area of grid holding Shift", () => clickAreaHolding(page, "row header 12", el("grid"), "Shift"));
      await session.step(168, "And user picks \"Pin > Pin Selected Rows\" from the context menu of the \"cell 1 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Selected Rows", "cell 1 of USUBJID", el("grid")));
      await session.step(169, "Then the \"pinned rows\" reading of grid should be 10", () => readingIs(page, "pinned rows", el("grid"), 10));
      await session.step(170, "When user picks \"Pin > Unpin All Rows\" from the context menu of the \"cell 1 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Unpin All Rows", "cell 1 of USUBJID", el("grid")));
      await session.step(171, "Then the \"pinned rows\" reading of grid should be 0", () => readingIs(page, "pinned rows", el("grid"), 0));
      await session.step(172, "When user saves the layout of the current table view", () => saveLayout(page));
      await session.step(173, "And user loads the saved layout", () => loadLayout(page));
      await session.step(174, "Then the \"pinned rows\" reading of grid should be 0", () => readingIs(page, "pinned rows", el("grid"), 0));
      await session.step(175, "And no errors should have been logged", () => noErrors(page));
    }, 174);
  });
  test("A pinned row deleted and brought back by Control+Z is pinned again", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:186", "@realizes:core#1196"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(179, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 3 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 3 of USUBJID", el("grid")));
      await session.step(180, "Then the \"pinned rows\" reading of grid should be 1", () => readingIs(page, "pinned rows", el("grid"), 1));
      await session.step(181, "When user clicks on the \"row header 3\" area of grid holding Control", () => clickAreaHolding(page, "row header 3", el("grid"), "Control"));
      await session.step(182, "And user presses Shift+Delete", () => pressKey(page, "Shift+Delete"));
      await session.step(183, "Then the table should have 999 rows", () => rowCount(page, 999));
      await session.step(184, "When user presses Control+Z", () => pressKey(page, "Control+Z"));
      await session.step(185, "Then the table should have 1000 rows", () => rowCount(page, 1000));
      await session.step(186, "And the \"pinned rows\" reading of grid should be 1", () => readingIs(page, "pinned rows", el("grid"), 1));
      await session.step(187, "And no errors should have been logged", () => noErrors(page));
    }, 186);
  });
  test("Default HTML Form with no column ticked cannot be confirmed", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:193", "@realizes:core#1197"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(191, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
      await session.step(192, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(193, "And OK button in \"Select columns...\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Select columns...\" dialog"), "disabled"));
    }, 193);
  });
  test("A malformed number format marks the format field invalid", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:200", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(197, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(198, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(199, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(200, "Then Custom input in Format AGE dialog should be invalid", () => shouldBe(page, el("Custom input in Format AGE dialog"), "invalid"));
    }, 200);
  });
  test("A malformed number format is not taken and the old format stays", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:208", "@realizes:core#1198"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(204, "When user picks \"Format > Custom...\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Format > Custom...", "header AGE", el("grid")));
      await session.step(205, "Then Format AGE dialog should be visible", () => shouldBe(page, el("Format AGE dialog"), "visible"));
      await session.step(206, "When user enters \"0.0.0\" into Custom input in Format AGE dialog", () => enterInto(page, "0.0.0", el("Custom input in Format AGE dialog")));
      await session.step(207, "And user clicks on OK button in Format AGE dialog", () => clickOn(page, el("OK button in Format AGE dialog")));
      await session.step(208, "Then \"AGE\" column should have tag \"format\" equal to \"\"", () => columnTag(page, "AGE", "format", ""));
      await session.step(209, "And the \"text of cell 1 of AGE\" reading of grid should be \"26\"", () => readingReads(page, "text of cell 1 of AGE", el("grid"), "26"));
      await session.step(210, "And no errors should have been logged", () => noErrors(page));
    }, 208);
  });
  test("The sort follows a renamed column", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(214, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(215, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(216, "When user renames \"AGE\" column to \"YEARS\"", () => renameColumn(page, "AGE", "YEARS"));
    await session.step(217, "Then the \"sort column\" reading of grid should be \"YEARS\"", () => readingReads(page, "sort column", el("grid"), "YEARS"));
    await session.step(218, "And no errors should have been logged", () => noErrors(page));
  });
  test("A pinned row follows a renamed column through a layout", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(222, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 3 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 3 of USUBJID", el("grid")));
    await session.step(223, "And user renames \"USUBJID\" column to \"ID\"", () => renameColumn(page, "USUBJID", "ID"));
    await session.step(224, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(225, "And user loads the saved layout", () => loadLayout(page));
    await session.step(226, "Then the \"pinned rows\" reading of grid should be 1", () => readingIs(page, "pinned rows", el("grid"), 1));
    await session.step(227, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column group follows a renamed column", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(231, "When user clicks on the \"header AGE\" area of grid holding Control", () => clickAreaHolding(page, "header AGE", el("grid"), "Control"));
    await session.step(232, "And user clicks on the \"header SEX\" area of grid holding Control", () => clickAreaHolding(page, "header SEX", el("grid"), "Control"));
    await session.step(233, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(234, "And Actions accordion header in context panel is expanded", () => isExpanded(page, el("Actions accordion header in context panel")));
    await session.step(235, "When user clicks on \"Group columns...\" text in context panel", () => clickOn(page, el("\"Group columns...\" text in context panel")));
    await session.step(236, "And user types \"Person\" into Group input in dialog", () => typeInto(page, "Person", el("Group input in dialog")));
    await session.step(237, "And user clicks on OK button in dialog", () => clickOn(page, el("OK button in dialog")));
    await session.step(238, "Then grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
    await session.step(239, "When user renames \"AGE\" column to \"YEARS\"", () => renameColumn(page, "AGE", "YEARS"));
    await session.step(240, "Then \"YEARS\" column should have tag \"group\" equal to \"Person\"", () => columnTag(page, "YEARS", "group", "Person"));
    await session.step(241, "And grid should have a \"group Person\" area", () => hasArea(page, el("grid"), "group Person"));
    await session.step(242, "And no errors should have been logged", () => noErrors(page));
  });
  test("Row Height typed in the context panel as 0 is not taken [typed=0]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:251", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(246, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(247, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(248, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(249, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(250, "When user enters \"0\" into \"Row Height\" property in context panel", () => enterInto(page, "0", el("\"Row Height\" property in context panel")));
      await session.step(251, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
      await session.step(252, "And no errors should have been logged", () => noErrors(page));
    }, 251);
  });
  test("Row Height typed in the context panel as 501 is not taken [typed=501]", {tag: ["@viewers", "@realizes:viewers.grid", "@edge-cases", "@known-failure:251", "@realizes:core#1200"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await knownFailure(async () => {
      await session.step(246, "When user hovers over the \"header AGE\" area of grid", () => hoverArea(page, "header AGE", el("grid")));
      await session.step(247, "And user clicks on \"Edit properties (F4)\" icon in grid", () => clickOn(page, el("\"Edit properties (F4)\" icon in grid")));
      await session.step(248, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(249, "Given \"Rows\" category in context panel is expanded", () => isExpanded(page, el("\"Rows\" category in context panel")));
      await session.step(250, "When user enters \"501\" into \"Row Height\" property in context panel", () => enterInto(page, "501", el("\"Row Height\" property in context panel")));
      await session.step(251, "Then \"Row Height\" property of grid should be \"28\"", () => propertyShouldBe(page, "Row Height", el("grid"), "28"));
      await session.step(252, "And no errors should have been logged", () => noErrors(page));
    }, 251);
  });
});
