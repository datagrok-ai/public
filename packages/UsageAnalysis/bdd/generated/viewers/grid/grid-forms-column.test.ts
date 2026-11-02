/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-forms-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#1197, core#1194]
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
import {clickOn, pressKey, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {hasNoColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterBetween, noneSelected, selectedPassFilter, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {contextPanelOpen, contextPanelShows, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hasArea, noErrors, pickFromAreaContextMenu, readingReads, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid form columns", () => {
  const session = feature(test, "features/viewers/grid/grid-forms-column.feature", import.meta.url);
  test("Design a Form... adds a form column and opens its designer, and Edit reopens it", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(37, "When user picks \"Add > Summary Columns > Design a Form...\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Design a Form...", "cell 2 of USUBJID", el("grid")));
    await session.step(38, "Then the \"Summary\" view should be current", () => viewIsCurrent(page, "Summary"));
    await session.step(39, "And \"CLOSE AND APPLY\" button should be visible", () => shouldBe(page, el("\"CLOSE AND APPLY\" button"), "visible"));
    await session.step(40, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
    await session.step(41, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(42, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form"));
    await session.step(43, "And the \"cell type of form\" reading of grid should be \"form\"", () => readingReads(page, "cell type of form", el("grid"), "form"));
    await session.step(44, "When user clicks on the \"header form\" area of grid", () => clickArea(page, "header form", el("grid")));
    await session.step(45, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(46, "Then the context panel should show \"form\"", () => contextPanelShows(page, "form"));
    await session.step(47, "And context panel should contain text \"Renderer\"", () => shouldContainText(page, el("context panel"), "Renderer"));
    await session.step(48, "And context panel should contain text \"Actions\"", () => shouldContainText(page, el("context panel"), "Actions"));
    await session.step(49, "When user clicks on EDIT button in context panel", () => clickOn(page, el("EDIT button in context panel")));
    await session.step(50, "Then the \"Summary\" view should be current", () => viewIsCurrent(page, "Summary"));
    await session.step(51, "When user clicks on \"CLOSE AND APPLY\" button", () => clickOn(page, el("\"CLOSE AND APPLY\" button")));
    await session.step(52, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(53, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, form"));
    await session.step(54, "And the \"cell type of form\" reading of grid should be \"form\"", () => readingReads(page, "cell type of form", el("grid"), "form"));
    await session.step(55, "And no errors should have been logged", () => noErrors(page));
  });
  test("Default HTML Form asks for the columns and adds an HTML column, which keeps rendering after a column it shows is removed", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(58, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
    await session.step(59, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(60, "When user clicks on All link in \"Select columns...\" dialog", () => clickOn(page, el("All link in \"Select columns...\" dialog")));
    await session.step(61, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(62, "Then \"Select columns...\" dialog should be hidden", () => shouldBe(page, el("\"Select columns...\" dialog"), "hidden"));
    await session.step(63, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html"));
    await session.step(64, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(65, "When user picks \"Remove\" from the context menu of the \"header CONTROL\" area of grid", () => pickFromAreaContextMenu(page, "Remove", "header CONTROL", el("grid")));
    await session.step(66, "Then the table should not have a column \"CONTROL\"", () => hasNoColumn(page, "CONTROL"));
    await session.step(67, "When user clicks on the \"cell 2 of AGE\" area of grid", () => clickArea(page, "cell 2 of AGE", el("grid")));
    await session.step(68, "And user presses Control+End", () => pressKey(page, "Control+End"));
    await session.step(69, "Then grid should have a \"cell 1000 of html\" area", () => hasArea(page, el("grid"), "cell 1000 of html"));
    await session.step(70, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(71, "And no errors should have been logged", () => noErrors(page));
  });
  test("Custom HTML Form... takes the markup and adds an HTML column", {tag: ["@viewers", "@realizes:viewers.grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await session.step(74, "When user picks \"Add > Summary Columns > Custom HTML Form...\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Custom HTML Form...", "cell 2 of USUBJID", el("grid")));
    await session.step(75, "Then \"Add Custom Form\" dialog should be visible", () => shouldBe(page, el("\"Add Custom Form\" dialog"), "visible"));
    await session.step(76, "When user clicks on OK button in \"Add Custom Form\" dialog", () => clickOn(page, el("OK button in \"Add Custom Form\" dialog")));
    await session.step(77, "Then \"Add Custom Form\" dialog should be hidden", () => shouldBe(page, el("\"Add Custom Form\" dialog"), "hidden"));
    await session.step(78, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY, html"));
    await session.step(79, "And the \"cell type of html\" reading of grid should be \"html\"", () => readingReads(page, "cell type of html", el("grid"), "html"));
    await session.step(80, "And no errors should have been logged", () => noErrors(page));
  });
  test("Default HTML Form with no column ticked cannot be confirmed", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:86", "@realizes:core#1197"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await knownFailure(async () => {
      await session.step(84, "When user picks \"Add > Summary Columns > Default HTML Form\" from the context menu of the \"cell 2 of USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Add > Summary Columns > Default HTML Form", "cell 2 of USUBJID", el("grid")));
      await session.step(85, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
      await session.step(86, "And OK button in \"Select columns...\" dialog should be disabled", () => shouldBe(page, el("OK button in \"Select columns...\" dialog"), "disabled"));
    }, 86);
  });
  test("Shift+Enter selects nothing while Allow Row Selection is off", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:93", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await knownFailure(async () => {
      await session.step(90, "When user sets \"Allow Row Selection\" property of grid to \"false\"", () => setProperty(page, "Allow Row Selection", el("grid"), "false"));
      await session.step(91, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(92, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(93, "Then no rows should be selected", () => noneSelected(page));
    }, 93);
  });
  test("Shift+Enter selects only the same rows that pass the filter", {tag: ["@viewers", "@realizes:viewers.grid", "@known-failure:100", "@realizes:core#1194"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(31, "Given user is logged in", () => loggedIn(page));
    await session.step(32, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(33, "Then grid should show 1000 rows", () => showsRows(page, el("grid"), 1000));
    await session.step(34, "And the \"column order\" reading of grid should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("grid"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
    await knownFailure(async () => {
      await session.step(97, "When user filters rows where \"AGE\" is between 20 and 30", () => filterBetween(page, "AGE", 20, 30));
      await session.step(98, "And user clicks on the \"cell 1 of SEX\" area of grid", () => clickArea(page, "cell 1 of SEX", el("grid")));
      await session.step(99, "And user presses Shift+Enter", () => pressKey(page, "Shift+Enter"));
      await session.step(100, "Then 82 rows should be selected", () => selectedRowCount(page, 82));
      await session.step(101, "And every selected row should pass the filter", () => selectedPassFilter(page));
    }, 100);
  });
});
