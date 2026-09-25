/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/add-new-column/formula-refreshing.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.cp.add-new-column-persists, GROK-17109]
--- */
import {test} from '@playwright/test';
import '../../bindings/enrichment.js';
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '../../bindings/io.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {everyValueEquals, everyValueLog, holdsFormula, noTablesOpen} from '../../bindings/add-new-column.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, pressKey, pressKeyIn, shouldBe, shouldBeSwitchedOn, shouldBecomeVisibleWithin, switchOn, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnTag, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, noProjectOnServer, projectsOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, doubleClickArea, noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("A chain of calculated columns recalculates on a formula edit and survives a project round trip", () => {
  const session = feature(test, "features/add-new-column/formula-refreshing.feature", import.meta.url);
  test("A chain of calculated columns recalculates on a formula edit and survives a project round trip", {tag: ["@journey", "@realizes:powerpack.cp.add-new-column-persists", "@realizes:GROK-17109"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And no project named \"bdd-anc-chain-{run}\" is on the server", () => noProjectOnServer(page, session.text("bdd-anc-chain-{run}")));
    await session.step(20, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(21, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(22, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(23, "When user double-clicks on Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(24, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(25, "And the table should have 5850 rows", () => rowCount(page, 5850));
    await run.scenario("Three columns, each computed from the one before", async () => {
      await session.step(28, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(29, "And user types \"Weight2\" into column name input", () => typeInto(page, "Weight2", el("column name input")));
      await session.step(30, "And user types \"${WEIGHT} + 100\" into formula editor", () => typeInto(page, "${WEIGHT} + 100", el("formula editor")));
      await session.step(31, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(32, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(33, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 100", () => everyValueEquals(page, "Weight2", "WEIGHT", 100));
      await session.step(34, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(35, "And user types \"Weight3\" into column name input", () => typeInto(page, "Weight3", el("column name input")));
      await session.step(36, "And user types \"${Weight2} + 100\" into formula editor", () => typeInto(page, "${Weight2} + 100", el("formula editor")));
      await session.step(37, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(38, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(39, "And every value of \"Weight3\" column should equal \"WEIGHT\" column plus 200", () => everyValueEquals(page, "Weight3", "WEIGHT", 200));
      await session.step(40, "When user clicks on \"Add New Column...\" icon", () => clickOn(page, el("\"Add New Column...\" icon")));
      await session.step(41, "And user types \"Weight4\" into column name input", () => typeInto(page, "Weight4", el("column name input")));
      await session.step(42, "And user types \"Log10(${Weight3}) - 0.2\" into formula editor", () => typeInto(page, "Log10(${Weight3}) - 0.2", el("formula editor")));
      await session.step(43, "And user clicks on OK button in \"Add New Column\" dialog", () => clickOn(page, el("OK button in \"Add New Column\" dialog")));
      await session.step(44, "Then \"Add New Column\" dialog should be hidden", () => shouldBe(page, el("\"Add New Column\" dialog"), "hidden"));
      await session.step(45, "And every value of \"Weight4\" column should be the decimal log of \"Weight3\" column minus 0.2", () => everyValueLog(page, "Weight4", "Weight3", 0.2));
      await session.step(46, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An edit of Weight2 in its Formula pane recalculates Weight3 and Weight4", async () => {
      await session.step(49, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(50, "When user clicks on the \"header Weight2\" area of grid", () => clickArea(page, "header Weight2", el("grid")));
      await session.step(51, "Then the context panel should show \"Weight2\"", () => contextPanelShows(page, "Weight2"));
      await session.step(52, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(53, "Then formula pane editor in Formula pane in context panel should hold the formula \"${WEIGHT} + 100\"", () => holdsFormula(page, el("formula pane editor in Formula pane in context panel"), "${WEIGHT} + 100"));
      await session.step(54, "When user types \"${WEIGHT} + 200\" into formula pane editor in Formula pane in context panel", () => typeInto(page, "${WEIGHT} + 200", el("formula pane editor in Formula pane in context panel")));
      await session.step(55, "And user clicks on Apply button in Formula pane in context panel", () => clickOn(page, el("Apply button in Formula pane in context panel")));
      await session.step(56, "Then \"Weight2\" column should have tag \"formula\" equal to \"${WEIGHT} + 200\"", () => columnTag(page, "Weight2", "formula", "${WEIGHT} + 200"));
      await session.step(57, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 200", () => everyValueEquals(page, "Weight2", "WEIGHT", 200));
      await session.step(58, "And every value of \"Weight3\" column should equal \"WEIGHT\" column plus 300", () => everyValueEquals(page, "Weight3", "WEIGHT", 300));
      await session.step(59, "And every value of \"Weight4\" column should be the decimal log of \"Weight3\" column minus 0.2", () => everyValueLog(page, "Weight4", "Weight3", 0.2));
      await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(61, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An edit of Weight3 recalculates Weight4 and leaves Weight2", async () => {
      await session.step(64, "When user clicks on the \"header Weight3\" area of grid", () => clickArea(page, "header Weight3", el("grid")));
      await session.step(65, "Then the context panel should show \"Weight3\"", () => contextPanelShows(page, "Weight3"));
      await session.step(66, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(67, "When user types \"${Weight2} + 50\" into formula pane editor in Formula pane in context panel", () => typeInto(page, "${Weight2} + 50", el("formula pane editor in Formula pane in context panel")));
      await session.step(68, "And user clicks on Apply button in Formula pane in context panel", () => clickOn(page, el("Apply button in Formula pane in context panel")));
      await session.step(69, "Then \"Weight3\" column should have tag \"formula\" equal to \"${Weight2} + 50\"", () => columnTag(page, "Weight3", "formula", "${Weight2} + 50"));
      await session.step(70, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 50", () => everyValueEquals(page, "Weight3", "Weight2", 50));
      await session.step(71, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 200", () => everyValueEquals(page, "Weight2", "WEIGHT", 200));
      await session.step(72, "And every value of \"Weight4\" column should be the decimal log of \"Weight3\" column minus 0.2", () => everyValueLog(page, "Weight4", "Weight3", 0.2));
      await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(74, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("An edit of Weight4 changes Weight4 only", async () => {
      await session.step(77, "When user clicks on the \"header Weight4\" area of grid", () => clickArea(page, "header Weight4", el("grid")));
      await session.step(78, "Then the context panel should show \"Weight4\"", () => contextPanelShows(page, "Weight4"));
      await session.step(79, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(80, "When user types \"Log10(${Weight3}) - 0.1\" into formula pane editor in Formula pane in context panel", () => typeInto(page, "Log10(${Weight3}) - 0.1", el("formula pane editor in Formula pane in context panel")));
      await session.step(81, "And user clicks on Apply button in Formula pane in context panel", () => clickOn(page, el("Apply button in Formula pane in context panel")));
      await session.step(82, "Then \"Weight4\" column should have tag \"formula\" equal to \"Log10(${Weight3}) - 0.1\"", () => columnTag(page, "Weight4", "formula", "Log10(${Weight3}) - 0.1"));
      await session.step(83, "And every value of \"Weight4\" column should be the decimal log of \"Weight3\" column minus 0.1", () => everyValueLog(page, "Weight4", "Weight3", 0.1));
      await session.step(84, "And every value of \"Weight3\" column should equal \"WEIGHT\" column plus 250", () => everyValueEquals(page, "Weight3", "WEIGHT", 250));
      await session.step(85, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 200", () => everyValueEquals(page, "Weight2", "WEIGHT", 200));
      await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(87, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Saved with Data sync, closed and reopened, the chain comes back as last edited", async () => {
      await session.step(90, "When user double-clicks on the \"cell 1 of WEIGHT\" area of grid", () => doubleClickArea(page, "cell 1 of WEIGHT", el("grid")));
      await session.step(91, "And user presses Control+A in cell editor", () => pressKeyIn(page, "Control+A", el("cell editor")));
      await session.step(92, "And user types \"500\" into cell editor", () => typeInto(page, "500", el("cell editor")));
      await session.step(93, "And user presses Enter", () => pressKey(page, "Enter"));
      await session.step(94, "Then the value of \"WEIGHT\" column in row 1 should be \"500\"", () => valueInRow(page, "WEIGHT", 1, "500"));
      await session.step(95, "And the value of \"Weight2\" column in row 1 should be \"700\"", () => valueInRow(page, "Weight2", 1, "700"));
      await session.step(96, "And the value of \"Weight3\" column in row 1 should be \"750\"", () => valueInRow(page, "Weight3", 1, "750"));
      await session.step(97, "When user clicks on Save button in toolbar", () => clickOn(page, el("Save button in toolbar")));
      await session.step(98, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(99, "When user enters \"bdd-anc-chain-{run}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("bdd-anc-chain-{run}"), el("Name text input in \"Save project\" dialog")));
      await session.step(100, "And user switches on Data sync input in \"Save project\" dialog", () => switchOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(101, "Then Data sync input in \"Save project\" dialog should be switched on", () => shouldBeSwitchedOn(page, el("Data sync input in \"Save project\" dialog")));
      await session.step(102, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(103, "Then \"Save project\" dialog should be hidden", () => shouldBe(page, el("\"Save project\" dialog"), "hidden"));
      await session.step(104, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(105, "And 1 project named \"bdd-anc-chain-{run}\" should be on the server", () => projectsOnServer(page, 1, session.text("bdd-anc-chain-{run}")));
      await session.step(106, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(107, "And user closes all views", () => closeAllViews(page));
      await session.step(108, "Then no table should be open", () => noTablesOpen(page));
      await session.step(109, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(110, "Then the \"Projects\" view should be current", () => viewIsCurrent(page, "Projects"));
      await session.step(111, "When user types \"bdd-anc-chain-{run}\" into gallery search", () => typeInto(page, session.text("bdd-anc-chain-{run}"), el("gallery search")));
      await session.step(112, "Then \"bdd-anc-chain-{run}\" project card should become visible within 60 seconds", () => shouldBecomeVisibleWithin(page, el(session.text("\"bdd-anc-chain-{run}\" project card")), 60));
      await session.step(113, "When user double-clicks on \"bdd-anc-chain-{run}\" project card", () => doubleClickOn(page, el(session.text("\"bdd-anc-chain-{run}\" project card"))));
      await session.step(114, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(115, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(117, "And the value of \"WEIGHT\" column in row 1 should be \"73.19999694824219\"", () => valueInRow(page, "WEIGHT", 1, "73.19999694824219"));
      await session.step(118, "And the value of \"Weight2\" column in row 1 should be \"273.20001220703125\"", () => valueInRow(page, "Weight2", 1, "273.20001220703125"));
      await session.step(119, "And \"Weight2\" column should have tag \"formula\" equal to \"${WEIGHT} + 200\"", () => columnTag(page, "Weight2", "formula", "${WEIGHT} + 200"));
      await session.step(120, "And \"Weight3\" column should have tag \"formula\" equal to \"${Weight2} + 50\"", () => columnTag(page, "Weight3", "formula", "${Weight2} + 50"));
      await session.step(121, "And \"Weight4\" column should have tag \"formula\" equal to \"Log10(${Weight3}) - 0.1\"", () => columnTag(page, "Weight4", "formula", "Log10(${Weight3}) - 0.1"));
      await session.step(122, "And every value of \"Weight2\" column should equal \"WEIGHT\" column plus 200", () => everyValueEquals(page, "Weight2", "WEIGHT", 200));
      await session.step(123, "And every value of \"Weight3\" column should equal \"Weight2\" column plus 50", () => everyValueEquals(page, "Weight3", "Weight2", 50));
      await session.step(124, "And every value of \"Weight4\" column should be the decimal log of \"Weight3\" column minus 0.1", () => everyValueLog(page, "Weight4", "Weight3", 0.1));
      await session.step(125, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(126, "When user clicks on the \"header Weight3\" area of grid", () => clickArea(page, "header Weight3", el("grid")));
      await session.step(127, "Then the context panel should show \"Weight3\"", () => contextPanelShows(page, "Weight3"));
      await session.step(128, "Given Formula pane in context panel is expanded", () => isExpanded(page, el("Formula pane in context panel")));
      await session.step(129, "Then formula pane editor in Formula pane in context panel should hold the formula \"${Weight2} + 50\"", () => holdsFormula(page, el("formula pane editor in Formula pane in context panel"), "${Weight2} + 50"));
      await session.step(130, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(131, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
