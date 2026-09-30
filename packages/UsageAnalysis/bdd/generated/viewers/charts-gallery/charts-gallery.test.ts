/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/charts-gallery/charts-gallery.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.surface-plot, charts.viewer.globe, charts.viewer.group-analysis]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, followingShouldBe, hoverOver, isExpanded, pressKey, selectIn, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, openDataset, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hasNoArea, loadLayout, noErrors, propertyShouldBe, readingIs, readingReads, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Add viewer gallery for Charts viewers, Surface plot, Globe and Group Analysis", () => {
  const session = feature(test, "features/viewers/charts-gallery/charts-gallery.feature", import.meta.url);
  test("The gallery disables the Charts viewers a table cannot feed and says why", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(22, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(23, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(24, "And the following elements should be enabled:", () => followingShouldBe(page, "enabled", [["first \"Sankey\" card in \"Add Viewer\" dialog"],["first \"Chord\" card in \"Add Viewer\" dialog"],["first \"Timelines\" card in \"Add Viewer\" dialog"],["first \"Radar\" card in \"Add Viewer\" dialog"],["first \"Sunburst\" card in \"Add Viewer\" dialog"],["first \"Tree\" card in \"Add Viewer\" dialog"],["first \"Surface plot\" card in \"Add Viewer\" dialog"],["first \"Globe\" card in \"Add Viewer\" dialog"],["first \"Group Analysis\" card in \"Add Viewer\" dialog"],["first \"Word cloud\" card in \"Add Viewer\" dialog"]]), [["first \"Sankey\" card in \"Add Viewer\" dialog"],["first \"Chord\" card in \"Add Viewer\" dialog"],["first \"Timelines\" card in \"Add Viewer\" dialog"],["first \"Radar\" card in \"Add Viewer\" dialog"],["first \"Sunburst\" card in \"Add Viewer\" dialog"],["first \"Tree\" card in \"Add Viewer\" dialog"],["first \"Surface plot\" card in \"Add Viewer\" dialog"],["first \"Globe\" card in \"Add Viewer\" dialog"],["first \"Group Analysis\" card in \"Add Viewer\" dialog"],["first \"Word cloud\" card in \"Add Viewer\" dialog"]]);
    await session.step(35, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(36, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(37, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(38, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(39, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(40, "And Files---App-Data---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Chem tree node inside browse tree")));
    await session.step(41, "When user double-clicks on Files---App-Data---Chem---chem_standards.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Chem---chem_standards.csv tree node inside browse tree")));
    await session.step(42, "Then the \"chem_standards\" view should be current", () => viewIsCurrent(page, "chem_standards"));
    await session.step(43, "And the table should have 2 columns", () => columnCount(page, 2));
    await session.step(44, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(45, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(46, "And first \"Radar\" card in \"Add Viewer\" dialog should be disabled", () => shouldBe(page, el("first \"Radar\" card in \"Add Viewer\" dialog"), "disabled"));
    await session.step(47, "And first \"Sankey\" card in \"Add Viewer\" dialog should be disabled", () => shouldBe(page, el("first \"Sankey\" card in \"Add Viewer\" dialog"), "disabled"));
    await session.step(48, "And first \"Sunburst\" card in \"Add Viewer\" dialog should be enabled", () => shouldBe(page, el("first \"Sunburst\" card in \"Add Viewer\" dialog"), "enabled"));
    await session.step(49, "And first \"Tree\" card in \"Add Viewer\" dialog should be enabled", () => shouldBe(page, el("first \"Tree\" card in \"Add Viewer\" dialog"), "enabled"));
    await session.step(50, "When user hovers over first \"Radar\" card in \"Add Viewer\" dialog", () => hoverOver(page, el("first \"Radar\" card in \"Add Viewer\" dialog")));
    await session.step(51, "Then tooltip should contain text \"Radar viewer needs at least 1 numerical column\"", () => shouldContainText(page, el("tooltip"), "Radar viewer needs at least 1 numerical column"));
    await session.step(52, "When user hovers over first \"Sankey\" card in \"Add Viewer\" dialog", () => hoverOver(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(53, "Then tooltip should contain text \"Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column\"", () => shouldContainText(page, el("tooltip"), "Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column"));
    await session.step(54, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(55, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
  });
  test("Globe draws on earthquakes", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(59, "Given user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(60, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(61, "And user clicks on first \"Globe\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Globe\" card in \"Add Viewer\" dialog")));
    await session.step(62, "Then globe viewer should be visible", () => shouldBe(page, el("globe viewer"), "visible"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
    await session.step(64, "When user clicks on close icon of globe viewer", () => clickOn(page, el("close icon of globe viewer")));
    await session.step(65, "Then globe viewer should be absent", () => shouldBe(page, el("globe viewer"), "absent"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
  });
  test("Group Analysis adds an analysed column and keeps it in a layout (GROK-19039, GROK-19047)", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(69, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(70, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(71, "And user clicks on first \"Group Analysis\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Group Analysis\" card in \"Add Viewer\" dialog")));
    await session.step(72, "Then group analysis viewer should be visible", () => shouldBe(page, el("group analysis viewer"), "visible"));
    await session.step(73, "When user clicks on first grid", () => clickOn(page, el("first grid")));
    await session.step(74, "And user clicks on settings icon of group analysis viewer", () => clickOn(page, el("settings icon of group analysis viewer")));
    await session.step(75, "Then \"Group By\" property in context panel should be visible", () => shouldBe(page, el("\"Group By\" property in context panel"), "visible"));
    await session.step(76, "When user clicks on \"...\" button in \"Group By\" property in context panel", () => clickOn(page, el("\"...\" button in \"Group By\" property in context panel")));
    await session.step(77, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(78, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(79, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(80, "When user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(81, "Then \"1 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"1 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(82, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(83, "Then \"Group By\" property of group analysis viewer should be \"SEX\"", () => propertyShouldBe(page, "Group By", el("group analysis viewer"), "SEX"));
    await session.step(84, "And the \"rows\" reading of grid in group analysis viewer should be 2", () => readingIs(page, "rows", el("grid in group analysis viewer"), 2));
    await session.step(85, "And the \"text of cell 1 of SEX\" reading of grid in group analysis viewer should be \"F\"", () => readingReads(page, "text of cell 1 of SEX", el("grid in group analysis viewer"), "F"));
    await session.step(86, "And the \"text of cell 2 of SEX\" reading of grid in group analysis viewer should be \"M\"", () => readingReads(page, "text of cell 2 of SEX", el("grid in group analysis viewer"), "M"));
    await session.step(87, "And grid in group analysis viewer should not have a \"header min(AGE)\" area", () => hasNoArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(88, "When user clicks on \"Add column to analyze\" icon in group analysis viewer", () => clickOn(page, el("\"Add column to analyze\" icon in group analysis viewer")));
    await session.step(89, "Then \"Add column\" dialog should be visible", () => shouldBe(page, el("\"Add column\" dialog"), "visible"));
    await session.step(90, "When user selects \"AGE\" in Column input in \"Add column\" dialog", () => selectIn(page, "AGE", el("Column input in \"Add column\" dialog")));
    await session.step(91, "And user clicks on OK button in \"Add column\" dialog", () => clickOn(page, el("OK button in \"Add column\" dialog")));
    await session.step(92, "Then \"Add column\" dialog should be absent", () => shouldBe(page, el("\"Add column\" dialog"), "absent"));
    await session.step(93, "And grid in group analysis viewer should have a \"header min(AGE)\" area", () => hasArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(94, "And the \"text of cell 1 of min(AGE)\" reading of grid in group analysis viewer should be \"18.00\"", () => readingReads(page, "text of cell 1 of min(AGE)", el("grid in group analysis viewer"), "18.00"));
    await session.step(95, "And no errors should have been logged", () => noErrors(page));
    await session.step(96, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(97, "And user clicks on close icon of group analysis viewer", () => clickOn(page, el("close icon of group analysis viewer")));
    await session.step(98, "Then group analysis viewer should be absent", () => shouldBe(page, el("group analysis viewer"), "absent"));
    await session.step(99, "When user loads the saved layout", () => loadLayout(page));
    await session.step(100, "Then group analysis viewer should be visible", () => shouldBe(page, el("group analysis viewer"), "visible"));
    await session.step(101, "And \"Group By\" property of group analysis viewer should be \"SEX\"", () => propertyShouldBe(page, "Group By", el("group analysis viewer"), "SEX"));
    await session.step(102, "And the \"rows\" reading of grid in group analysis viewer should be 2", () => readingIs(page, "rows", el("grid in group analysis viewer"), 2));
    await session.step(103, "And grid in group analysis viewer should have a \"header min(AGE)\" area", () => hasArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(104, "And no errors should have been logged", () => noErrors(page));
  });
});
