/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/charts-gallery/charts-gallery.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.surface-plot, charts.viewer.globe, charts.viewer.group-analysis]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, followingShouldBe, hoverOver, isExpanded, pressKey, selectIn, shouldBe, shouldContainText, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, browsePanelOpen, layoutsDeleted, openDataset, openDatasetRowsAs, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hasArea, hasNoArea, noErrors, propertyShouldBe, propertyShouldNotBe, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Add viewer gallery for Charts viewers, Surface plot, Globe and Group Analysis", () => {
  const session = feature(test, "features/viewers/charts-gallery/charts-gallery.feature", import.meta.url);
  test("The gallery disables the Charts viewers a table cannot feed and says why", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(17, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(18, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(19, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(20, "And the following elements should be enabled:", () => followingShouldBe(page, "enabled", [["first \"Sankey\" card in \"Add Viewer\" dialog"],["first \"Chord\" card in \"Add Viewer\" dialog"],["first \"Timelines\" card in \"Add Viewer\" dialog"],["first \"Radar\" card in \"Add Viewer\" dialog"],["first \"Sunburst\" card in \"Add Viewer\" dialog"],["first \"Tree\" card in \"Add Viewer\" dialog"],["first \"Surface plot\" card in \"Add Viewer\" dialog"],["first \"Globe\" card in \"Add Viewer\" dialog"],["first \"Group Analysis\" card in \"Add Viewer\" dialog"],["first \"Word cloud\" card in \"Add Viewer\" dialog"]]), [["first \"Sankey\" card in \"Add Viewer\" dialog"],["first \"Chord\" card in \"Add Viewer\" dialog"],["first \"Timelines\" card in \"Add Viewer\" dialog"],["first \"Radar\" card in \"Add Viewer\" dialog"],["first \"Sunburst\" card in \"Add Viewer\" dialog"],["first \"Tree\" card in \"Add Viewer\" dialog"],["first \"Surface plot\" card in \"Add Viewer\" dialog"],["first \"Globe\" card in \"Add Viewer\" dialog"],["first \"Group Analysis\" card in \"Add Viewer\" dialog"],["first \"Word cloud\" card in \"Add Viewer\" dialog"]]);
    await session.step(31, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(32, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(33, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(36, "And Files---App-Data---Chem tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Chem tree node inside browse tree")));
    await session.step(37, "When user double-clicks on Files---App-Data---Chem---chem_standards.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Chem---chem_standards.csv tree node inside browse tree")));
    await session.step(38, "Then the \"chem_standards\" view should be current", () => viewIsCurrent(page, "chem_standards"));
    await session.step(39, "And the table should have 2 columns", () => columnCount(page, 2));
    await session.step(40, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(41, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(42, "And first \"Radar\" card in \"Add Viewer\" dialog should be disabled", () => shouldBe(page, el("first \"Radar\" card in \"Add Viewer\" dialog"), "disabled"));
    await session.step(43, "And first \"Sankey\" card in \"Add Viewer\" dialog should be disabled", () => shouldBe(page, el("first \"Sankey\" card in \"Add Viewer\" dialog"), "disabled"));
    await session.step(44, "And first \"Sunburst\" card in \"Add Viewer\" dialog should be enabled", () => shouldBe(page, el("first \"Sunburst\" card in \"Add Viewer\" dialog"), "enabled"));
    await session.step(45, "And first \"Tree\" card in \"Add Viewer\" dialog should be enabled", () => shouldBe(page, el("first \"Tree\" card in \"Add Viewer\" dialog"), "enabled"));
    await session.step(46, "When user hovers over first \"Radar\" card in \"Add Viewer\" dialog", () => hoverOver(page, el("first \"Radar\" card in \"Add Viewer\" dialog")));
    await session.step(47, "Then tooltip should contain text \"Radar viewer needs at least 1 numerical column\"", () => shouldContainText(page, el("tooltip"), "Radar viewer needs at least 1 numerical column"));
    await session.step(48, "When user hovers over first \"Sankey\" card in \"Add Viewer\" dialog", () => hoverOver(page, el("first \"Sankey\" card in \"Add Viewer\" dialog")));
    await session.step(49, "Then tooltip should contain text \"Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column\"", () => shouldContainText(page, el("tooltip"), "Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column"));
    await session.step(50, "When user presses Escape", () => pressKey(page, "Escape"));
    await session.step(51, "Then \"Add Viewer\" dialog should be absent", () => shouldBe(page, el("\"Add Viewer\" dialog"), "absent"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("Globe draws on earthquakes", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(55, "Given user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(56, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(57, "And user clicks on first \"Globe\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Globe\" card in \"Add Viewer\" dialog")));
    await session.step(58, "Then globe viewer should be visible", () => shouldBe(page, el("globe viewer"), "visible"));
    await session.step(59, "And the \"points\" reading of globe viewer should be 2426", () => readingIs(page, "points", el("globe viewer"), 2426));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
    await session.step(61, "When user clicks on close icon of globe viewer", () => clickOn(page, el("close icon of globe viewer")));
    await session.step(62, "Then globe viewer should be absent", () => shouldBe(page, el("globe viewer"), "absent"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
  test("Surface plot on demog takes its columns, Projection and Wireframe", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(66, "Given user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(67, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(68, "And user clicks on first \"Surface plot\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Surface plot\" card in \"Add Viewer\" dialog")));
    await session.step(69, "Then surface plot viewer should be visible", () => shouldBe(page, el("surface plot viewer"), "visible"));
    await session.step(70, "And \"XColumnName\" property of surface plot viewer should not be \"\"", () => propertyShouldNotBe(page, "XColumnName", el("surface plot viewer"), ""));
    await session.step(71, "And \"YColumnName\" property of surface plot viewer should not be \"\"", () => propertyShouldNotBe(page, "YColumnName", el("surface plot viewer"), ""));
    await session.step(72, "And \"ZColumnName\" property of surface plot viewer should not be \"\"", () => propertyShouldNotBe(page, "ZColumnName", el("surface plot viewer"), ""));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
    await session.step(74, "When user clicks on grid", () => clickOn(page, el("grid")));
    await session.step(75, "And user clicks on settings icon of surface plot viewer", () => clickOn(page, el("settings icon of surface plot viewer")));
    await session.step(76, "Given \"Misc\" category in context panel is expanded", () => isExpanded(page, el("\"Misc\" category in context panel")));
    await session.step(77, "Then \"Projection\" property in context panel should be visible", () => shouldBe(page, el("\"Projection\" property in context panel"), "visible"));
    await session.step(78, "When user selects \"orthographic\" in \"Projection\" property in context panel", () => selectIn(page, "orthographic", el("\"Projection\" property in context panel")));
    await session.step(79, "Then \"Projection\" property of surface plot viewer should be \"orthographic\"", () => propertyShouldBe(page, "Projection", el("surface plot viewer"), "orthographic"));
    await session.step(80, "When user unchecks \"Wireframe\" property in context panel", () => uncheck(page, el("\"Wireframe\" property in context panel")));
    await session.step(81, "Then \"Wireframe\" property of surface plot viewer should be \"false\"", () => propertyShouldBe(page, "Wireframe", el("surface plot viewer"), "false"));
    await session.step(82, "And surface plot viewer should be visible", () => shouldBe(page, el("surface plot viewer"), "visible"));
    await session.step(83, "And no errors should have been logged", () => noErrors(page));
  });
  test("Group Analysis adds an analysed column and keeps it in a layout (GROK-19039, GROK-19047)", {tag: ["@viewers", "@realizes:charts.viewer.surface-plot", "@realizes:charts.viewer.globe", "@realizes:charts.viewer.group-analysis"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(86, "Given user opens demog dataset keeping the first 5850 rows as \"GroupAnalysis-{time}\"", () => openDatasetRowsAs(page, ds("demog"), 5850, session.text("GroupAnalysis-{time}")));
    await session.step(87, "And the layouts named \"GroupAnalysis-{time}\" are deleted when the feature ends", () => layoutsDeleted(page, session.text("GroupAnalysis-{time}")));
    await session.step(88, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(89, "And user clicks on first \"Group Analysis\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Group Analysis\" card in \"Add Viewer\" dialog")));
    await session.step(90, "Then group analysis viewer should be visible", () => shouldBe(page, el("group analysis viewer"), "visible"));
    await session.step(91, "When user clicks on first grid", () => clickOn(page, el("first grid")));
    await session.step(92, "And user clicks on settings icon of group analysis viewer", () => clickOn(page, el("settings icon of group analysis viewer")));
    await session.step(93, "Then \"Group By\" property in context panel should be visible", () => shouldBe(page, el("\"Group By\" property in context panel"), "visible"));
    await session.step(94, "When user clicks on \"...\" button in \"Group By\" property in context panel", () => clickOn(page, el("\"...\" button in \"Group By\" property in context panel")));
    await session.step(95, "Then \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"Select columns...\" dialog"), "visible"));
    await session.step(96, "When user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(97, "Then \"0 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"0 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(98, "When user toggles the \"SEX\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "SEX", el("\"Select columns...\" dialog")));
    await session.step(99, "Then \"1 checked\" text in \"Select columns...\" dialog should be visible", () => shouldBe(page, el("\"1 checked\" text in \"Select columns...\" dialog"), "visible"));
    await session.step(100, "When user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(101, "Then \"Group By\" property of group analysis viewer should be \"SEX\"", () => propertyShouldBe(page, "Group By", el("group analysis viewer"), "SEX"));
    await session.step(102, "And the \"rows\" reading of grid in group analysis viewer should be 2", () => readingIs(page, "rows", el("grid in group analysis viewer"), 2));
    await session.step(103, "And the \"text of cell 1 of SEX\" reading of grid in group analysis viewer should be \"F\"", () => readingReads(page, "text of cell 1 of SEX", el("grid in group analysis viewer"), "F"));
    await session.step(104, "And the \"text of cell 2 of SEX\" reading of grid in group analysis viewer should be \"M\"", () => readingReads(page, "text of cell 2 of SEX", el("grid in group analysis viewer"), "M"));
    await session.step(105, "And grid in group analysis viewer should not have a \"header min(AGE)\" area", () => hasNoArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(106, "When user clicks on \"Add column to analyze\" icon in group analysis viewer", () => clickOn(page, el("\"Add column to analyze\" icon in group analysis viewer")));
    await session.step(107, "Then \"Add column\" dialog should be visible", () => shouldBe(page, el("\"Add column\" dialog"), "visible"));
    await session.step(108, "When user selects \"AGE\" in Column input in \"Add column\" dialog", () => selectIn(page, "AGE", el("Column input in \"Add column\" dialog")));
    await session.step(109, "And user clicks on OK button in \"Add column\" dialog", () => clickOn(page, el("OK button in \"Add column\" dialog")));
    await session.step(110, "Then \"Add column\" dialog should be absent", () => shouldBe(page, el("\"Add column\" dialog"), "absent"));
    await session.step(111, "And grid in group analysis viewer should have a \"header min(AGE)\" area", () => hasArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(112, "And the \"text of cell 1 of min(AGE)\" reading of grid in group analysis viewer should be \"18.00\"", () => readingReads(page, "text of cell 1 of min(AGE)", el("grid in group analysis viewer"), "18.00"));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
    await session.step(114, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(115, "And Layouts accordion header in toolbox is expanded", () => isExpanded(page, el("Layouts accordion header in toolbox")));
    await session.step(116, "When user clicks on Save button in layouts pane", () => clickOn(page, el("Save button in layouts pane")));
    await session.step(117, "Then \"GroupAnalysis-{time}\" layout card should be visible", () => shouldBe(page, el(session.text("\"GroupAnalysis-{time}\" layout card")), "visible"));
    await session.step(118, "When user clicks on close icon of group analysis viewer", () => clickOn(page, el("close icon of group analysis viewer")));
    await session.step(119, "Then group analysis viewer should be absent", () => shouldBe(page, el("group analysis viewer"), "absent"));
    await session.step(120, "When user clicks on \"GroupAnalysis-{time}\" layout card", () => clickOn(page, el(session.text("\"GroupAnalysis-{time}\" layout card"))));
    await session.step(121, "Then group analysis viewer should be visible", () => shouldBe(page, el("group analysis viewer"), "visible"));
    await session.step(122, "And \"Group By\" property of group analysis viewer should be \"SEX\"", () => propertyShouldBe(page, "Group By", el("group analysis viewer"), "SEX"));
    await session.step(123, "And the \"rows\" reading of grid in group analysis viewer should be 2", () => readingIs(page, "rows", el("grid in group analysis viewer"), 2));
    await session.step(124, "And grid in group analysis viewer should have a \"header min(AGE)\" area", () => hasArea(page, el("grid in group analysis viewer"), "header min(AGE)"));
    await session.step(125, "And no errors should have been logged", () => noErrors(page));
  });
});
