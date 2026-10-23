/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/analyses.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [jmp.analyses]
--- */
import {test} from '@playwright/test';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {above} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, noTableViewOpen, tableViewOpened} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, noViewerError, propertyShouldBe, viewerCount, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Saved JMP analyses open as Datagrok viewers, and row states come across", () => {
  const session = feature(test, "features/analyses.feature", import.meta.url);
  test("A saved Graph Builder opens as a bar chart with the table", {tag: ["@realizes:jmp.analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "When user drops the \"fixtures/saved-bar-chart.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/saved-bar-chart.jmp", el("browse tree")));
    await session.step(24, "Then the \"saved-bar-chart\" table view should open with 130 rows", () => tableViewOpened(page, "saved-bar-chart", 130));
    await session.step(25, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
    await session.step(26, "And \"Split\" property of bar chart viewer should be \"A\"", () => propertyShouldBe(page, "Split", el("bar chart viewer"), "A"));
    await session.step(27, "And \"Title\" property of bar chart viewer should be \"Bar chart of A\"", () => propertyShouldBe(page, "Title", el("bar chart viewer"), "Bar chart of A"));
    await session.step(28, "And no viewer of the current view should report an error", () => noViewerError(page));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
    await session.step(30, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An analysis with no Datagrok counterpart is named, and the table still opens", {tag: ["@realizes:jmp.analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(34, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(35, "When user drops the \"fixtures/saved-fit-model.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/saved-fit-model.jmp", el("browse tree")));
    await session.step(36, "Then the \"saved-fit-model\" table view should open with 130 rows", () => tableViewOpened(page, "saved-fit-model", 130));
    await session.step(37, "And the open tableview should have 0 bar chart viewers", () => viewerCount(page, 0, "bar chart"));
    await session.step(38, "And a warning balloon containing \"1 saved analysis not converted: Fit A\" should have been shown", () => warningBalloonText(page, "1 saved analysis not converted: Fit A"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
  });
  test("The preview lists the saved analyses, and Open opens them", {tag: ["@realizes:jmp.analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(42, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(43, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(44, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(45, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(46, "And Files---App-Data---Jmp tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Jmp tree node inside browse tree")));
    await session.step(47, "And Files---App-Data---Jmp---tests tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Jmp---tests tree node inside browse tree")));
    await session.step(48, "When user double-clicks on Files---App-Data---Jmp---tests---saved-bar-chart.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Jmp---tests---saved-bar-chart.jmp tree node inside browse tree")));
    await session.step(49, "Then the \"saved-bar-chart.jmp\" view should be current", () => viewIsCurrent(page, "saved-bar-chart.jmp"));
    await session.step(50, "And \"Saved analyses\" section should be visible", () => shouldBe(page, el("\"Saved analyses\" section"), "visible"));
    await session.step(51, "And \"Saved analyses\" section should contain text \"1 of 1 open as viewers\"", () => shouldContainText(page, el("\"Saved analyses\" section"), "1 of 1 open as viewers"));
    await session.step(52, "And \"Saved analyses\" section should contain text \"Bar chart of A\"", () => shouldContainText(page, el("\"Saved analyses\" section"), "Bar chart of A"));
    await session.step(53, "And \"Saved analyses\" section should be above grid", () => above(page, el("\"Saved analyses\" section"), el("grid")));
    await session.step(54, "When user clicks on Open button", () => clickOn(page, el("Open button")));
    await session.step(55, "Then the \"saved-bar-chart\" table view should open with 130 rows", () => tableViewOpened(page, "saved-bar-chart", 130));
    await session.step(56, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
  });
  test("The preview says why an analysis doesn't open", {tag: ["@realizes:jmp.analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(60, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(61, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(62, "And Files---App-Data tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data tree node inside browse tree")));
    await session.step(63, "And Files---App-Data---Jmp tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Jmp tree node inside browse tree")));
    await session.step(64, "And Files---App-Data---Jmp---tests tree node inside browse tree is expanded", () => isExpanded(page, el("Files---App-Data---Jmp---tests tree node inside browse tree")));
    await session.step(65, "When user double-clicks on Files---App-Data---Jmp---tests---saved-fit-model.jmp tree node inside browse tree", () => doubleClickOn(page, el("Files---App-Data---Jmp---tests---saved-fit-model.jmp tree node inside browse tree")));
    await session.step(66, "Then the \"saved-fit-model.jmp\" view should be current", () => viewIsCurrent(page, "saved-fit-model.jmp"));
    await session.step(67, "And \"Saved analyses\" section should contain text \"0 of 1 open as viewers\"", () => shouldContainText(page, el("\"Saved analyses\" section"), "0 of 1 open as viewers"));
    await session.step(68, "And \"Saved analyses\" section should contain text \"Not converted: JMP's Fit Model isn't converted\"", () => shouldContainText(page, el("\"Saved analyses\" section"), "Not converted: JMP's Fit Model isn't converted"));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
  test("Row states come across as marker and color columns", {tag: ["@realizes:jmp.analyses"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(72, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(73, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(74, "When user drops the \"fixtures/rowstate.jmp\" file of the project onto browse tree", () => dropFile(page, "fixtures/rowstate.jmp", el("browse tree")));
    await session.step(75, "Then the \"rowstate\" table view should open with 3 rows", () => tableViewOpened(page, "rowstate", 3));
    await session.step(76, "And the table should have the columns \"rowstate1 marker, rowstate1 color, rowstate2 marker, rowstate2 color, rowstate3 marker, rowstate3 color\"", () => columnsExactly(page, "rowstate1 marker, rowstate1 color, rowstate2 marker, rowstate2 color, rowstate3 marker, rowstate3 color"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
    await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
