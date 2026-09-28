/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-join.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20013]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/minimized-viewers.js';
import '../../bindings/projects-copies.js';
import '../../bindings/projects-derived.js';
import '../../bindings/projects-sources.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noErrorsButPreviewNoise, tableViewTabs} from '../../bindings/projects-regressions.js';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {taskBarFinished, watchTaskBar} from '@datagrok-libraries/bdd/bindings/platform/events';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, savedWithDataSync, switchTableView, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTablesExactly, tableViewsOpen} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a copy saved after an in-place join", () => {
  const session = feature(test, "features/projects/projects-regressions-join.feature", import.meta.url);
  test("Projects regressions: a copy saved after an in-place join", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:GROK-20013"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("A copy is saved after an in-place join onto a reopened join result", async () => {
      await session.step(29, "Given no project named \"BDDRegJoin{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegJoin{time}")));
      await session.step(30, "And no project named \"BDDRegJoinCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegJoinCopy{time}")));
      await session.step(31, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(32, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(33, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(34, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(35, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(36, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(37, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(38, "Then the \"demog (2)\" view should be current", () => viewIsCurrent(page, "demog (2)"));
      await session.step(39, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
      await session.step(40, "Then \"Join Tables\" dialog should be visible", () => shouldBe(page, el("\"Join Tables\" dialog"), "visible"));
      await session.step(41, "And join left table selector should have value \"demog\"", () => shouldHaveValue(page, el("join left table selector"), "demog"));
      await session.step(42, "And join right table selector should have value \"demog (2)\"", () => shouldHaveValue(page, el("join right table selector"), "demog (2)"));
      await session.step(43, "When user clicks on OK button in \"Join Tables\" dialog", () => clickOn(page, el("OK button in \"Join Tables\" dialog")));
      await session.step(44, "Then the \"Join Tables\" dialog should close", () => dialogCloses(page, "Join Tables"));
      await session.step(45, "And the \"result\" view should be current", () => viewIsCurrent(page, "result"));
      await session.step(46, "And table \"result\" should have 5850 rows", () => tableRows(page, "result", 5850));
      await session.step(47, "And the table should have 22 columns", () => columnCount(page, 22));
      await session.step(48, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(49, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(50, "When user enters \"BDDRegJoin{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegJoin{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(51, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(52, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(53, "And \"Share BDDRegJoin{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegJoin{time}\" dialog")), "visible"));
      await session.step(54, "When user clicks on CANCEL button in \"Share BDDRegJoin{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegJoin{time}\" dialog"))));
      await session.step(55, "Then the \"Share BDDRegJoin{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegJoin{time}")));
      await session.step(56, "And the \"result\" table of the \"BDDRegJoin{time}\" project should be saved with data sync", () => savedWithDataSync(page, "result", session.text("BDDRegJoin{time}")));
      await session.step(57, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(58, "When user closes all views", () => closeAllViews(page));
      await session.step(59, "And the browse panel is open", () => browsePanelOpen(page));
      await session.step(60, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(61, "And user enters \"BDDRegJoin{time}\" into gallery search", () => enterInto(page, session.text("BDDRegJoin{time}"), el("gallery search")));
      await session.step(62, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(63, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(64, "When user double-clicks on BDDRegJoin{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegJoin{time} gallery card"))));
      await session.step(65, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
      await session.step(66, "Then the table views \"demog, demog (2), result\" should be open", () => tableViewsOpen(page, "demog, demog (2), result"));
      await session.step(67, "And no errors should have been logged", () => noErrors(page));
      await session.step(68, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(69, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
      await session.step(70, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(71, "When user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(72, "And user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
      await session.step(73, "Then \"Join Tables\" dialog should be visible", () => shouldBe(page, el("\"Join Tables\" dialog"), "visible"));
      await session.step(74, "When user selects \"result\" in join left table selector", () => selectIn(page, "result", el("join left table selector")));
      await session.step(75, "And user selects \"demog-1000\" in join right table selector", () => selectIn(page, "demog-1000", el("join right table selector")));
      await session.step(76, "And user selects \"left\" in \"Join Type\" input in \"Join Tables\" dialog", () => selectIn(page, "left", el("\"Join Type\" input in \"Join Tables\" dialog")));
      await session.step(77, "And user checks \"In-place\" input in \"Join Tables\" dialog", () => check(page, el("\"In-place\" input in \"Join Tables\" dialog")));
      await session.step(78, "And user clicks on OK button in \"Join Tables\" dialog", () => clickOn(page, el("OK button in \"Join Tables\" dialog")));
      await session.step(79, "Then the \"Join Tables\" dialog should close", () => dialogCloses(page, "Join Tables"));
      await session.step(80, "And the open tables should be exactly \"demog, demog (2), demog-1000, result+demog-1000\"", () => openTablesExactly(page, "demog, demog (2), demog-1000, result+demog-1000"));
      await session.step(81, "And table \"result+demog-1000\" should have 5850 rows", () => tableRows(page, "result+demog-1000", 5850));
      await session.step(82, "When user switches to the \"result+demog-1000\" table view", () => switchTableView(page, "result+demog-1000"));
      await session.step(83, "Then the table should have 33 columns", () => columnCount(page, 33));
      await session.step(84, "And the table should have a column \"demog-1000.AGE\"", () => hasColumn(page, "demog-1000.AGE"));
      await session.step(85, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(86, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(87, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(88, "And user enters \"BDDRegJoinCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegJoinCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(89, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(90, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(91, "And an info balloon containing 'Project \"BDDRegJoinCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDRegJoinCopy{time}\" uploaded")));
      await session.step(92, "And the \"result\" table of the \"BDDRegJoinCopy{time}\" project should be saved with data sync", () => savedWithDataSync(page, "result", session.text("BDDRegJoinCopy{time}")));
      await session.step(93, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(94, "When user picks \"Close All\" from the context menu of left sidebar", () => pickFromContextMenu(page, "Close All", el("left sidebar")));
      await session.step(95, "Then the table view tabs should read \"\"", () => tableViewTabs(page, ""));
    });
    await run.scenario("The copy reopens in a reloaded page with all its tables, the join done in place", async () => {
      await session.step(99, "When user reloads the page", () => reloadPage(page));
      await session.step(100, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(101, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(102, "And user enters \"BDDRegJoinCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDRegJoinCopy{time}"), el("gallery search")));
      await session.step(103, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(104, "Given user watches the task bar", () => watchTaskBar(page));
      await session.step(105, "When user double-clicks on BDDRegJoinCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegJoinCopy{time} gallery card"))));
      await session.step(106, "Then the task bar should have finished \"Opening project\"", () => taskBarFinished(page, "Opening project"));
      await session.step(107, "Then the open tables should be exactly \"demog, demog (2), demog-1000, result+demog-1000\"", () => openTablesExactly(page, "demog, demog (2), demog-1000, result+demog-1000"));
      await session.step(108, "And table \"result+demog-1000\" should have 5850 rows", () => tableRows(page, "result+demog-1000", 5850));
      await session.step(109, "When user switches to the \"result+demog-1000\" table view", () => switchTableView(page, "result+demog-1000"));
      await session.step(110, "Then the table should have 33 columns", () => columnCount(page, 33));
      await session.step(111, "And the table should have a column \"demog-1000.AGE\"", () => hasColumn(page, "demog-1000.AGE"));
      await session.step(112, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(113, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
