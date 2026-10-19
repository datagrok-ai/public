/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-join.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, data.menu.join-tables, GROK-20013]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/file-formats.js';
import '../../bindings/flow.js';
import '../../bindings/grid.js';
import '../../bindings/queries.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn, reloadPage} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, doubleClickOn, enterInto, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, dialogCloses, doubleClickOpensProject, noProjectOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noTableLeft, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {infoBalloonText, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a copy saved after an in-place join", () => {
  const session = feature(test, "features/projects/projects-regressions-join.feature", import.meta.url);
  test("Projects regressions: a copy saved after an in-place join", {tag: ["@journey", "@realizes:views.projects", "@realizes:data.menu.join-tables", "@realizes:GROK-20013", "@known-failure:109"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And no project named \"BDDRegJoin{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegJoin{time}")));
    await session.step(31, "And no project named \"BDDRegJoinCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegJoinCopy{time}")));
    await run.scenario("A copy is saved after an in-place join onto a reopened join result", async () => {
      await session.step(35, "Given simple mode is off", () => simpleModeOff(page));
      await session.step(36, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(37, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(38, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(39, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(40, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
      await session.step(41, "And user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(42, "Then the \"demog (2)\" view should be current", () => viewIsCurrent(page, "demog (2)"));
      await session.step(43, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
      await session.step(44, "Then \"Join Tables\" dialog should be visible", () => shouldBe(page, el("\"Join Tables\" dialog"), "visible"));
      await session.step(45, "When user selects \"demog\" in Tables input in \"Join Tables\" dialog", () => selectIn(page, "demog", el("Tables input in \"Join Tables\" dialog")));
      await session.step(46, "And user selects \"demog (2)\" in second choice input in \"Join Tables\" dialog", () => selectIn(page, "demog (2)", el("second choice input in \"Join Tables\" dialog")));
      await session.step(47, "Then \"Preview result columns\" dialog should be visible", () => shouldBe(page, el("\"Preview result columns\" dialog"), "visible"));
      await session.step(48, "When user clicks on OK button in \"Join Tables\" dialog", () => clickOn(page, el("OK button in \"Join Tables\" dialog")));
      await session.step(49, "Then the \"Join Tables\" dialog should close", () => dialogCloses(page, "Join Tables"));
      await session.step(50, "And the \"result\" view should be current", () => viewIsCurrent(page, "result"));
      await session.step(51, "And table \"result\" should have 5850 rows", () => tableRows(page, "result", 5850));
      await session.step(52, "And the table should have 22 columns", () => columnCount(page, 22));
      await session.step(53, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(54, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(55, "When user enters \"BDDRegJoin{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegJoin{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(56, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(57, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(58, "And \"Share BDDRegJoin{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegJoin{time}\" dialog")), "visible"));
      await session.step(59, "When user clicks on CANCEL button in \"Share BDDRegJoin{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegJoin{time}\" dialog"))));
      await session.step(60, "Then the \"Share BDDRegJoin{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegJoin{time}")));
      await session.step(61, "When user closes all views", () => closeAllViews(page));
      await session.step(62, "Then no table should be left in the workspace", () => noTableLeft(page));
      await session.step(63, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(64, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(65, "And user enters \"BDDRegJoin{time}\" into gallery search", () => enterInto(page, session.text("BDDRegJoin{time}"), el("gallery search")));
      await session.step(66, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(67, "And user double-clicks on BDDRegJoin{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDRegJoin{time} gallery card"))));
      await session.step(68, "Then the open table views should be exactly \"demog (2), demog, result\"", () => openTableViewsExactly(page, "demog (2), demog, result"));
      await session.step(69, "And no errors should have been logged", () => noErrors(page));
      await session.step(70, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(71, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(72, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(73, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
      await session.step(74, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(75, "When user clicks on result view", () => clickOn(page, el("result view")));
      await session.step(76, "Then the \"result\" view should be current", () => viewIsCurrent(page, "result"));
      await session.step(77, "When user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
      await session.step(78, "Then \"Join Tables\" dialog should be visible", () => shouldBe(page, el("\"Join Tables\" dialog"), "visible"));
      await session.step(79, "When user selects \"result\" in Tables input in \"Join Tables\" dialog", () => selectIn(page, "result", el("Tables input in \"Join Tables\" dialog")));
      await session.step(80, "And user selects \"demog-1000\" in second choice input in \"Join Tables\" dialog", () => selectIn(page, "demog-1000", el("second choice input in \"Join Tables\" dialog")));
      await session.step(81, "Then \"Preview result columns\" dialog should be visible", () => shouldBe(page, el("\"Preview result columns\" dialog"), "visible"));
      await session.step(82, "When user selects \"left\" in \"Join Type\" input in \"Join Tables\" dialog", () => selectIn(page, "left", el("\"Join Type\" input in \"Join Tables\" dialog")));
      await session.step(83, "And user checks \"In-place\" input in \"Join Tables\" dialog", () => check(page, el("\"In-place\" input in \"Join Tables\" dialog")));
      await session.step(84, "And user clicks on OK button in \"Join Tables\" dialog", () => clickOn(page, el("OK button in \"Join Tables\" dialog")));
      await session.step(85, "Then the \"Join Tables\" dialog should close", () => dialogCloses(page, "Join Tables"));
      await session.step(86, "And the open table views should be exactly \"demog (2), demog, result+demog-1000, demog-1000\"", () => openTableViewsExactly(page, "demog (2), demog, result+demog-1000, demog-1000"));
      await session.step(87, "And table \"result+demog-1000\" should have 5850 rows", () => tableRows(page, "result+demog-1000", 5850));
      await session.step(88, "When user clicks on result+demog-1000 view", () => clickOn(page, el("result+demog-1000 view")));
      await session.step(89, "Then the table should have 33 columns", () => columnCount(page, 33));
      await session.step(90, "And the table should have a column \"demog-1000.AGE\"", () => hasColumn(page, "demog-1000.AGE"));
      await session.step(91, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(92, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(93, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(94, "And user enters \"BDDRegJoinCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegJoinCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(95, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(96, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(97, "And an info balloon containing 'Project \"BDDRegJoinCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDRegJoinCopy{time}\" uploaded")));
    });
    await run.scenario("The copy reopened in a new session comes back with the join applied once", async () => {
      await session.step(102, "When user closes all views", () => closeAllViews(page));
      await session.step(103, "And user reloads the page", () => reloadPage(page));
      await session.step(104, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(105, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(106, "And user enters \"BDDRegJoinCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDRegJoinCopy{time}"), el("gallery search")));
      await session.step(107, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(108, "And user double-clicks on BDDRegJoinCopy{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDRegJoinCopy{time} gallery card"))));
      await session.step(109, "Then the open table views should be exactly \"demog (2), demog, result+demog-1000, demog-1000\"", () => openTableViewsExactly(page, "demog (2), demog, result+demog-1000, demog-1000"));
      await session.step(110, "And table \"result+demog-1000\" should have 5850 rows", () => tableRows(page, "result+demog-1000", 5850));
      await session.step(111, "And \"Data loading error\" dialog should be absent", () => shouldBe(page, el("\"Data loading error\" dialog"), "absent"));
    }, {knownFailure: true, line: 109});
    run.finish();
  });
});
