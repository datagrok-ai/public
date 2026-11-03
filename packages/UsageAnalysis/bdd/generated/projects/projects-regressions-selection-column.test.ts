/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-selection-column.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#602]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {selectFirstRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, dialogCloses, doubleClickOpensProject, loadedAsSnapshot, noProjectOnServer, reloadedByDataSync, savedAsSnapshot, savedWithDataSync, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromOpenMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a column saved from the row selection", () => {
  const session = feature(test, "features/projects/projects-regressions-selection-column.feature", import.meta.url);
  test("Saved as a snapshot, the table comes back with the selection column", {tag: ["@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(18, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(19, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(20, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(21, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(22, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(23, "When user selects the first 10 rows", () => selectFirstRows(page, 10));
    await session.step(24, "And user clicks on \"Selected: 10\" text in status bar", () => clickOn(page, el("\"Selected: 10\" text in status bar")));
    await session.step(25, "And user picks \"Save as Column\" from the open menu", () => pickFromOpenMenu(page, "Save as Column"));
    await session.step(26, "Then \"A\" column should have type \"bool\"", () => columnType(page, "A", "bool"));
    await session.step(29, "Given no project named \"BDDRegSelColSnap{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegSelColSnap{time}")));
    await session.step(30, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(31, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(32, "When user enters \"BDDRegSelColSnap{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegSelColSnap{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(33, "And user unchecks Data sync switch in \"demog-1000\" project table in \"Save project\" dialog", () => uncheck(page, el("Data sync switch in \"demog-1000\" project table in \"Save project\" dialog")));
    await session.step(34, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(35, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(36, "And \"Share BDDRegSelColSnap{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegSelColSnap{time}\" dialog")), "visible"));
    await session.step(37, "When user clicks on CANCEL button in \"Share BDDRegSelColSnap{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegSelColSnap{time}\" dialog"))));
    await session.step(38, "Then the \"Share BDDRegSelColSnap{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegSelColSnap{time}")));
    await session.step(39, "And the \"demog-1000\" table of the \"BDDRegSelColSnap{time}\" project should be saved as a snapshot", () => savedAsSnapshot(page, "demog-1000", session.text("BDDRegSelColSnap{time}")));
    await session.step(40, "When user closes all views", () => closeAllViews(page));
    await session.step(41, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(42, "And user enters \"BDDRegSelColSnap{time}\" into gallery search", () => enterInto(page, session.text("BDDRegSelColSnap{time}"), el("gallery search")));
    await session.step(43, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(44, "And user double-clicks on BDDRegSelColSnap{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDRegSelColSnap{time} gallery card"))));
    await session.step(45, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(46, "And the table should have been loaded as a snapshot", () => loadedAsSnapshot(page));
    await session.step(47, "And \"A\" column should have type \"bool\"", () => columnType(page, "A", "bool"));
  });
  test("Saved with Data sync, the table comes back with the selection column", {tag: ["@realizes:views.projects", "@known-failure:70", "@realizes:core#602"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(18, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(19, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(20, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(21, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(22, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(23, "When user selects the first 10 rows", () => selectFirstRows(page, 10));
    await session.step(24, "And user clicks on \"Selected: 10\" text in status bar", () => clickOn(page, el("\"Selected: 10\" text in status bar")));
    await session.step(25, "And user picks \"Save as Column\" from the open menu", () => pickFromOpenMenu(page, "Save as Column"));
    await session.step(26, "Then \"A\" column should have type \"bool\"", () => columnType(page, "A", "bool"));
    await knownFailure(async () => {
      await session.step(52, "Given no project named \"BDDRegSelColSync{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegSelColSync{time}")));
      await session.step(53, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(54, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(55, "And Data sync switch in \"demog-1000\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog-1000\" project table in \"Save project\" dialog"), "checked"));
      await session.step(56, "When user enters \"BDDRegSelColSync{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegSelColSync{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(57, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(58, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(59, "And \"Share BDDRegSelColSync{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegSelColSync{time}\" dialog")), "visible"));
      await session.step(60, "When user clicks on CANCEL button in \"Share BDDRegSelColSync{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegSelColSync{time}\" dialog"))));
      await session.step(61, "Then the \"Share BDDRegSelColSync{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegSelColSync{time}")));
      await session.step(62, "And the \"demog-1000\" table of the \"BDDRegSelColSync{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog-1000", session.text("BDDRegSelColSync{time}")));
      await session.step(63, "When user closes all views", () => closeAllViews(page));
      await session.step(64, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(65, "And user enters \"BDDRegSelColSync{time}\" into gallery search", () => enterInto(page, session.text("BDDRegSelColSync{time}"), el("gallery search")));
      await session.step(66, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(67, "And user double-clicks on BDDRegSelColSync{time} gallery card and the project opens", () => doubleClickOpensProject(page, el(session.text("BDDRegSelColSync{time} gallery card"))));
      await session.step(68, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(69, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
      await session.step(70, "And \"A\" column should have type \"bool\"", () => columnType(page, "A", "bool"));
    }, 70);
  });
});
