/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-pinned.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-20607]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/nx.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, reloadedByDataSync, savedWithDataSync, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, pickFromAreaContextMenu, readingIs, readingReads, warningBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {readingExcludes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a project with pinned rows", () => {
  const session = feature(test, "features/projects/projects-regressions-pinned.feature", import.meta.url);
  test("The project with pinned rows reopens without errors, sorted and without HEIGHT", {tag: ["@serial", "@realizes:views.projects", "@realizes:GROK-20607"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "And no project named \"BDDRegPinned{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegPinned{time}")));
    await session.step(24, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(25, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(26, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(27, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(28, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(29, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(30, "When user picks \"Sort > Ascending\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Sort > Ascending", "header AGE", el("grid")));
    await session.step(31, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(32, "When user picks \"Hide\" from the context menu of the \"header HEIGHT\" area of grid", () => pickFromAreaContextMenu(page, "Hide", "header HEIGHT", el("grid")));
    await session.step(33, "Then the \"column order\" reading of grid should not include the text \"HEIGHT\"", () => readingExcludes(page, "column order", el("grid"), "HEIGHT"));
    await session.step(34, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 199 of SEX\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 199 of SEX", el("grid")));
    await session.step(35, "Then the \"pinned rows\" reading of grid should be 1", () => readingIs(page, "pinned rows", el("grid"), 1));
    await session.step(36, "When user picks \"Pin > Pin Row\" from the context menu of the \"cell 344 of RACE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Row", "cell 344 of RACE", el("grid")));
    await session.step(37, "Then the \"pinned rows\" reading of grid should be 2", () => readingIs(page, "pinned rows", el("grid"), 2));
    await session.step(38, "And a warning balloon containing \"pinned a non-unique value\" should have been shown", () => warningBalloonText(page, "pinned a non-unique value"));
    await session.step(39, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(40, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(41, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
    await session.step(42, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
    await session.step(43, "When user enters \"BDDRegPinned{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegPinned{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(44, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(45, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(46, "And \"Share BDDRegPinned{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegPinned{time}\" dialog")), "visible"));
    await session.step(47, "When user clicks on CANCEL button in \"Share BDDRegPinned{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegPinned{time}\" dialog"))));
    await session.step(48, "Then the \"Share BDDRegPinned{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegPinned{time}")));
    await session.step(49, "And the \"demog\" table of the \"BDDRegPinned{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog", session.text("BDDRegPinned{time}")));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
    await session.step(51, "When user closes all views", () => closeAllViews(page));
    await session.step(52, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(53, "And user enters \"BDDRegPinned{time}\" into gallery search", () => enterInto(page, session.text("BDDRegPinned{time}"), el("gallery search")));
    await session.step(54, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(55, "And user double-clicks on BDDRegPinned{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegPinned{time} gallery card"))));
    await session.step(56, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(57, "And the table should have 5850 rows", () => rowCount(page, 5850));
    await session.step(58, "And the table should have been reloaded by data sync", () => reloadedByDataSync(page));
    await session.step(62, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(63, "And the \"column order\" reading of grid should not include the text \"HEIGHT\"", () => readingExcludes(page, "column order", el("grid"), "HEIGHT"));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
