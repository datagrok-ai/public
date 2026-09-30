/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-derived.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-19580]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {clearSavedParameters} from '../../bindings/pivot-table.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, downloadContains, enterInto, fileDownloaded, isExpanded, shouldBe, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, savedWithDataSync, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noBalloons, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: a data-sync project with a published pivot saved as a zip", () => {
  const session = feature(test, "features/projects/projects-regressions-derived.feature", import.meta.url);
  test("A data-sync project with a published pivot is saved as a zip from the gallery", {tag: ["@serial", "@realizes:views.projects", "@realizes:GROK-19580"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(27, "Given no project named \"BDDRegZip{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegZip{time}")));
    await session.step(28, "And user clears the saved pivot table parameters", () => clearSavedParameters(page));
    await session.step(29, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(31, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(32, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(33, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(34, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(35, "When user clicks on \"pivot table\" icon in toolbox", () => clickOn(page, el("\"pivot table\" icon in toolbox")));
    await session.step(36, "Then pivot table viewer should be visible", () => shouldBe(page, el("pivot table viewer"), "visible"));
    await session.step(37, "When user clicks on the \"add to workspace\" area of pivot table viewer", () => clickArea(page, "add to workspace", el("pivot table viewer")));
    await session.step(38, "Then the \"demog aggregation\" view should be current", () => viewIsCurrent(page, "demog aggregation"));
    await session.step(39, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(40, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(41, "And \"Creation script\" button in \"demog\" project table in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Creation script\" button in \"demog\" project table in \"Save project\" dialog"), "visible"));
    await session.step(42, "And Data sync switch in \"demog\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog\" project table in \"Save project\" dialog"), "checked"));
    await session.step(43, "And Data sync switch in \"demog aggregation\" project table in \"Save project\" dialog should be checked", () => shouldBe(page, el("Data sync switch in \"demog aggregation\" project table in \"Save project\" dialog"), "checked"));
    await session.step(44, "When user enters \"BDDRegZip{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegZip{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(45, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(46, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(47, "And \"Share BDDRegZip{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegZip{time}\" dialog")), "visible"));
    await session.step(48, "When user clicks on CANCEL button in \"Share BDDRegZip{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegZip{time}\" dialog"))));
    await session.step(49, "Then the \"Share BDDRegZip{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegZip{time}")));
    await session.step(50, "And the \"demog aggregation\" table of the \"BDDRegZip{time}\" project should be saved with data sync", () => savedWithDataSync(page, "demog aggregation", session.text("BDDRegZip{time}")));
    await session.step(51, "When user closes all views", () => closeAllViews(page));
    await session.step(52, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(53, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(54, "And user enters \"BDDRegZip{time}\" into gallery search", () => enterInto(page, session.text("BDDRegZip{time}"), el("gallery search")));
    await session.step(55, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(56, "Given user watches downloads", () => watchDownloads(page));
    await session.step(57, "When user picks \"Save as Zip\" from the context menu of BDDRegZip{time} gallery card", () => pickFromContextMenu(page, "Save as Zip", el(session.text("BDDRegZip{time} gallery card"))));
    await session.step(58, "Then a file \"BDDRegZip{time}.zip\" should have been downloaded", () => fileDownloaded(page, session.text("BDDRegZip{time}.zip")));
    await session.step(59, "And the downloaded file \"BDDRegZip{time}.zip\" should contain the text \"DemogAggregation = Aggregate(\"", () => downloadContains(page, session.text("BDDRegZip{time}.zip"), "DemogAggregation = Aggregate("));
    await session.step(60, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
