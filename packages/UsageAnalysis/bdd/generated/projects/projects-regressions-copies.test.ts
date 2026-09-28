/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-copies.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-19792]
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
import {galleryPictures, noErrorsButPreviewNoise} from '../../bindings/projects-regressions.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, openDataset, projectsOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: two copies of a project saved without renaming", () => {
  const session = feature(test, "features/projects/projects-regressions-copies.feature", import.meta.url);
  test("Projects regressions: two copies of a project saved without renaming", {tag: ["@journey", "@serial", "@realizes:views.projects", "@realizes:GROK-19792", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await run.scenario("Two copies saved without renaming each keep a preview of their own", async () => {
      await session.step(32, "Given no project named \"BDDRegCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegCopy{time}")));
      await session.step(33, "And no project named \"Copy of BDDRegCopy{time}\" is on the server", () => noProjectOnServer(page, session.text("Copy of BDDRegCopy{time}")));
      await session.step(34, "And no project named \"Copy of BDDRegCopy{time} (2)\" is on the server", () => noProjectOnServer(page, session.text("Copy of BDDRegCopy{time} (2)")));
      await session.step(35, "And user opens demog dataset", () => openDataset(page, ds("demog")));
      await session.step(36, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(37, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(38, "When user enters \"BDDRegCopy{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegCopy{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(39, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(40, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(41, "And \"Share BDDRegCopy{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegCopy{time}\" dialog")), "visible"));
      await session.step(42, "When user clicks on CANCEL button in \"Share BDDRegCopy{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegCopy{time}\" dialog"))));
      await session.step(43, "Then the \"Share BDDRegCopy{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegCopy{time}")));
      await session.step(44, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(45, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(46, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(47, "Then Name text input in \"Save project\" dialog should have value \"Copy of BDDRegCopy{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("Copy of BDDRegCopy{time}")));
      await session.step(48, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(49, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(50, "And an info balloon containing 'Project \"Copy of BDDRegCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"Copy of BDDRegCopy{time}\" uploaded")));
      await session.step(51, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(52, "When user closes all views", () => closeAllViews(page));
      await session.step(53, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(54, "And user enters \"BDDRegCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDRegCopy{time}"), el("gallery search")));
      await session.step(55, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(56, "And user double-clicks on BDDRegCopy{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegCopy{time} gallery card"))));
      await session.step(57, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(58, "And the table should have 5850 rows", () => rowCount(page, 5850));
      await session.step(59, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(60, "When user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")));
      await session.step(61, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(62, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(63, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(64, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(65, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(66, "And an info balloon containing 'Project \"BDDRegCopy{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDRegCopy{time}\" uploaded")));
      await session.step(67, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(68, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(69, "When user selects \"Save a copy\" in radio input in \"Save project\" dialog", () => selectIn(page, "Save a copy", el("radio input in \"Save project\" dialog")));
      await session.step(70, "Then Name text input in \"Save project\" dialog should have value \"Copy of BDDRegCopy{time}\"", () => shouldHaveValue(page, el("Name text input in \"Save project\" dialog"), session.text("Copy of BDDRegCopy{time}")));
      await session.step(71, "When user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(72, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(73, "And an info balloon containing 'Project \"Copy of BDDRegCopy' should have been shown", () => infoBalloonText(page, "Project \"Copy of BDDRegCopy"));
      await session.step(74, "And no errors but the project preview's should have been logged", () => noErrorsButPreviewNoise(page));
      await session.step(75, "When user closes all views", () => closeAllViews(page));
      await session.step(76, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(77, "And user enters \"BDDRegCopy{time}\" into gallery search", () => enterInto(page, session.text("BDDRegCopy{time}"), el("gallery search")));
      await session.step(78, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(79, "Then the 3 gallery cards should show 3 different pictures", () => galleryPictures(page, 3, 3));
    });
    await run.scenario("The second copy saved without renaming gets a name of its own", async () => {
      await session.step(84, "Then 1 project named \"Copy of BDDRegCopy{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("Copy of BDDRegCopy{time}")));
      await session.step(85, "And 1 project named \"Copy of BDDRegCopy{time} (2)\" should be on the server", () => projectsOnServer(page, 1, session.text("Copy of BDDRegCopy{time} (2)")));
    }, {knownFailure: true});
    run.finish();
  });
});
