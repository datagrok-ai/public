/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-views.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, GROK-18454]
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
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, simpleModeOff, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: the view a reopened project makes current", () => {
  const session = feature(test, "features/projects/projects-regressions-views.feature", import.meta.url);
  test("A reopened project makes current the view that was current when it was saved", {tag: ["@serial", "@realizes:views.projects", "@realizes:GROK-18454"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "Given no project named \"BDDRegActive{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegActive{time}")));
    await session.step(27, "And simple mode is off", () => simpleModeOff(page));
    await session.step(28, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(29, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(30, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(31, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(32, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(33, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(34, "When user double-clicks Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(35, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(36, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(37, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(38, "When user double-clicks Files---Demo---iris.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---iris.csv tree node inside browse tree")));
    await session.step(39, "Then the \"iris\" view should be current", () => viewIsCurrent(page, "iris"));
    await session.step(40, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(41, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(42, "When user clicks on cars view", () => clickOn(page, el("cars view")));
    await session.step(43, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(44, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(45, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(46, "When user enters \"BDDRegActive{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegActive{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(47, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(48, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(49, "And \"Share BDDRegActive{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegActive{time}\" dialog")), "visible"));
    await session.step(50, "When user clicks on CANCEL button in \"Share BDDRegActive{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegActive{time}\" dialog"))));
    await session.step(51, "Then the \"Share BDDRegActive{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegActive{time}")));
    await session.step(52, "When user closes all views", () => closeAllViews(page));
    await session.step(53, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(54, "And user enters \"BDDRegActive{time}\" into gallery search", () => enterInto(page, session.text("BDDRegActive{time}"), el("gallery search")));
    await session.step(55, "And user double-clicks on BDDRegActive{time} gallery card", () => doubleClickOn(page, el(session.text("BDDRegActive{time} gallery card"))));
    await session.step(56, "Then demog view should be visible", () => shouldBe(page, el("demog view"), "visible"));
    await session.step(57, "And iris view should be visible", () => shouldBe(page, el("iris view"), "visible"));
    await session.step(58, "And the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(59, "And cars view should be selected", () => shouldBe(page, el("cars view"), "selected"));
  });
});
