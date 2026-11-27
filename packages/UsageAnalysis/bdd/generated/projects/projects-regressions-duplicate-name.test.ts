/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-regressions-duplicate-name.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects, core#546]
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
import {clickOn, doubleClickOn, enterInto, isExpanded, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openRememberedAddress, rememberAddress} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, urlShouldContain, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Projects regressions: two projects saved under one name", () => {
  const session = feature(test, "features/projects/projects-regressions-duplicate-name.feature", import.meta.url);
  test("The second project saved under a taken name opens from its own address", {tag: ["@realizes:views.projects", "@realizes:core#546"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(14, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(19, "Given no project named \"BDDRegDupName{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDRegDupName{time}")));
    await session.step(20, "When user double-clicks Files---Demo---demog-1000.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(21, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(22, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(23, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(24, "When user enters \"BDDRegDupName{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegDupName{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(25, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(26, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(27, "And \"Share BDDRegDupName{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegDupName{time}\" dialog")), "visible"));
    await session.step(28, "When user clicks on CANCEL button in \"Share BDDRegDupName{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegDupName{time}\" dialog"))));
    await session.step(29, "Then the \"Share BDDRegDupName{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegDupName{time}")));
    await session.step(30, "When user closes all views", () => closeAllViews(page));
    await session.step(31, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(32, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(33, "And user double-clicks Files---Demo---cars.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---cars.csv tree node inside browse tree")));
    await session.step(34, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
    await session.step(35, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
    await session.step(36, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
    await session.step(37, "When user enters \"BDDRegDupName{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDRegDupName{time}"), el("Name text input in \"Save project\" dialog")));
    await session.step(38, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
    await session.step(39, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    await session.step(40, "And \"Share BDDRegDupName{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDRegDupName{time}\" dialog")), "visible"));
    await session.step(41, "When user clicks on CANCEL button in \"Share BDDRegDupName{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDRegDupName{time}\" dialog"))));
    await session.step(42, "Then the \"Share BDDRegDupName{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDRegDupName{time}")));
    await session.step(43, "And the page address should contain \"/p/\"", () => urlShouldContain(page, "/p/"));
    await session.step(44, "When user remembers the page address", () => rememberAddress(page));
    await session.step(45, "And user closes all views", () => closeAllViews(page));
    await session.step(46, "And user opens the remembered address", () => openRememberedAddress(page));
    await session.step(47, "Then the \"cars\" view should be current", () => viewIsCurrent(page, "cars"));
  });
});
