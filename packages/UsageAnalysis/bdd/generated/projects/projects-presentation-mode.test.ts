/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/projects/projects-presentation-mode.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.projects]
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
import {clickOn, doubleClickOn, enterInto, hoverOver, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, noProjectOnServer, projectsOnServer, toolboxPaneShown, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {infoBalloonText, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("The Presentation mode switch of a saved project's Save dialog", () => {
  const session = feature(test, "features/projects/projects-presentation-mode.feature", import.meta.url);
  test("The Presentation mode switch of a saved project's Save dialog", {tag: ["@journey", "@serial", "@realizes:views.projects"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And no project named \"BDDPresentProj{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDPresentProj{time}")));
    await run.scenario("demog with a scatter plot is saved as a project", async () => {
      await session.step(27, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(28, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(29, "When user double-clicks Files---Demo---demog.csv tree node inside browse tree", () => doubleClickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
      await session.step(30, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(31, "Given the toolbox pane is shown", () => toolboxPaneShown(page));
      await session.step(32, "When user clicks on \"scatter plot\" icon in toolbox", () => clickOn(page, el("\"scatter plot\" icon in toolbox")));
      await session.step(33, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(34, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(35, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(36, "When user enters \"BDDPresentProj{time}\" into Name text input in \"Save project\" dialog", () => enterInto(page, session.text("BDDPresentProj{time}"), el("Name text input in \"Save project\" dialog")));
      await session.step(37, "And user clicks on OK button in \"Save project\" dialog", () => clickOn(page, el("OK button in \"Save project\" dialog")));
      await session.step(38, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
      await session.step(39, "And an info balloon containing 'Project \"BDDPresentProj{time}\" uploaded' should have been shown", () => infoBalloonText(page, session.text("Project \"BDDPresentProj{time}\" uploaded")));
      await session.step(40, "And \"Share BDDPresentProj{time}\" dialog should be visible", () => shouldBe(page, el(session.text("\"Share BDDPresentProj{time}\" dialog")), "visible"));
      await session.step(41, "When user clicks on CANCEL button in \"Share BDDPresentProj{time}\" dialog", () => clickOn(page, el(session.text("CANCEL button in \"Share BDDPresentProj{time}\" dialog"))));
      await session.step(42, "Then the \"Share BDDPresentProj{time}\" dialog should close", () => dialogCloses(page, session.text("Share BDDPresentProj{time}")));
    });
    await run.scenario("The project opens in design mode", async () => {
      await session.step(45, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(46, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(47, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(48, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(49, "And user enters \"BDDPresent\" into gallery search", () => enterInto(page, "BDDPresent", el("gallery search")));
      await session.step(50, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(51, "And user double-clicks on BDDPresentProj{time} gallery card", () => doubleClickOn(page, el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(52, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
      await session.step(53, "And scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
      await session.step(54, "And \"Data\" menu item should be visible", () => shouldBe(page, el("\"Data\" menu item"), "visible"));
      await session.step(55, "And browse tab should be visible", () => shouldBe(page, el("browse tab"), "visible"));
      await session.step(56, "And status bar should be visible", () => shouldBe(page, el("status bar"), "visible"));
      await session.step(57, "And toolbox tab should be visible", () => shouldBe(page, el("toolbox tab"), "visible"));
    });
    await run.scenario("The Save dialog offers Presentation mode, and its tooltip tells how to get back", async () => {
      await session.step(60, "When user clicks on Save ribbon item", () => clickOn(page, el("Save ribbon item")));
      await session.step(61, "Then \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Save project\" dialog"), "visible"));
      await session.step(62, "And \"Presentation mode\" input in \"Save project\" dialog should be visible", () => shouldBe(page, el("\"Presentation mode\" input in \"Save project\" dialog"), "visible"));
      await session.step(63, "When user hovers over \"Presentation mode\" input in \"Save project\" dialog", () => hoverOver(page, el("\"Presentation mode\" input in \"Save project\" dialog")));
      await session.step(64, "Then tooltip should contain text \"visualization\"", () => shouldContainText(page, el("tooltip"), "visualization"));
      await session.step(65, "And tooltip should contain text \"Design mode\"", () => shouldContainText(page, el("tooltip"), "Design mode"));
      await session.step(66, "When user clicks on CANCEL button in \"Save project\" dialog", () => clickOn(page, el("CANCEL button in \"Save project\" dialog")));
      await session.step(67, "Then the \"Save project\" dialog should close", () => dialogCloses(page, "Save project"));
    });
    await run.scenario("The project is deleted from its card", async () => {
      await session.step(70, "When user picks \"Close All\" from the context menu of browse tab", () => pickFromContextMenu(page, "Close All", el("browse tab")));
      await session.step(71, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
      await session.step(72, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(73, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(74, "And user enters \"BDDPresentProj{time}\" into gallery search", () => enterInto(page, session.text("BDDPresentProj{time}"), el("gallery search")));
      await session.step(75, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(76, "And user picks \"Delete Project\" from the context menu of BDDPresentProj{time} gallery card", () => pickFromContextMenu(page, "Delete Project", el(session.text("BDDPresentProj{time} gallery card"))));
      await session.step(77, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(78, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(79, "And 0 projects named \"BDDPresentProj{time}\" should be on the server", () => projectsOnServer(page, 0, session.text("BDDPresentProj{time}")));
    });
    run.finish();
  });
});
