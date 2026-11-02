/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-known-failures.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space, core#976, core#1042]
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
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, openDataset, projectsOnServer, saveAsProject, spaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInSpace} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("The known failures of spaces", () => {
  const session = feature(test, "features/spaces/spaces-known-failures.feature", import.meta.url);
  test("Deleting a space keeps a project that was only linked into it", {tag: ["@serial", "@spaces", "@realizes:views.space", "@edge-cases", "@known-failure:50", "@realizes:core#976"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await knownFailure(async () => {
      await session.step(26, "Given no project named \"bddedgelinkproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgelinkproj{time}")));
      await session.step(27, "And a space named \"BDDEdgeLink{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeLink{time}")));
      await session.step(28, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(29, "When user saves the current view as project \"bddedgelinkproj{time}\"", () => saveAsProject(page, session.text("bddedgelinkproj{time}")));
      await session.step(30, "And user closes all views", () => closeAllViews(page));
      await session.step(31, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(32, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(33, "And user enters \"bddedgelinkproj{time}\" into gallery search", () => enterInto(page, session.text("bddedgelinkproj{time}"), el("gallery search")));
      await session.step(34, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(35, "Then bddedgelinkproj{time} gallery card should be visible", () => shouldBe(page, el(session.text("bddedgelinkproj{time} gallery card")), "visible"));
      await session.step(36, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(37, "When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree", () => dragTo(page, el(session.text("bddedgelinkproj{time} gallery card")), el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(38, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(39, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(40, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(41, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(42, "When user double-clicks on BDDEdgeLink{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(43, "Then the \"BDDEdgeLink{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeLink{time}")));
      await session.step(44, "And \"bddedgelinkproj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bddedgelinkproj{time}\" link in gallery")), "visible"));
      await session.step(45, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
      await session.step(46, "When user picks \"Delete Space\" from the context menu of BDDEdgeLink{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(47, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(48, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(49, "And 0 spaces named \"BDDEdgeLink{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeLink{time}")));
      await session.step(50, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
    }, 50);
  });
  test("A file linked from one space into another shows in the target space", {tag: ["@serial", "@spaces", "@realizes:views.space", "@edge-cases", "@known-failure:68", "@realizes:core#1042"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And the browse panel is open", () => browsePanelOpen(page));
    await knownFailure(async () => {
      await session.step(57, "Given a space named \"BDDEdgeFileSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeFileSrc{time}")));
      await session.step(58, "And a space named \"BDDEdgeFileDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeFileDst{time}")));
      await session.step(59, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeFileSrc{time}\" as \"edgelinkcars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeFileSrc{time}"), session.text("edgelinkcars{time}.csv")));
      await session.step(60, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(61, "When user double-clicks on BDDEdgeFileSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeFileSrc{time} tree node inside browse tree"))));
      await session.step(62, "Then the \"BDDEdgeFileSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeFileSrc{time}")));
      await session.step(63, "And \"edgelinkcars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgelinkcars{time}.csv\" link in gallery")), "visible"));
      await session.step(64, "When user drags \"edgelinkcars{time}.csv\" link in gallery to BDDEdgeFileDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"edgelinkcars{time}.csv\" link in gallery")), el(session.text("BDDEdgeFileDst{time} tree node inside browse tree"))));
      await session.step(65, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(66, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(67, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(68, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
    }, 68);
  });
});
