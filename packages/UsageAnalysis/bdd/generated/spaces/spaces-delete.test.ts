/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space, core#976]
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
import {clickOn, doubleClickOn, dragTo, enterInto, isExpanded, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, dialogCloses, noProjectOnServer, openDataset, projectsOnServer, saveAsProject, spaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a space", () => {
  const session = feature(test, "features/spaces/spaces-delete.feature", import.meta.url);
  test("Deleting a space", {tag: ["@journey", "@serial", "@spaces", "@realizes:views.space", "@edge-cases", "@known-failure:74", "@realizes:core#976"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(22, "Given user is logged in", () => loggedIn(page));
    await session.step(23, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(24, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(25, "And a space named \"BDD-Del\" is on the server", () => spaceOnServer(page, "BDD-Del"));
    await session.step(26, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("Deleting asks first", async () => {
      await session.step(29, "When user picks \"Delete Space\" from the context menu of Spaces---BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("Spaces---BDD-Del tree node inside browse tree")));
      await session.step(30, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(31, "And \"Are you sure?\" dialog should contain text \"BDD-Del\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "BDD-Del"));
    });
    await run.scenario("Cancelling the confirmation keeps the space", async () => {
      await session.step(34, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(35, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(36, "And 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(37, "And Spaces---BDD-Del tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Del tree node inside browse tree"), "visible"));
    });
    await run.scenario("Confirming removes it from the server and the tree", async () => {
      await session.step(40, "When user picks \"Delete Space\" from the context menu of Spaces---BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("Spaces---BDD-Del tree node inside browse tree")));
      await session.step(41, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(42, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(43, "And 0 spaces named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del"));
      await session.step(44, "And Spaces---BDD-Del tree node inside browse tree should be absent", () => shouldBe(page, el("Spaces---BDD-Del tree node inside browse tree"), "absent"));
    });
    await run.scenario("Deleting a space keeps a project that was only linked into it", async () => {
      await session.step(50, "Given no project named \"bddedgelinkproj{time}\" is on the server", () => noProjectOnServer(page, session.text("bddedgelinkproj{time}")));
      await session.step(51, "And a space named \"BDDEdgeLink{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeLink{time}")));
      await session.step(52, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(53, "When user saves the current view as project \"bddedgelinkproj{time}\"", () => saveAsProject(page, session.text("bddedgelinkproj{time}")));
      await session.step(54, "And user closes all views", () => closeAllViews(page));
      await session.step(55, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(56, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(57, "And user enters \"bddedgelinkproj{time}\" into gallery search", () => enterInto(page, session.text("bddedgelinkproj{time}"), el("gallery search")));
      await session.step(58, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(59, "Then bddedgelinkproj{time} gallery card should be visible", () => shouldBe(page, el(session.text("bddedgelinkproj{time} gallery card")), "visible"));
      await session.step(60, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(61, "When user drags bddedgelinkproj{time} gallery card to BDDEdgeLink{time} tree node inside browse tree", () => dragTo(page, el(session.text("bddedgelinkproj{time} gallery card")), el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(62, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(63, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(64, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(65, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(66, "When user double-clicks on BDDEdgeLink{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(67, "Then the \"BDDEdgeLink{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeLink{time}")));
      await session.step(68, "And \"bddedgelinkproj{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"bddedgelinkproj{time}\" link in gallery")), "visible"));
      await session.step(69, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
      await session.step(70, "When user picks \"Delete Space\" from the context menu of BDDEdgeLink{time} tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el(session.text("BDDEdgeLink{time} tree node inside browse tree"))));
      await session.step(71, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(72, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(73, "And 0 spaces named \"BDDEdgeLink{time}\" should be on the server", () => spacesOnServer(page, 0, session.text("BDDEdgeLink{time}")));
      await session.step(74, "And 1 project named \"bddedgelinkproj{time}\" should be on the server", () => projectsOnServer(page, 1, session.text("bddedgelinkproj{time}")));
    }, {knownFailure: true, line: 74});
    run.finish();
  });
});
