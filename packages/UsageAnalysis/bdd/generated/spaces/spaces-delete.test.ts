/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-delete.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/biostructure.js';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, isExpanded, shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, spaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Deleting a space", () => {
  const session = feature(test, "features/spaces/spaces-delete.feature", import.meta.url);
  test("Deleting a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(13, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(14, "And a space named \"BDD-Del\" is on the server", () => spaceOnServer(page, "BDD-Del"));
    await session.step(15, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("Deleting asks first", async () => {
      await session.step(18, "When user picks \"Delete Space\" from the context menu of Spaces---BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("Spaces---BDD-Del tree node inside browse tree")));
      await session.step(19, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(20, "And \"Are you sure?\" dialog should contain text \"BDD-Del\"", () => shouldContainText(page, el("\"Are you sure?\" dialog"), "BDD-Del"));
    });
    await run.scenario("Cancelling the confirmation keeps the space", async () => {
      await session.step(23, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(24, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(25, "And 1 space named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 1, "BDD-Del"));
      await session.step(26, "And Spaces---BDD-Del tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Del tree node inside browse tree"), "visible"));
    });
    await run.scenario("Confirming removes it from the server and the tree", async () => {
      await session.step(29, "When user picks \"Delete Space\" from the context menu of Spaces---BDD-Del tree node inside browse tree", () => pickFromContextMenu(page, "Delete Space", el("Spaces---BDD-Del tree node inside browse tree")));
      await session.step(30, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(31, "Then the \"Are you sure?\" dialog should close", () => dialogCloses(page, "Are you sure?"));
      await session.step(32, "And 0 spaces named \"BDD-Del\" should be on the server", () => spacesOnServer(page, 0, "BDD-Del"));
      await session.step(33, "And Spaces---BDD-Del tree node inside browse tree should be absent", () => shouldBe(page, el("Spaces---BDD-Del tree node inside browse tree"), "absent"));
    });
    run.finish();
  });
});
