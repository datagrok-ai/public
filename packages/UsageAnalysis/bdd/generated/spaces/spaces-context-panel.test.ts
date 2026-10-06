/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-context-panel.feature
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
import {clickOn, doubleClickOn, followingShouldBe, isExpanded, shouldBe, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, contextPanelOpen, contextPanelShows, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("What the context panel says about a space", () => {
  const session = feature(test, "features/spaces/spaces-context-panel.feature", import.meta.url);
  test("What the context panel says about a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(18, "And a space named \"BDD-CP-Root\" is on the server", () => spaceOnServer(page, "BDD-CP-Root"));
    await session.step(19, "And a space named \"BDD-CP-One\" under \"BDD-CP-Root\" is on the server", () => childSpaceOnServer(page, "BDD-CP-One", "BDD-CP-Root"));
    await session.step(20, "And a space named \"BDD-CP-Two\" under \"BDD-CP-Root\" is on the server", () => childSpaceOnServer(page, "BDD-CP-Two", "BDD-CP-Root"));
    await session.step(21, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The root's view lists its two children", async () => {
      await session.step(24, "When user double-clicks on Spaces---BDD-CP-Root tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-CP-Root tree node inside browse tree")));
      await session.step(25, "Then the \"BDD-CP-Root\" view should be current", () => viewIsCurrent(page, "BDD-CP-Root"));
      await session.step(26, "And BDD-CP-One link in gallery should be visible", () => shouldBe(page, el("BDD-CP-One link in gallery"), "visible"));
      await session.step(27, "And BDD-CP-Two link in gallery should be visible", () => shouldBe(page, el("BDD-CP-Two link in gallery"), "visible"));
    });
    await run.scenario("Selecting a space shows its details", async () => {
      await session.step(30, "When user clicks on Spaces---BDD-CP-Root tree node inside browse tree", () => clickOn(page, el("Spaces---BDD-CP-Root tree node inside browse tree")));
      await session.step(31, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(32, "And the context panel should show \"BDD-CP-Root\"", () => contextPanelShows(page, "BDD-CP-Root"));
      await session.step(33, "And \"Details\" accordion header in context panel should be visible", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "visible"));
    });
    await run.scenario("The panel carries the sections a space has", async () => {
      await session.step(40, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Details\" accordion header in context panel"],["\"Content\" accordion header in context panel"],["\"Sharing\" accordion header in context panel"],["\"Chats\" accordion header in context panel"]]), [["\"Details\" accordion header in context panel"],["\"Content\" accordion header in context panel"],["\"Sharing\" accordion header in context panel"],["\"Chats\" accordion header in context panel"]]);
      await session.step(45, "And \"Activity\" accordion header in context panel should be present", () => shouldBe(page, el("\"Activity\" accordion header in context panel"), "present"));
    });
    await run.scenario("Clicking one child, then the other, switches the panel", async () => {
      await session.step(48, "When user double-clicks on Spaces---BDD-CP-Root tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-CP-Root tree node inside browse tree")));
      await session.step(49, "Then the \"BDD-CP-Root\" view should be current", () => viewIsCurrent(page, "BDD-CP-Root"));
      await session.step(50, "When user clicks on BDD-CP-One link in gallery", () => clickOn(page, el("BDD-CP-One link in gallery")));
      await session.step(51, "Then the context panel should show \"BDD-CP-One\"", () => contextPanelShows(page, "BDD-CP-One"));
      await session.step(52, "When user clicks on BDD-CP-Two link in gallery", () => clickOn(page, el("BDD-CP-Two link in gallery")));
      await session.step(53, "Then the context panel should show \"BDD-CP-Two\"", () => contextPanelShows(page, "BDD-CP-Two"));
      await session.step(54, "And context panel should not contain text \"BDD-CP-One\"", () => shouldNotContainText(page, el("context panel"), "BDD-CP-One"));
      await session.step(55, "When user clicks on BDD-CP-One link in gallery", () => clickOn(page, el("BDD-CP-One link in gallery")));
      await session.step(56, "Then the context panel should show \"BDD-CP-One\"", () => contextPanelShows(page, "BDD-CP-One"));
      await session.step(57, "And context panel should not contain text \"BDD-CP-Two\"", () => shouldNotContainText(page, el("context panel"), "BDD-CP-Two"));
    });
    run.finish();
  });
});
