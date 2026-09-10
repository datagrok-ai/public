/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-context-panel.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {createDialogCloses, noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, enterInto, followingShouldBe, shouldBe, shouldContainText, shouldNotContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("What the context panel says about a space", () => {
  const session = feature(test, "features/spaces/spaces-context-panel.feature", import.meta.url);
  test("What the context panel says about a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And no space named \"BDD-CP-Root, BDD-CP-One, BDD-CP-Two\" is on the server", () => noSpaceOnServer(page, "BDD-CP-Root, BDD-CP-One, BDD-CP-Two"));
    await run.scenario("Two children to switch between", async () => {
      await session.step(19, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(20, "And user enters \"BDD-CP-Root\" into Name input in Create Space dialog", () => enterInto(page, "BDD-CP-Root", el("Name input in Create Space dialog")));
      await session.step(21, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(22, "Then 1 space named \"BDD-CP-Root\" should be on the server", () => spacesOnServer(page, 1, "BDD-CP-Root"));
      await session.step(23, "When user picks \"Create Child Space...\" from the context menu of BDD-CP-Root tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-CP-Root tree node inside browse tree")));
      await session.step(24, "And user enters \"BDD-CP-One\" into Name input in Create Space dialog", () => enterInto(page, "BDD-CP-One", el("Name input in Create Space dialog")));
      await session.step(25, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(26, "Then the Create Space dialog should close", () => createDialogCloses(page));
      await session.step(27, "When user picks \"Create Child Space...\" from the context menu of BDD-CP-Root tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-CP-Root tree node inside browse tree")));
      await session.step(28, "And user enters \"BDD-CP-Two\" into Name input in Create Space dialog", () => enterInto(page, "BDD-CP-Two", el("Name input in Create Space dialog")));
      await session.step(29, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(30, "Then the Create Space dialog should close", () => createDialogCloses(page));
      await session.step(31, "When user double-clicks on BDD-CP-Root tree node inside browse tree", () => doubleClickOn(page, el("BDD-CP-Root tree node inside browse tree")));
      await session.step(32, "Then the \"BDD-CP-Root\" view should be current", () => viewIsCurrent(page, "BDD-CP-Root"));
      await session.step(33, "And BDD-CP-One link in space gallery should be visible", () => shouldBe(page, el("BDD-CP-One link in space gallery"), "visible"));
      await session.step(34, "And BDD-CP-Two link in space gallery should be visible", () => shouldBe(page, el("BDD-CP-Two link in space gallery"), "visible"));
    });
    await run.scenario("Selecting a space shows its details", async () => {
      await session.step(37, "When user clicks on BDD-CP-Root tree node inside browse tree", () => clickOn(page, el("BDD-CP-Root tree node inside browse tree")));
      await session.step(38, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(39, "And context panel should contain text \"BDD-CP-Root\"", () => shouldContainText(page, el("context panel"), "BDD-CP-Root"));
      await session.step(40, "And \"Details\" accordion header in context panel should be visible", () => shouldBe(page, el("\"Details\" accordion header in context panel"), "visible"));
    });
    await run.scenario("The panel carries the sections a space has", async () => {
      await session.step(43, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["\"Details\" accordion header in context panel"],["\"Content\" accordion header in context panel"],["\"Activity\" accordion header in context panel"],["\"Sharing\" accordion header in context panel"],["\"Chats\" accordion header in context panel"]]));
    });
    await run.scenario("Clicking one child, then the other, switches the panel", async () => {
      await session.step(51, "When user double-clicks on BDD-CP-Root tree node inside browse tree", () => doubleClickOn(page, el("BDD-CP-Root tree node inside browse tree")));
      await session.step(52, "Then the \"BDD-CP-Root\" view should be current", () => viewIsCurrent(page, "BDD-CP-Root"));
      await session.step(53, "When user clicks on BDD-CP-One link in space gallery", () => clickOn(page, el("BDD-CP-One link in space gallery")));
      await session.step(54, "Then context panel should contain text \"BDD-CP-One\"", () => shouldContainText(page, el("context panel"), "BDD-CP-One"));
      await session.step(55, "When user clicks on BDD-CP-Two link in space gallery", () => clickOn(page, el("BDD-CP-Two link in space gallery")));
      await session.step(56, "Then context panel should contain text \"BDD-CP-Two\"", () => shouldContainText(page, el("context panel"), "BDD-CP-Two"));
      await session.step(57, "And context panel should not contain text \"BDD-CP-One\"", () => shouldNotContainText(page, el("context panel"), "BDD-CP-One"));
      await session.step(58, "When user clicks on BDD-CP-One link in space gallery", () => clickOn(page, el("BDD-CP-One link in space gallery")));
      await session.step(59, "Then context panel should contain text \"BDD-CP-One\"", () => shouldContainText(page, el("context panel"), "BDD-CP-One"));
      await session.step(60, "And context panel should not contain text \"BDD-CP-Two\"", () => shouldNotContainText(page, el("context panel"), "BDD-CP-Two"));
    });
    run.finish();
  });
});
