/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-hierarchy.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, enterInto, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, noSpaceOnServer, spacesOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Nested spaces and moving between them", () => {
  const session = feature(test, "features/spaces/spaces-hierarchy.feature", import.meta.url);
  test("Nested spaces and moving between them", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(16, "And no space named \"BDD-Hier-Root, BDD-Hier-Child, BDD-Hier-Grand\" is on the server", () => noSpaceOnServer(page, "BDD-Hier-Root, BDD-Hier-Child, BDD-Hier-Grand"));
    await run.scenario("A child is created from the tree", async () => {
      await session.step(19, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(20, "And user enters \"BDD-Hier-Root\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Hier-Root", el("Name input in Create Space dialog")));
      await session.step(21, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(22, "Then 1 space named \"BDD-Hier-Root\" should be on the server", () => spacesOnServer(page, 1, "BDD-Hier-Root"));
      await session.step(23, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(24, "When user picks \"Create Child Space...\" from the context menu of BDD-Hier-Root tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Hier-Root tree node inside browse tree")));
      await session.step(25, "And user enters \"BDD-Hier-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Hier-Child", el("Name input in Create Space dialog")));
      await session.step(26, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(27, "Then BDD-Hier-Child tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Hier-Child tree node inside browse tree"), "visible"));
    });
    await run.scenario("The parent's view lists the child", async () => {
      await session.step(30, "When user double-clicks on BDD-Hier-Root tree node inside browse tree", () => doubleClickOn(page, el("BDD-Hier-Root tree node inside browse tree")));
      await session.step(31, "Then the \"BDD-Hier-Root\" view should be current", () => viewIsCurrent(page, "BDD-Hier-Root"));
      await session.step(32, "And space gallery should be visible", () => shouldBe(page, el("space gallery"), "visible"));
      await session.step(33, "And BDD-Hier-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Hier-Child link in space gallery"), "visible"));
    });
    await run.scenario("A grandchild is created from the child's card", async () => {
      await session.step(36, "When user picks \"Create Child Space...\" from the context menu of BDD-Hier-Child link in space gallery", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Hier-Child link in space gallery")));
      await session.step(37, "And user enters \"BDD-Hier-Grand\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Hier-Grand", el("Name input in Create Space dialog")));
      await session.step(38, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(39, "Then the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(40, "And BDD-Hier-Grand tree node inside browse tree should be present", () => shouldBe(page, el("BDD-Hier-Grand tree node inside browse tree"), "present"));
    });
    await run.scenario("Opening the child shows the grandchild", async () => {
      await session.step(43, "When user double-clicks on BDD-Hier-Child link in space gallery", () => doubleClickOn(page, el("BDD-Hier-Child link in space gallery")));
      await session.step(44, "Then the \"BDD-Hier-Child\" view should be current", () => viewIsCurrent(page, "BDD-Hier-Child"));
      await session.step(45, "And BDD-Hier-Grand link in space gallery should be visible", () => shouldBe(page, el("BDD-Hier-Grand link in space gallery"), "visible"));
      await session.step(46, "And BDD-Hier-Child link in space gallery should be absent", () => shouldBe(page, el("BDD-Hier-Child link in space gallery"), "absent"));
    });
    await run.scenario("Opening the grandchild leaves an empty space", async () => {
      await session.step(49, "When user double-clicks on BDD-Hier-Grand link in space gallery", () => doubleClickOn(page, el("BDD-Hier-Grand link in space gallery")));
      await session.step(50, "Then the \"BDD-Hier-Grand\" view should be current", () => viewIsCurrent(page, "BDD-Hier-Grand"));
      await session.step(51, "And BDD-Hier-Grand link in space gallery should be absent", () => shouldBe(page, el("BDD-Hier-Grand link in space gallery"), "absent"));
    });
    await run.scenario("Going back up the tree finds the content again", async () => {
      await session.step(54, "When user double-clicks on BDD-Hier-Root tree node inside browse tree", () => doubleClickOn(page, el("BDD-Hier-Root tree node inside browse tree")));
      await session.step(55, "Then the \"BDD-Hier-Root\" view should be current", () => viewIsCurrent(page, "BDD-Hier-Root"));
      await session.step(56, "And BDD-Hier-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Hier-Child link in space gallery"), "visible"));
      await session.step(57, "When user double-clicks on BDD-Hier-Child link in space gallery", () => doubleClickOn(page, el("BDD-Hier-Child link in space gallery")));
      await session.step(58, "Then the \"BDD-Hier-Child\" view should be current", () => viewIsCurrent(page, "BDD-Hier-Child"));
      await session.step(59, "And BDD-Hier-Grand link in space gallery should be visible", () => shouldBe(page, el("BDD-Hier-Grand link in space gallery"), "visible"));
    });
    await run.scenario("Search still works after the walk", async () => {
      await session.step(62, "When user enters \"zzz-no-such-space\" into space search", () => enterInto(page, "zzz-no-such-space", el("space search")));
      await session.step(63, "Then BDD-Hier-Grand link in space gallery should be absent", () => shouldBe(page, el("BDD-Hier-Grand link in space gallery"), "absent"));
      await session.step(64, "When user clears space search", () => clearField(page, el("space search")));
      await session.step(65, "Then BDD-Hier-Grand link in space gallery should be visible", () => shouldBe(page, el("BDD-Hier-Grand link in space gallery"), "visible"));
    });
    run.finish();
  });
});
