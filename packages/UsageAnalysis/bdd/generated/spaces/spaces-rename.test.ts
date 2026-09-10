/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-rename.feature
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
import {noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, close, doubleClickOn, enterInto, expand, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {errorBalloonText, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Renaming a space", () => {
  const session = feature(test, "features/spaces/spaces-rename.feature", import.meta.url);
  test("Renaming a space", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(17, "And no space named \"BDD-Ren, BDD-Ren-New, BDD-Other, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew\" is on the server", () => noSpaceOnServer(page, "BDD-Ren, BDD-Ren-New, BDD-Other, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew"));
    await run.scenario("The Rename dialog opens on the current name", async () => {
      await session.step(20, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(21, "And user enters \"BDD-Ren\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ren", el("Name input in Create Space dialog")));
      await session.step(22, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(23, "Then 1 space named \"BDD-Ren\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren"));
      await session.step(24, "When user picks \"Rename...\" from the context menu of BDD-Ren tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("BDD-Ren tree node inside browse tree")));
      await session.step(25, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(26, "And Name input in Rename project dialog should have value \"BDD-Ren\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), "BDD-Ren"));
    });
    await run.scenario("Cancelling the rename keeps the old name", async () => {
      await session.step(29, "When user enters \"BDD-Ren-New\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-New", el("Name input in Rename project dialog")));
      await session.step(30, "And user clicks on CANCEL button in Rename project dialog", () => clickOn(page, el("CANCEL button in Rename project dialog")));
      await session.step(31, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
      await session.step(32, "And 1 space named \"BDD-Ren\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren"));
      await session.step(33, "And 0 spaces named \"BDD-Ren-New\" should be on the server", () => spacesOnServer(page, 0, "BDD-Ren-New"));
      await session.step(34, "And BDD-Ren tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren tree node inside browse tree"), "visible"));
    });
    await run.scenario("A rename reaches the server and the tree", async () => {
      await session.step(37, "When user picks \"Rename...\" from the context menu of BDD-Ren tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("BDD-Ren tree node inside browse tree")));
      await session.step(38, "And user enters \"BDD-Ren-New\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-New", el("Name input in Rename project dialog")));
      await session.step(39, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(40, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
      await session.step(41, "And 1 space named \"BDD-Ren-New\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren-New"));
      await session.step(42, "And 0 spaces named \"BDD-Ren\" should be on the server", () => spacesOnServer(page, 0, "BDD-Ren"));
      await session.step(43, "And BDD-Ren-New tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren-New tree node inside browse tree"), "visible"));
      await session.step(44, "And BDD-Ren tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Ren tree node inside browse tree"), "absent"));
    });
    await run.scenario("Renaming onto an existing name is refused", async () => {
      await session.step(47, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(48, "And user enters \"BDD-Other\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Other", el("Name input in Create Space dialog")));
      await session.step(49, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(50, "Then 1 space named \"BDD-Other\" should be on the server", () => spacesOnServer(page, 1, "BDD-Other"));
      await session.step(51, "When user picks \"Rename...\" from the context menu of BDD-Ren-New tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("BDD-Ren-New tree node inside browse tree")));
      await session.step(52, "And user enters \"BDD-Other\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Other", el("Name input in Rename project dialog")));
      await session.step(53, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(54, "Then an error balloon containing \"already exists\" should have been shown", () => errorBalloonText(page, "already exists"));
      await session.step(55, "And 1 space named \"BDD-Ren-New\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren-New"));
      await session.step(56, "And 1 space named \"BDD-Other\" should be on the server", () => spacesOnServer(page, 1, "BDD-Other"));
      await session.step(57, "When user closes Rename project dialog", () => close(page, el("Rename project dialog")));
      await session.step(58, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
    });
    await run.scenario("A child space is renamed from its card in the parent", async () => {
      await session.step(61, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(62, "And user enters \"BDD-Ren-Parent\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ren-Parent", el("Name input in Create Space dialog")));
      await session.step(63, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(64, "Then 1 space named \"BDD-Ren-Parent\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren-Parent"));
      await session.step(65, "When user picks \"Create Child Space...\" from the context menu of BDD-Ren-Parent tree node inside browse tree", () => pickFromContextMenu(page, "Create Child Space...", el("BDD-Ren-Parent tree node inside browse tree")));
      await session.step(66, "And user enters \"BDD-Ren-Child\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ren-Child", el("Name input in Create Space dialog")));
      await session.step(67, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(68, "Then BDD-Ren-Child tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren-Child tree node inside browse tree"), "visible"));
      await session.step(69, "When user double-clicks on BDD-Ren-Parent tree node inside browse tree", () => doubleClickOn(page, el("BDD-Ren-Parent tree node inside browse tree")));
      await session.step(70, "Then BDD-Ren-Child link in space gallery should be visible", () => shouldBe(page, el("BDD-Ren-Child link in space gallery"), "visible"));
      await session.step(71, "When user picks \"Rename...\" from the context menu of BDD-Ren-Child link in space gallery", () => pickFromContextMenu(page, "Rename...", el("BDD-Ren-Child link in space gallery")));
      await session.step(72, "Then Name input in Rename project dialog should have value \"BDD-Ren-Child\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), "BDD-Ren-Child"));
      await session.step(73, "When user enters \"BDD-Ren-ChildNew\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-ChildNew", el("Name input in Rename project dialog")));
      await session.step(74, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(75, "Then BDD-Ren-ChildNew link in space gallery should be visible", () => shouldBe(page, el("BDD-Ren-ChildNew link in space gallery"), "visible"));
      await session.step(76, "And BDD-Ren-Child link in space gallery should be absent", () => shouldBe(page, el("BDD-Ren-Child link in space gallery"), "absent"));
      await session.step(77, "When user expands BDD-Ren-Parent tree node inside browse tree", () => expand(page, el("BDD-Ren-Parent tree node inside browse tree")));
      await session.step(78, "Then BDD-Ren-ChildNew tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren-ChildNew tree node inside browse tree"), "visible"));
      await session.step(79, "And BDD-Ren-Child tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Ren-Child tree node inside browse tree"), "absent"));
    });
    run.finish();
  });
});
