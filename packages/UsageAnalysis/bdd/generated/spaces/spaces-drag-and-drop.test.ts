/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-drag-and-drop.feature
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
import {browsePanelOpen, noSpaceOnServer, nodeExpanded, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, dragTo, enterInto, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Putting files into a space by dragging them", () => {
  const session = feature(test, "features/spaces/spaces-drag-and-drop.feature", import.meta.url);
  test("Putting files into a space by dragging them", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And no space named \"BDD-DnD, BDD-DnD-Src\" is on the server", () => noSpaceOnServer(page, "BDD-DnD, BDD-DnD-Src"));
    await run.scenario("Two spaces and the demo files", async () => {
      await session.step(25, "When user picks \"Create Space...\" from the context menu of Spaces tree node", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node")));
      await session.step(26, "And user enters \"BDD-DnD\" into Name input in Create Space dialog", () => enterInto(page, "BDD-DnD", el("Name input in Create Space dialog")));
      await session.step(27, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(28, "Then 1 space named \"BDD-DnD\" should be on the server", () => spacesOnServer(page, 1, "BDD-DnD"));
      await session.step(29, "When user picks \"Create Space...\" from the context menu of Spaces tree node", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node")));
      await session.step(30, "And user enters \"BDD-DnD-Src\" into Name input in Create Space dialog", () => enterInto(page, "BDD-DnD-Src", el("Name input in Create Space dialog")));
      await session.step(31, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(32, "Then 1 space named \"BDD-DnD-Src\" should be on the server", () => spacesOnServer(page, 1, "BDD-DnD-Src"));
      await session.step(33, "Given the \"Files\" tree node is expanded", () => nodeExpanded(page, "Files"));
      await session.step(34, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(35, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(36, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The dialog offers Link, Copy and Move, and names the target", async () => {
      await session.step(39, "When user drags demog.csv link in gallery to BDD-DnD tree node", () => dragTo(page, el("demog.csv link in gallery"), el("BDD-DnD tree node")));
      await session.step(40, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(41, "And Move entity dialog should contain text \"BDDDnD\"", () => shouldContainText(page, el("Move entity dialog"), "BDDDnD"));
      await session.step(42, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
    });
    await run.scenario("Cancelling leaves the space empty", async () => {
      await session.step(45, "When user clicks on CANCEL button in Move entity dialog", () => clickOn(page, el("CANCEL button in Move entity dialog")));
      await session.step(46, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(47, "When user double-clicks on BDD-DnD tree node", () => doubleClickOn(page, el("BDD-DnD tree node")));
      await session.step(48, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(49, "And demog.csv link in gallery should be absent", () => shouldBe(page, el("demog.csv link in gallery"), "absent"));
    });
    await run.scenario("A copied file lands in the space", async () => {
      await session.step(52, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(53, "And user drags demog.csv link in gallery to BDD-DnD tree node", () => dragTo(page, el("demog.csv link in gallery"), el("BDD-DnD tree node")));
      await session.step(54, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(55, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(56, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(57, "When user double-clicks on BDD-DnD tree node", () => doubleClickOn(page, el("BDD-DnD tree node")));
      await session.step(58, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(59, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file's details show in the context panel", async () => {
      await session.step(62, "When user clicks on demog.csv link in gallery", () => clickOn(page, el("demog.csv link in gallery")));
      await session.step(63, "Then context panel should be visible", () => shouldBe(page, el("context panel"), "visible"));
      await session.step(64, "And context panel should contain text \"demog\"", () => shouldContainText(page, el("context panel"), "demog"));
    });
    await run.scenario("A second file joins the first", async () => {
      await session.step(67, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(68, "And user drags TSLA.csv link in gallery to BDD-DnD tree node", () => dragTo(page, el("TSLA.csv link in gallery"), el("BDD-DnD tree node")));
      await session.step(69, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(70, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(71, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(72, "When user double-clicks on BDD-DnD tree node", () => doubleClickOn(page, el("BDD-DnD tree node")));
      await session.step(73, "Then demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(74, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy is made in the source space to move later", async () => {
      await session.step(77, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(78, "And user drags beer.csv link in gallery to BDD-DnD-Src tree node", () => dragTo(page, el("beer.csv link in gallery"), el("BDD-DnD-Src tree node")));
      await session.step(79, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(80, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(81, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(82, "When user double-clicks on BDD-DnD-Src tree node", () => doubleClickOn(page, el("BDD-DnD-Src tree node")));
      await session.step(83, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(84, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("Moving takes the file out of the source space", async () => {
      await session.step(87, "When user drags beer.csv link in gallery to BDD-DnD tree node", () => dragTo(page, el("beer.csv link in gallery"), el("BDD-DnD tree node")));
      await session.step(88, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(89, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(90, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(91, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(92, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(93, "When user double-clicks on BDD-DnD tree node", () => doubleClickOn(page, el("BDD-DnD tree node")));
      await session.step(94, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(95, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
      await session.step(96, "When user double-clicks on BDD-DnD-Src tree node", () => doubleClickOn(page, el("BDD-DnD-Src tree node")));
      await session.step(97, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(98, "And beer.csv link in gallery should be absent", () => shouldBe(page, el("beer.csv link in gallery"), "absent"));
    });
    await run.scenario("The demo files kept their originals", async () => {
      await session.step(101, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(102, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(103, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(104, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(105, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    run.finish();
  });
});
