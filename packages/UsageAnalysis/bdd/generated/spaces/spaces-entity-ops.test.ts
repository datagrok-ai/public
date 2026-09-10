/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-entity-ops.feature
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
import {clickOn, doubleClickOn, dragTo, enterInto, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Working with what a space holds", () => {
  const session = feature(test, "features/spaces/spaces-entity-ops.feature", import.meta.url);
  test("Working with what a space holds", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 7, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(14, "And no space named \"BDD-Ops, BDD-Ops-Copy\" is on the server", () => noSpaceOnServer(page, "BDD-Ops, BDD-Ops-Copy"));
    await run.scenario("A space with two files", async () => {
      await session.step(17, "When user picks \"Create Space...\" from the context menu of Spaces tree node", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node")));
      await session.step(18, "And user enters \"BDD-Ops\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ops", el("Name input in Create Space dialog")));
      await session.step(19, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(20, "Then 1 space named \"BDD-Ops\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ops"));
      await session.step(21, "When user picks \"Create Space...\" from the context menu of Spaces tree node", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node")));
      await session.step(22, "And user enters \"BDD-Ops-Copy\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ops-Copy", el("Name input in Create Space dialog")));
      await session.step(23, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(24, "Then 1 space named \"BDD-Ops-Copy\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ops-Copy"));
      await session.step(25, "Given the \"Files\" tree node is expanded", () => nodeExpanded(page, "Files"));
      await session.step(26, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(27, "And user drags TSLA.csv link in gallery to BDD-Ops tree node", () => dragTo(page, el("TSLA.csv link in gallery"), el("BDD-Ops tree node")));
      await session.step(28, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(29, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(30, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(31, "When user clicks on \"Files > Demo\" tree node", () => clickOn(page, el("\"Files > Demo\" tree node")));
      await session.step(32, "And user drags acidiq.csv link in gallery to BDD-Ops tree node", () => dragTo(page, el("acidiq.csv link in gallery"), el("BDD-Ops tree node")));
      await session.step(33, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(34, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(35, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(36, "When user double-clicks on BDD-Ops tree node", () => doubleClickOn(page, el("BDD-Ops tree node")));
      await session.step(37, "Then the \"BDD-Ops\" view should be current", () => viewIsCurrent(page, "BDD-Ops"));
      await session.step(38, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(39, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A file offers open, rename and delete", async () => {
      await session.step(42, "When user opens the context menu of TSLA.csv link in gallery", () => openContextMenu(page, el("TSLA.csv link in gallery")));
      await session.step(43, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
      await session.step(44, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(45, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
      await session.step(46, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("A cancelled rename changes nothing", async () => {
      await session.step(49, "When user picks \"Rename...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Rename...", el("acidiq.csv link in gallery")));
      await session.step(50, "Then Rename dialog should be visible", () => shouldBe(page, el("Rename dialog"), "visible"));
      await session.step(51, "When user enters \"BDD-should-not-appear\" into \"File name\" input in Rename dialog", () => enterInto(page, "BDD-should-not-appear", el("\"File name\" input in Rename dialog")));
      await session.step(52, "And user clicks on CANCEL button in Rename dialog", () => clickOn(page, el("CANCEL button in Rename dialog")));
      await session.step(53, "Then Rename dialog should be hidden", () => shouldBe(page, el("Rename dialog"), "hidden"));
      await session.step(54, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(55, "And BDD-should-not-appear link in gallery should be absent", () => shouldBe(page, el("BDD-should-not-appear link in gallery"), "absent"));
    });
    await run.scenario("A file is renamed", async () => {
      await session.step(58, "When user picks \"Rename...\" from the context menu of TSLA.csv link in gallery", () => pickFromContextMenu(page, "Rename...", el("TSLA.csv link in gallery")));
      await session.step(59, "And user enters \"BDD-Ops-renamed\" into \"File name\" input in Rename dialog", () => enterInto(page, "BDD-Ops-renamed", el("\"File name\" input in Rename dialog")));
      await session.step(60, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
      await session.step(61, "Then Rename dialog should be hidden", () => shouldBe(page, el("Rename dialog"), "hidden"));
      await session.step(62, "And BDD-Ops-renamed link in gallery should be visible", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "visible"));
      await session.step(63, "And TSLA.csv link in gallery should be absent", () => shouldBe(page, el("TSLA.csv link in gallery"), "absent"));
    });
    await run.scenario("A cancelled delete keeps the file", async () => {
      await session.step(66, "When user picks \"Delete...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Delete...", el("acidiq.csv link in gallery")));
      await session.step(67, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(68, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(69, "Then \"Are you sure?\" dialog should be hidden", () => shouldBe(page, el("\"Are you sure?\" dialog"), "hidden"));
      await session.step(70, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy in another space survives the original being deleted", async () => {
      await session.step(73, "When user drags acidiq.csv link in gallery to BDD-Ops-Copy tree node", () => dragTo(page, el("acidiq.csv link in gallery"), el("BDD-Ops-Copy tree node")));
      await session.step(74, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(75, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(76, "When user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(77, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(78, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(79, "When user double-clicks on BDD-Ops-Copy tree node", () => doubleClickOn(page, el("BDD-Ops-Copy tree node")));
      await session.step(80, "Then the \"BDD-Ops-Copy\" view should be current", () => viewIsCurrent(page, "BDD-Ops-Copy"));
      await session.step(81, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(82, "When user double-clicks on BDD-Ops tree node", () => doubleClickOn(page, el("BDD-Ops tree node")));
      await session.step(83, "And user picks \"Delete...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Delete...", el("acidiq.csv link in gallery")));
      await session.step(84, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(85, "Then acidiq.csv link in gallery should be absent", () => shouldBe(page, el("acidiq.csv link in gallery"), "absent"));
      await session.step(86, "And BDD-Ops-renamed link in gallery should be visible", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "visible"));
      await session.step(87, "When user double-clicks on BDD-Ops-Copy tree node", () => doubleClickOn(page, el("BDD-Ops-Copy tree node")));
      await session.step(88, "Then the \"BDD-Ops-Copy\" view should be current", () => viewIsCurrent(page, "BDD-Ops-Copy"));
      await session.step(89, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A file opens as a table", async () => {
      await session.step(92, "When user double-clicks on acidiq.csv link in gallery", () => doubleClickOn(page, el("acidiq.csv link in gallery")));
      await session.step(93, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    });
    run.finish();
  });
});
