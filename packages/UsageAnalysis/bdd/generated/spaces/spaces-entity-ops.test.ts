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
import {noSpaceOnServer, spacesOnServer} from '../../bindings/spaces.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clearField, clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, dialogCloses, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {closeContextMenu, menuLists, openContextMenu, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Working with what a space holds", () => {
  const session = feature(test, "features/spaces/spaces-entity-ops.feature", import.meta.url);
  test("Working with what a space holds", {tag: ["@journey", "@spaces", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(14, "And no space named \"BDD-Ops, BDD-Ops-Copy\" is on the server", () => noSpaceOnServer(page, "BDD-Ops, BDD-Ops-Copy"));
    await run.scenario("A space with two files", async () => {
      await session.step(17, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(18, "And user enters \"BDD-Ops\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ops", el("Name input in Create Space dialog")));
      await session.step(19, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(20, "Then 1 space named \"BDD-Ops\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ops"));
      await session.step(21, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(22, "When user picks \"Create Space...\" from the context menu of Spaces tree node inside browse tree", () => pickFromContextMenu(page, "Create Space...", el("Spaces tree node inside browse tree")));
      await session.step(23, "And user enters \"BDD-Ops-Copy\" into Name input in Create Space dialog", () => enterInto(page, "BDD-Ops-Copy", el("Name input in Create Space dialog")));
      await session.step(24, "And user clicks on OK button in Create Space dialog", () => clickOn(page, el("OK button in Create Space dialog")));
      await session.step(25, "Then 1 space named \"BDD-Ops-Copy\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ops-Copy"));
      await session.step(26, "And the \"Create Space\" dialog should close", () => dialogCloses(page, "Create Space"));
      await session.step(27, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(28, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(29, "And user drags TSLA.csv link in gallery to BDD-Ops tree node inside browse tree", () => dragTo(page, el("TSLA.csv link in gallery"), el("BDD-Ops tree node inside browse tree")));
      await session.step(30, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(31, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(32, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(33, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(34, "And user drags acidiq.csv link in gallery to BDD-Ops tree node inside browse tree", () => dragTo(page, el("acidiq.csv link in gallery"), el("BDD-Ops tree node inside browse tree")));
      await session.step(35, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(36, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(37, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(38, "When user double-clicks on BDD-Ops tree node inside browse tree", () => doubleClickOn(page, el("BDD-Ops tree node inside browse tree")));
      await session.step(39, "Then the \"BDD-Ops\" view should be current", () => viewIsCurrent(page, "BDD-Ops"));
      await session.step(40, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(41, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A file offers open, rename and delete", async () => {
      await session.step(44, "When user opens the context menu of TSLA.csv link in gallery", () => openContextMenu(page, el("TSLA.csv link in gallery")));
      await session.step(45, "Then the open menu should list \"Open\"", () => menuLists(page, "Open"));
      await session.step(46, "And the open menu should list \"Rename...\"", () => menuLists(page, "Rename..."));
      await session.step(47, "And the open menu should list \"Delete...\"", () => menuLists(page, "Delete..."));
      await session.step(48, "When user closes the context menu", () => closeContextMenu(page));
    });
    await run.scenario("A cancelled rename changes nothing", async () => {
      await session.step(51, "When user picks \"Rename...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Rename...", el("acidiq.csv link in gallery")));
      await session.step(52, "Then Rename dialog should be visible", () => shouldBe(page, el("Rename dialog"), "visible"));
      await session.step(53, "When user enters \"BDD-should-not-appear\" into \"File name\" input in Rename dialog", () => enterInto(page, "BDD-should-not-appear", el("\"File name\" input in Rename dialog")));
      await session.step(54, "And user clicks on CANCEL button in Rename dialog", () => clickOn(page, el("CANCEL button in Rename dialog")));
      await session.step(55, "Then Rename dialog should be hidden", () => shouldBe(page, el("Rename dialog"), "hidden"));
      await session.step(56, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(57, "And BDD-should-not-appear link in gallery should be absent", () => shouldBe(page, el("BDD-should-not-appear link in gallery"), "absent"));
    });
    await run.scenario("A file is renamed", async () => {
      await session.step(60, "When user picks \"Rename...\" from the context menu of TSLA.csv link in gallery", () => pickFromContextMenu(page, "Rename...", el("TSLA.csv link in gallery")));
      await session.step(61, "And user enters \"BDD-Ops-renamed\" into \"File name\" input in Rename dialog", () => enterInto(page, "BDD-Ops-renamed", el("\"File name\" input in Rename dialog")));
      await session.step(62, "And user clicks on OK button in Rename dialog", () => clickOn(page, el("OK button in Rename dialog")));
      await session.step(63, "Then Rename dialog should be hidden", () => shouldBe(page, el("Rename dialog"), "hidden"));
      await session.step(64, "And BDD-Ops-renamed link in gallery should be visible", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "visible"));
      await session.step(65, "And TSLA.csv link in gallery should be absent", () => shouldBe(page, el("TSLA.csv link in gallery"), "absent"));
    });
    await run.scenario("The search inside a space filters what it holds", async () => {
      await session.step(68, "When user enters \"acidiq\" into space search", () => enterInto(page, "acidiq", el("space search")));
      await session.step(69, "Then acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(70, "And BDD-Ops-renamed link in gallery should be absent", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "absent"));
      await session.step(71, "When user enters \"aci\" into space search", () => enterInto(page, "aci", el("space search")));
      await session.step(72, "Then acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(73, "And BDD-Ops-renamed link in gallery should be absent", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "absent"));
      await session.step(74, "When user enters \"zzz-no-such-file\" into space search", () => enterInto(page, "zzz-no-such-file", el("space search")));
      await session.step(75, "Then acidiq.csv link in gallery should be absent", () => shouldBe(page, el("acidiq.csv link in gallery"), "absent"));
      await session.step(76, "And BDD-Ops-renamed link in gallery should be absent", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "absent"));
      await session.step(77, "When user clears space search", () => clearField(page, el("space search")));
      await session.step(78, "Then acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(79, "And BDD-Ops-renamed link in gallery should be visible", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "visible"));
    });
    await run.scenario("A cancelled delete keeps the file", async () => {
      await session.step(82, "When user picks \"Delete...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Delete...", el("acidiq.csv link in gallery")));
      await session.step(83, "Then \"Are you sure?\" dialog should be visible", () => shouldBe(page, el("\"Are you sure?\" dialog"), "visible"));
      await session.step(84, "When user clicks on CANCEL button in \"Are you sure?\" dialog", () => clickOn(page, el("CANCEL button in \"Are you sure?\" dialog")));
      await session.step(85, "Then \"Are you sure?\" dialog should be hidden", () => shouldBe(page, el("\"Are you sure?\" dialog"), "hidden"));
      await session.step(86, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy in another space survives the original being deleted", async () => {
      await session.step(89, "When user drags acidiq.csv link in gallery to BDD-Ops-Copy tree node inside browse tree", () => dragTo(page, el("acidiq.csv link in gallery"), el("BDD-Ops-Copy tree node inside browse tree")));
      await session.step(90, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(91, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(92, "When user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(93, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(94, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(95, "When user double-clicks on BDD-Ops-Copy tree node inside browse tree", () => doubleClickOn(page, el("BDD-Ops-Copy tree node inside browse tree")));
      await session.step(96, "Then the \"BDD-Ops-Copy\" view should be current", () => viewIsCurrent(page, "BDD-Ops-Copy"));
      await session.step(97, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
      await session.step(98, "When user double-clicks on BDD-Ops tree node inside browse tree", () => doubleClickOn(page, el("BDD-Ops tree node inside browse tree")));
      await session.step(99, "And user picks \"Delete...\" from the context menu of acidiq.csv link in gallery", () => pickFromContextMenu(page, "Delete...", el("acidiq.csv link in gallery")));
      await session.step(100, "And user clicks on DELETE button in \"Are you sure?\" dialog", () => clickOn(page, el("DELETE button in \"Are you sure?\" dialog")));
      await session.step(101, "Then acidiq.csv link in gallery should be absent", () => shouldBe(page, el("acidiq.csv link in gallery"), "absent"));
      await session.step(102, "And BDD-Ops-renamed link in gallery should be visible", () => shouldBe(page, el("BDD-Ops-renamed link in gallery"), "visible"));
      await session.step(103, "When user double-clicks on BDD-Ops-Copy tree node inside browse tree", () => doubleClickOn(page, el("BDD-Ops-Copy tree node inside browse tree")));
      await session.step(104, "Then the \"BDD-Ops-Copy\" view should be current", () => viewIsCurrent(page, "BDD-Ops-Copy"));
      await session.step(105, "And acidiq.csv link in gallery should be visible", () => shouldBe(page, el("acidiq.csv link in gallery"), "visible"));
    });
    await run.scenario("A single click previews the file in place", async () => {
      await session.step(108, "When user clicks on acidiq.csv link in gallery", () => clickOn(page, el("acidiq.csv link in gallery")));
      await session.step(109, "Then \"Toggle entity preview\" icon should be visible", () => shouldBe(page, el("\"Toggle entity preview\" icon"), "visible"));
      await session.step(110, "And grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    });
    await run.scenario("A file opens as a table", async () => {
      await session.step(113, "When user double-clicks on acidiq.csv link in gallery", () => doubleClickOn(page, el("acidiq.csv link in gallery")));
      await session.step(114, "Then grid should be visible", () => shouldBe(page, el("grid"), "visible"));
    });
    run.finish();
  });
});
