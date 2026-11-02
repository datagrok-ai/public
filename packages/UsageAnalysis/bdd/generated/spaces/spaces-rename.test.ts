/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-rename.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space]
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
import {clearField, clickOn, doubleClickOn, enterInto, expand, isExpanded, shouldBe, shouldHaveValue, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, childSpaceOnServer, noSpaceOnServer, spaceOnServer, spacesOnServer} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Renaming a space", () => {
  const session = feature(test, "features/spaces/spaces-rename.feature", import.meta.url);
  test("Renaming a space", {tag: ["@journey", "@spaces", "@realizes:views.space", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(22, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(23, "And no space named \"BDD-Ren, BDD-Ren-New, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew\" is on the server", () => noSpaceOnServer(page, "BDD-Ren, BDD-Ren-New, BDD-Ren-Parent, BDD-Ren-Child, BDD-Ren-ChildNew"));
    await run.scenario("The Rename dialog opens on the current name", async () => {
      await session.step(26, "Given a space named \"BDD-Ren\" is on the server", () => spaceOnServer(page, "BDD-Ren"));
      await session.step(27, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(28, "When user picks \"Rename...\" from the context menu of Spaces---BDD-Ren tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("Spaces---BDD-Ren tree node inside browse tree")));
      await session.step(29, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(30, "And Name input in Rename project dialog should have value \"BDD-Ren\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), "BDD-Ren"));
    });
    await run.scenario("Cancelling the rename keeps the old name", async () => {
      await session.step(33, "When user enters \"BDD-Ren-New\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-New", el("Name input in Rename project dialog")));
      await session.step(34, "And user clicks on CANCEL button in Rename project dialog", () => clickOn(page, el("CANCEL button in Rename project dialog")));
      await session.step(35, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
      await session.step(36, "And 1 space named \"BDD-Ren\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren"));
      await session.step(37, "And 0 spaces named \"BDD-Ren-New\" should be on the server", () => spacesOnServer(page, 0, "BDD-Ren-New"));
      await session.step(38, "And Spaces---BDD-Ren tree node inside browse tree should be visible", () => shouldBe(page, el("Spaces---BDD-Ren tree node inside browse tree"), "visible"));
    });
    await run.scenario("A rename reaches the server and the tree", async () => {
      await session.step(41, "When user picks \"Rename...\" from the context menu of Spaces---BDD-Ren tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el("Spaces---BDD-Ren tree node inside browse tree")));
      await session.step(42, "And user enters \"BDD-Ren-New\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-New", el("Name input in Rename project dialog")));
      await session.step(43, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(44, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
      await session.step(45, "And 1 space named \"BDD-Ren-New\" should be on the server", () => spacesOnServer(page, 1, "BDD-Ren-New"));
      await session.step(46, "And 0 spaces named \"BDD-Ren\" should be on the server", () => spacesOnServer(page, 0, "BDD-Ren"));
      await session.step(47, "And BDD-Ren-New tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren-New tree node inside browse tree"), "visible"));
      await session.step(48, "And BDD-Ren tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Ren tree node inside browse tree"), "absent"));
    });
    await run.scenario("A child space is renamed from its card in the parent", async () => {
      await session.step(51, "Given a space named \"BDD-Ren-Parent\" is on the server", () => spaceOnServer(page, "BDD-Ren-Parent"));
      await session.step(52, "Given a space named \"BDD-Ren-Child\" under \"BDD-Ren-Parent\" is on the server", () => childSpaceOnServer(page, "BDD-Ren-Child", "BDD-Ren-Parent"));
      await session.step(53, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(54, "When user double-clicks on Spaces---BDD-Ren-Parent tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-Ren-Parent tree node inside browse tree")));
      await session.step(55, "Then BDD-Ren-Child link in gallery should be visible", () => shouldBe(page, el("BDD-Ren-Child link in gallery"), "visible"));
      await session.step(56, "When user picks \"Rename...\" from the context menu of BDD-Ren-Child link in gallery", () => pickFromContextMenu(page, "Rename...", el("BDD-Ren-Child link in gallery")));
      await session.step(57, "Then Name input in Rename project dialog should have value \"BDD-Ren-Child\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), "BDD-Ren-Child"));
      await session.step(58, "When user enters \"BDD-Ren-ChildNew\" into Name input in Rename project dialog", () => enterInto(page, "BDD-Ren-ChildNew", el("Name input in Rename project dialog")));
      await session.step(59, "And user clicks on OK button in Rename project dialog", () => clickOn(page, el("OK button in Rename project dialog")));
      await session.step(60, "Then BDD-Ren-ChildNew link in gallery should be visible", () => shouldBe(page, el("BDD-Ren-ChildNew link in gallery"), "visible"));
      await session.step(61, "And BDD-Ren-Child link in gallery should be absent", () => shouldBe(page, el("BDD-Ren-Child link in gallery"), "absent"));
      await session.step(62, "When user expands Spaces---BDD-Ren-Parent tree node inside browse tree", () => expand(page, el("Spaces---BDD-Ren-Parent tree node inside browse tree")));
      await session.step(63, "Then BDD-Ren-ChildNew tree node inside browse tree should be visible", () => shouldBe(page, el("BDD-Ren-ChildNew tree node inside browse tree"), "visible"));
      await session.step(64, "And BDD-Ren-Child tree node inside browse tree should be absent", () => shouldBe(page, el("BDD-Ren-Child tree node inside browse tree"), "absent"));
    });
    await run.scenario("An emptied name disables OK until a name is typed again", async () => {
      await session.step(68, "Given a space named \"BDDEdgeRen{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeRen{time}")));
      await session.step(69, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(70, "When user picks \"Rename...\" from the context menu of Spaces---BDDEdgeRen{time} tree node inside browse tree", () => pickFromContextMenu(page, "Rename...", el(session.text("Spaces---BDDEdgeRen{time} tree node inside browse tree"))));
      await session.step(71, "Then Rename project dialog should be visible", () => shouldBe(page, el("Rename project dialog"), "visible"));
      await session.step(72, "And Name input in Rename project dialog should have value \"BDDEdgeRen{time}\"", () => shouldHaveValue(page, el("Name input in Rename project dialog"), session.text("BDDEdgeRen{time}")));
      await session.step(73, "And OK button in Rename project dialog should be enabled", () => shouldBe(page, el("OK button in Rename project dialog"), "enabled"));
      await session.step(74, "When user clears Name input in Rename project dialog", () => clearField(page, el("Name input in Rename project dialog")));
      await session.step(75, "Then Name input in Rename project dialog should be invalid", () => shouldBe(page, el("Name input in Rename project dialog"), "invalid"));
      await session.step(76, "And OK button in Rename project dialog should be disabled", () => shouldBe(page, el("OK button in Rename project dialog"), "disabled"));
      await session.step(77, "When user types \"x\" into Name input in Rename project dialog", () => typeInto(page, "x", el("Name input in Rename project dialog")));
      await session.step(78, "Then OK button in Rename project dialog should be enabled", () => shouldBe(page, el("OK button in Rename project dialog"), "enabled"));
      await session.step(79, "And no errors should have been logged", () => noErrors(page));
      await session.step(80, "When user clicks on CANCEL button in Rename project dialog", () => clickOn(page, el("CANCEL button in Rename project dialog")));
      await session.step(81, "Then Rename project dialog should be hidden", () => shouldBe(page, el("Rename project dialog"), "hidden"));
      await session.step(82, "And 1 space named \"BDDEdgeRen{time}\" should be on the server", () => spacesOnServer(page, 1, session.text("BDDEdgeRen{time}")));
    });
    run.finish();
  });
});
