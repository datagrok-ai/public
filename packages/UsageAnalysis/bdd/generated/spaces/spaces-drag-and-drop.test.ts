/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-drag-and-drop.feature
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
import {clearField, clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, openDataset, saveAsProject, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInSpace} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Putting files into a space by dragging them", () => {
  const session = feature(test, "features/spaces/spaces-drag-and-drop.feature", import.meta.url);
  test("Putting files into a space by dragging them", {tag: ["@journey", "@spaces", "@realizes:views.space", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 10, page);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "And a space named \"BDD-DnD\" is on the server", () => spaceOnServer(page, "BDD-DnD"));
    await session.step(32, "And a space named \"BDD-DnD-Src\" is on the server", () => spaceOnServer(page, "BDD-DnD-Src"));
    await session.step(33, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The demo files", async () => {
      await session.step(36, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(37, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(38, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(39, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The dialog offers Link, Copy and Move, and names the target", async () => {
      await session.step(42, "When user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(43, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(44, "And Move entity dialog should contain text \"BDDDnD\"", () => shouldContainText(page, el("Move entity dialog"), "BDDDnD"));
      await session.step(45, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
    });
    await run.scenario("Cancelling leaves the space empty", async () => {
      await session.step(48, "When user clicks on CANCEL button in Move entity dialog", () => clickOn(page, el("CANCEL button in Move entity dialog")));
      await session.step(49, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(50, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(51, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(52, "And demog.csv link in gallery should be absent", () => shouldBe(page, el("demog.csv link in gallery"), "absent"));
    });
    await run.scenario("A copied file lands in the space", async () => {
      await session.step(55, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(56, "And user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(57, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(58, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(59, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(60, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(61, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(62, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file's details show in the context panel", async () => {
      await session.step(65, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(66, "When user clicks on demog.csv link in gallery", () => clickOn(page, el("demog.csv link in gallery")));
      await session.step(67, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
      await session.step(68, "And context panel should contain text \"demog\"", () => shouldContainText(page, el("context panel"), "demog"));
    });
    await run.scenario("A second file joins the first", async () => {
      await session.step(71, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(72, "And user drags TSLA.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("TSLA.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(73, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(74, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(75, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(76, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(77, "Then demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(78, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy is made in the source space to move later", async () => {
      await session.step(81, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(82, "And user drags beer.csv link in gallery to Spaces---BDD-DnD-Src tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(83, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(84, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(85, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(86, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(87, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(88, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("Moving takes the file out of the source space", async () => {
      await session.step(91, "When user drags beer.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(92, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(93, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(94, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(95, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(96, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(97, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(98, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(99, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
      await session.step(100, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(101, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(102, "And beer.csv link in gallery should be absent", () => shouldBe(page, el("beer.csv link in gallery"), "absent"));
    });
    await run.scenario("The demo files kept their originals", async () => {
      await session.step(105, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(106, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(107, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(108, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(109, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("A linked project follows its original moved to a third space", async () => {
      await session.step(113, "Given a space named \"BDDEdgeSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSrc{time}")));
      await session.step(114, "And a space named \"BDDEdgeDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeDst{time}")));
      await session.step(115, "And a space named \"BDDEdgeMoveTo{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeMoveTo{time}")));
      await session.step(116, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeSrc{time}\" as \"edgecars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeSrc{time}"), session.text("edgecars{time}.csv")));
      await session.step(117, "And no project named \"BDDEdgeMoved{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDEdgeMoved{time}")));
      await session.step(118, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(119, "When user saves the current view as project \"BDDEdgeMoved{time}\"", () => saveAsProject(page, session.text("BDDEdgeMoved{time}")));
      await session.step(120, "And user closes all views", () => closeAllViews(page));
      await session.step(121, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(122, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(123, "And user enters \"BDDEdgeMoved{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeMoved{time}"), el("gallery search")));
      await session.step(124, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(125, "And user picks \"Move to Space...\" from the context menu of BDDEdgeMoved{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDEdgeMoved{time} gallery card"))));
      await session.step(126, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
      await session.step(127, "When user selects \"BDDEdgeSrc{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDEdgeSrc{time}"), el("Space input in \"Move to space\" dialog")));
      await session.step(128, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
      await session.step(129, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(130, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(131, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(132, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
      await session.step(133, "Then \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(134, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(135, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(136, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(137, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(138, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(139, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeMoveTo{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
      await session.step(140, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(141, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(142, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(143, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(144, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(145, "When user double-clicks on BDDEdgeMoveTo{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
      await session.step(146, "Then the \"BDDEdgeMoveTo{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeMoveTo{time}")));
      await session.step(147, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(148, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
      await session.step(149, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
      await session.step(150, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
      await session.step(151, "And \"BDDEdgeMoved{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "absent"));
      await session.step(152, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(153, "When user double-clicks on BDDEdgeDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(154, "Then the \"BDDEdgeDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeDst{time}")));
      await session.step(155, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(156, "When user clicks on \"BDDEdgeMoved{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery"))));
      await session.step(157, "Then the context panel should show \"BDDEdgeMoved{time}\"", () => contextPanelShows(page, session.text("BDDEdgeMoved{time}")));
      await session.step(158, "And no errors should have been logged", () => noErrors(page));
      await session.step(159, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    run.finish();
  });
});
