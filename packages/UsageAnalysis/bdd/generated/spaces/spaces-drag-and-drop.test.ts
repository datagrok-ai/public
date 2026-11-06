/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/spaces-drag-and-drop.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.space, core#1042]
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
  test("Putting files into a space by dragging them", {tag: ["@journey", "@spaces", "@realizes:views.space", "@edge-cases", "@realizes:core#1042"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 11, page);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "And a space named \"BDD-DnD\" is on the server", () => spaceOnServer(page, "BDD-DnD"));
    await session.step(31, "And a space named \"BDD-DnD-Src\" is on the server", () => spaceOnServer(page, "BDD-DnD-Src"));
    await session.step(32, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The demo files", async () => {
      await session.step(35, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(36, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(37, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(38, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The dialog offers Link, Copy and Move, and names the target", async () => {
      await session.step(41, "When user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(42, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(43, "And Move entity dialog should contain text \"BDDDnD\"", () => shouldContainText(page, el("Move entity dialog"), "BDDDnD"));
      await session.step(44, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
    });
    await run.scenario("Cancelling leaves the space empty", async () => {
      await session.step(47, "When user clicks on CANCEL button in Move entity dialog", () => clickOn(page, el("CANCEL button in Move entity dialog")));
      await session.step(48, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(49, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(50, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(51, "And demog.csv link in gallery should be absent", () => shouldBe(page, el("demog.csv link in gallery"), "absent"));
    });
    await run.scenario("A copied file lands in the space", async () => {
      await session.step(54, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(55, "And user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(56, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(57, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(58, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(59, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(60, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(61, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file's details show in the context panel", async () => {
      await session.step(64, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(65, "When user clicks on demog.csv link in gallery", () => clickOn(page, el("demog.csv link in gallery")));
      await session.step(66, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
      await session.step(67, "And context panel should contain text \"demog\"", () => shouldContainText(page, el("context panel"), "demog"));
    });
    await run.scenario("A second file joins the first", async () => {
      await session.step(70, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(71, "And user drags TSLA.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("TSLA.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(72, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(73, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(74, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(75, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(76, "Then demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(77, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy is made in the source space to move later", async () => {
      await session.step(80, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(81, "And user drags beer.csv link in gallery to Spaces---BDD-DnD-Src tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(82, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(83, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(84, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(85, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(86, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(87, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("Moving takes the file out of the source space", async () => {
      await session.step(90, "When user drags beer.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(91, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(92, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(93, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(94, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(95, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(96, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(97, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(98, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
      await session.step(99, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(100, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(101, "And beer.csv link in gallery should be absent", () => shouldBe(page, el("beer.csv link in gallery"), "absent"));
    });
    await run.scenario("The demo files kept their originals", async () => {
      await session.step(104, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(105, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(106, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(107, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(108, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("A linked project follows its original moved to a third space", async () => {
      await session.step(112, "Given a space named \"BDDEdgeSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSrc{time}")));
      await session.step(113, "And a space named \"BDDEdgeDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeDst{time}")));
      await session.step(114, "And a space named \"BDDEdgeMoveTo{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeMoveTo{time}")));
      await session.step(115, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeSrc{time}\" as \"edgecars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeSrc{time}"), session.text("edgecars{time}.csv")));
      await session.step(116, "And no project named \"BDDEdgeMoved{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDEdgeMoved{time}")));
      await session.step(117, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(118, "When user saves the current view as project \"BDDEdgeMoved{time}\"", () => saveAsProject(page, session.text("BDDEdgeMoved{time}")));
      await session.step(119, "And user closes all views", () => closeAllViews(page));
      await session.step(120, "Given the browse panel is open", () => browsePanelOpen(page));
      await session.step(121, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(122, "And user enters \"BDDEdgeMoved{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeMoved{time}"), el("gallery search")));
      await session.step(123, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
      await session.step(124, "And user picks \"Move to Space...\" from the context menu of BDDEdgeMoved{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDEdgeMoved{time} gallery card"))));
      await session.step(125, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
      await session.step(126, "When user selects \"BDDEdgeSrc{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDEdgeSrc{time}"), el("Space input in \"Move to space\" dialog")));
      await session.step(127, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
      await session.step(128, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
      await session.step(129, "When user clears gallery search", () => clearField(page, el("gallery search")));
      await session.step(130, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(131, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
      await session.step(132, "Then \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(133, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(134, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(135, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(136, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(137, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(138, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeMoveTo{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
      await session.step(139, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(140, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(141, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(142, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(143, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(144, "When user double-clicks on BDDEdgeMoveTo{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
      await session.step(145, "Then the \"BDDEdgeMoveTo{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeMoveTo{time}")));
      await session.step(146, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(147, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
      await session.step(148, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
      await session.step(149, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
      await session.step(150, "And \"BDDEdgeMoved{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "absent"));
      await session.step(151, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(152, "When user double-clicks on BDDEdgeDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(153, "Then the \"BDDEdgeDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeDst{time}")));
      await session.step(154, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
      await session.step(155, "When user clicks on \"BDDEdgeMoved{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery"))));
      await session.step(156, "Then the context panel should show \"BDDEdgeMoved{time}\"", () => contextPanelShows(page, session.text("BDDEdgeMoved{time}")));
      await session.step(157, "And no errors should have been logged", () => noErrors(page));
      await session.step(158, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("A file linked from one space into another shows in the target space", async () => {
      await session.step(162, "Given a space named \"BDDEdgeFileSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeFileSrc{time}")));
      await session.step(163, "And a space named \"BDDEdgeFileDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeFileDst{time}")));
      await session.step(164, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeFileSrc{time}\" as \"edgelinkcars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeFileSrc{time}"), session.text("edgelinkcars{time}.csv")));
      await session.step(165, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
      await session.step(166, "When user double-clicks on BDDEdgeFileSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeFileSrc{time} tree node inside browse tree"))));
      await session.step(167, "Then the \"BDDEdgeFileSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeFileSrc{time}")));
      await session.step(168, "And \"edgelinkcars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgelinkcars{time}.csv\" link in gallery")), "visible"));
      await session.step(169, "When user drags \"edgelinkcars{time}.csv\" link in gallery to BDDEdgeFileDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"edgelinkcars{time}.csv\" link in gallery")), el(session.text("BDDEdgeFileDst{time} tree node inside browse tree"))));
      await session.step(170, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(171, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(172, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(173, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(174, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(175, "When user double-clicks on BDDEdgeFileDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeFileDst{time} tree node inside browse tree"))));
      await session.step(176, "Then the \"BDDEdgeFileDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeFileDst{time}")));
      await session.step(177, "And \"edgelinkcars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgelinkcars{time}.csv\" link in gallery")), "visible"));
      await session.step(178, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
