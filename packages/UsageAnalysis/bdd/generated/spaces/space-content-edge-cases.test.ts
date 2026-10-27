/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/spaces/space-content-edge-cases.feature
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
import {clearField, clickOn, doubleClickOn, dragTo, enterInto, isExpanded, selectIn, shouldBe, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, closeAllViews, contextPanelOpen, contextPanelShows, dialogCloses, noProjectOnServer, openDataset, saveAsProject, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {fixtureInSpace} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors, pickFromContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Space content edge cases: a file linked between spaces, a link whose original moved, a package script linked", () => {
  const session = feature(test, "features/spaces/space-content-edge-cases.feature", import.meta.url);
  test("A file linked from one space into another shows in the target space", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space", "@known-failure:45", "@realizes:core#1042"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And a space named \"BDDEdgeSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSrc{time}")));
    await session.step(30, "And a space named \"BDDEdgeDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeDst{time}")));
    await session.step(31, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeSrc{time}\" as \"edgecars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeSrc{time}"), session.text("edgecars{time}.csv")));
    await session.step(32, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(33, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
    await session.step(34, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
    await session.step(35, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
    await knownFailure(async () => {
      await session.step(41, "When user drags \"edgecars{time}.csv\" link in gallery to BDDEdgeDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(42, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(43, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(44, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(45, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(46, "And no error or warning balloon should have been shown", () => noBalloons(page));
      await session.step(47, "When user double-clicks on BDDEdgeDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
      await session.step(48, "Then the \"BDDEdgeDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeDst{time}")));
      await session.step(49, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
    }, 45);
  });
  test("A package script dropped on a space is linked with no dialog", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And a space named \"BDDEdgeSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSrc{time}")));
    await session.step(30, "And a space named \"BDDEdgeDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeDst{time}")));
    await session.step(31, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeSrc{time}\" as \"edgecars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeSrc{time}"), session.text("edgecars{time}.csv")));
    await session.step(32, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(33, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
    await session.step(34, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
    await session.step(35, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
    await session.step(52, "Given Platform tree node inside browse tree is expanded", () => isExpanded(page, el("Platform tree node inside browse tree")));
    await session.step(53, "And Platform---Functions tree node inside browse tree is expanded", () => isExpanded(page, el("Platform---Functions tree node inside browse tree")));
    await session.step(54, "When user clicks on Platform---Functions---Scripts tree node inside browse tree", () => clickOn(page, el("Platform---Functions---Scripts tree node inside browse tree")));
    await session.step(55, "And user enters \"Protein-Ligand\" into gallery search", () => enterInto(page, "Protein-Ligand", el("gallery search")));
    await session.step(56, "Then \"Protein-Ligand Interaction Diagram\" gallery card should be visible", () => shouldBe(page, el("\"Protein-Ligand Interaction Diagram\" gallery card"), "visible"));
    await session.step(57, "When user drags \"Protein-Ligand Interaction Diagram\" gallery card to BDDEdgeDst{time} tree node inside browse tree", () => dragTo(page, el("\"Protein-Ligand Interaction Diagram\" gallery card"), el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
    await session.step(58, "Then Move entity dialog should be absent", () => shouldBe(page, el("Move entity dialog"), "absent"));
    await session.step(59, "When user clears gallery search", () => clearField(page, el("gallery search")));
    await session.step(60, "And user double-clicks on BDDEdgeDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
    await session.step(61, "Then the \"BDDEdgeDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeDst{time}")));
    await session.step(62, "And \"Protein-Ligand Interaction Diagram\" link in gallery should be visible", () => shouldBe(page, el("\"Protein-Ligand Interaction Diagram\" link in gallery"), "visible"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
    await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("A linked project still opens its details after its original moved to another space", {tag: ["@serial", "@spaces", "@edge-cases", "@realizes:views.space"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(29, "And a space named \"BDDEdgeSrc{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeSrc{time}")));
    await session.step(30, "And a space named \"BDDEdgeDst{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeDst{time}")));
    await session.step(31, "And the \"fixtures/cars-small.csv\" file of the project is in the space \"BDDEdgeSrc{time}\" as \"edgecars{time}.csv\"", () => fixtureInSpace(page, "fixtures/cars-small.csv", session.text("BDDEdgeSrc{time}"), session.text("edgecars{time}.csv")));
    await session.step(32, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(33, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
    await session.step(34, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
    await session.step(35, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
    await session.step(67, "Given no project named \"BDDEdgeMoved{time}\" is on the server", () => noProjectOnServer(page, session.text("BDDEdgeMoved{time}")));
    await session.step(68, "And a space named \"BDDEdgeMoveTo{time}\" is on the server", () => spaceOnServer(page, session.text("BDDEdgeMoveTo{time}")));
    await session.step(69, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(70, "When user saves the current view as project \"BDDEdgeMoved{time}\"", () => saveAsProject(page, session.text("BDDEdgeMoved{time}")));
    await session.step(71, "And user closes all views", () => closeAllViews(page));
    await session.step(72, "Given the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "When user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
    await session.step(74, "And user enters \"BDDEdgeMoved{time}\" into gallery search", () => enterInto(page, session.text("BDDEdgeMoved{time}"), el("gallery search")));
    await session.step(75, "And user clicks on \"Refresh\" icon inside gallery toolbar", () => clickOn(page, el("\"Refresh\" icon inside gallery toolbar")));
    await session.step(76, "And user picks \"Move to Space...\" from the context menu of BDDEdgeMoved{time} gallery card", () => pickFromContextMenu(page, "Move to Space...", el(session.text("BDDEdgeMoved{time} gallery card"))));
    await session.step(77, "Then \"Move to space\" dialog should be visible", () => shouldBe(page, el("\"Move to space\" dialog"), "visible"));
    await session.step(78, "When user selects \"BDDEdgeSrc{time}\" in Space input in \"Move to space\" dialog", () => selectIn(page, session.text("BDDEdgeSrc{time}"), el("Space input in \"Move to space\" dialog")));
    await session.step(79, "And user clicks on OK button in \"Move to space\" dialog", () => clickOn(page, el("OK button in \"Move to space\" dialog")));
    await session.step(80, "Then the \"Move to space\" dialog should close", () => dialogCloses(page, "Move to space"));
    await session.step(81, "Given Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await session.step(82, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
    await session.step(83, "Then \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
    await session.step(84, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeDst{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
    await session.step(85, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
    await session.step(86, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
    await session.step(87, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
    await session.step(88, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
    await session.step(89, "When user drags \"BDDEdgeMoved{time}\" link in gallery to BDDEdgeMoveTo{time} tree node inside browse tree", () => dragTo(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
    await session.step(90, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
    await session.step(91, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
    await session.step(92, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
    await session.step(93, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
    await session.step(94, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
    await session.step(95, "When user double-clicks on BDDEdgeMoveTo{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeMoveTo{time} tree node inside browse tree"))));
    await session.step(96, "Then the \"BDDEdgeMoveTo{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeMoveTo{time}")));
    await session.step(97, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
    await session.step(98, "When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeSrc{time} tree node inside browse tree"))));
    await session.step(99, "Then the \"BDDEdgeSrc{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeSrc{time}")));
    await session.step(100, "And \"edgecars{time}.csv\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"edgecars{time}.csv\" link in gallery")), "visible"));
    await session.step(101, "And \"BDDEdgeMoved{time}\" link in gallery should be absent", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "absent"));
    await session.step(102, "Given the context panel is open", () => contextPanelOpen(page));
    await session.step(103, "When user double-clicks on BDDEdgeDst{time} tree node inside browse tree", () => doubleClickOn(page, el(session.text("BDDEdgeDst{time} tree node inside browse tree"))));
    await session.step(104, "Then the \"BDDEdgeDst{time}\" view should be current", () => viewIsCurrent(page, session.text("BDDEdgeDst{time}")));
    await session.step(105, "And \"BDDEdgeMoved{time}\" link in gallery should be visible", () => shouldBe(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery")), "visible"));
    await session.step(106, "When user clicks on \"BDDEdgeMoved{time}\" link in gallery", () => clickOn(page, el(session.text("\"BDDEdgeMoved{time}\" link in gallery"))));
    await session.step(107, "Then the context panel should show \"BDDEdgeMoved{time}\"", () => contextPanelShows(page, session.text("BDDEdgeMoved{time}")));
    await session.step(108, "And no errors should have been logged", () => noErrors(page));
    await session.step(109, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
