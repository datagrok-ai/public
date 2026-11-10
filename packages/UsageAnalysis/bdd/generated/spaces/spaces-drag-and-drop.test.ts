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
import {clickOn, doubleClickOn, dragTo, isExpanded, selectIn, shouldBe, shouldContainText, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, contextPanelOpen, contextPanelShows, spaceOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Putting files into a space by dragging them", () => {
  const session = feature(test, "features/spaces/spaces-drag-and-drop.feature", import.meta.url);
  test("Putting files into a space by dragging them", {tag: ["@journey", "@spaces", "@realizes:views.space", "@realizes:core#1042"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(24, "Given user is logged in", () => loggedIn(page));
    await session.step(25, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(26, "And a space named \"BDD-DnD\" is on the server", () => spaceOnServer(page, "BDD-DnD"));
    await session.step(27, "And a space named \"BDD-DnD-Src\" is on the server", () => spaceOnServer(page, "BDD-DnD-Src"));
    await session.step(28, "And Spaces tree node inside browse tree is expanded", () => isExpanded(page, el("Spaces tree node inside browse tree")));
    await run.scenario("The demo files", async () => {
      await session.step(31, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(32, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(33, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(34, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The dialog offers Link, Copy and Move, and names the target", async () => {
      await session.step(37, "When user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(38, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(39, "And Move entity dialog should contain text \"BDDDnD\"", () => shouldContainText(page, el("Move entity dialog"), "BDDDnD"));
      await session.step(40, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
    });
    await run.scenario("Cancelling leaves the space empty", async () => {
      await session.step(43, "When user clicks on CANCEL button in Move entity dialog", () => clickOn(page, el("CANCEL button in Move entity dialog")));
      await session.step(44, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(45, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(46, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(47, "And demog.csv link in gallery should be absent", () => shouldBe(page, el("demog.csv link in gallery"), "absent"));
    });
    await run.scenario("A copied file lands in the space", async () => {
      await session.step(50, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(51, "And user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("demog.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(52, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(53, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(54, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(55, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(56, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(57, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
    });
    await run.scenario("The file's details show in the context panel", async () => {
      await session.step(60, "Given the context panel is open", () => contextPanelOpen(page));
      await session.step(61, "When user clicks on demog.csv link in gallery", () => clickOn(page, el("demog.csv link in gallery")));
      await session.step(62, "Then the context panel should show \"demog.csv\"", () => contextPanelShows(page, "demog.csv"));
      await session.step(63, "And context panel should contain text \"demog\"", () => shouldContainText(page, el("context panel"), "demog"));
    });
    await run.scenario("A second file joins the first", async () => {
      await session.step(66, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(67, "And user drags TSLA.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("TSLA.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(68, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(69, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(70, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(71, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(72, "Then demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(73, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
    });
    await run.scenario("A copy is made in the source space to move later", async () => {
      await session.step(76, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(77, "And user drags beer.csv link in gallery to Spaces---BDD-DnD-Src tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(78, "And user selects \"Copy\" in Move entity dialog", () => selectIn(page, "Copy", el("Move entity dialog")));
      await session.step(79, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(80, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(81, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(82, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(83, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    await run.scenario("Moving takes the file out of the source space, and a link puts it back there", async () => {
      await session.step(87, "When user drags beer.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(88, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(89, "And choice input in Move entity dialog should be visible", () => shouldBe(page, el("choice input in Move entity dialog"), "visible"));
      await session.step(90, "When user selects \"Move\" in Move entity dialog", () => selectIn(page, "Move", el("Move entity dialog")));
      await session.step(91, "And user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(92, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(93, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(94, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(95, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
      await session.step(96, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(97, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(98, "And beer.csv link in gallery should be absent", () => shouldBe(page, el("beer.csv link in gallery"), "absent"));
      await session.step(99, "When user double-clicks on Spaces---BDD-DnD tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD tree node inside browse tree")));
      await session.step(100, "Then the \"BDD-DnD\" view should be current", () => viewIsCurrent(page, "BDD-DnD"));
      await session.step(101, "When user drags beer.csv link in gallery to Spaces---BDD-DnD-Src tree node inside browse tree", () => dragTo(page, el("beer.csv link in gallery"), el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(102, "Then Move entity dialog should be visible", () => shouldBe(page, el("Move entity dialog"), "visible"));
      await session.step(103, "And choice input in Move entity dialog should have value \"Link\"", () => shouldHaveValue(page, el("choice input in Move entity dialog"), "Link"));
      await session.step(104, "When user clicks on YES button in Move entity dialog", () => clickOn(page, el("YES button in Move entity dialog")));
      await session.step(105, "Then Move entity dialog should be hidden", () => shouldBe(page, el("Move entity dialog"), "hidden"));
      await session.step(106, "When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree", () => doubleClickOn(page, el("Spaces---BDD-DnD-Src tree node inside browse tree")));
      await session.step(107, "Then the \"BDD-DnD-Src\" view should be current", () => viewIsCurrent(page, "BDD-DnD-Src"));
      await session.step(108, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
      await session.step(109, "And no errors should have been logged", () => noErrors(page));
      await session.step(110, "And no error or warning balloon should have been shown", () => noBalloons(page));
    });
    await run.scenario("The demo files kept their originals", async () => {
      await session.step(113, "When user clicks on \"Files > Demo\" tree node inside browse tree", () => clickOn(page, el("\"Files > Demo\" tree node inside browse tree")));
      await session.step(114, "Then the \"Demo\" view should be current", () => viewIsCurrent(page, "Demo"));
      await session.step(115, "And demog.csv link in gallery should be visible", () => shouldBe(page, el("demog.csv link in gallery"), "visible"));
      await session.step(116, "And TSLA.csv link in gallery should be visible", () => shouldBe(page, el("TSLA.csv link in gallery"), "visible"));
      await session.step(117, "And beer.csv link in gallery should be visible", () => shouldBe(page, el("beer.csv link in gallery"), "visible"));
    });
    run.finish();
  });
});
