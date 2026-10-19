/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-tree.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
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
import {clickOn, collapse, expand, isExpanded, pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browseTreeLoaded} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Working with the nodes of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-tree.feature", import.meta.url);
  test("A node opens and closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(34, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(36, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(37, "Then Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(38, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(39, "When user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(40, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(41, "And Files---App-Data tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "hidden"));
    await session.step(42, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
    await session.step(44, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The arrow keys walk the tree and open a node", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(49, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(50, "When user clicks on Files tree node inside browse tree", () => clickOn(page, el("Files tree node inside browse tree")));
    await session.step(51, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(52, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(53, "Then Dashboards tree node inside browse tree should be selected", () => shouldBe(page, el("Dashboards tree node inside browse tree"), "selected"));
    await session.step(54, "And Files tree node inside browse tree should not be selected", () => shouldNotBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(55, "When user presses ArrowUp", () => pressKey(page, "ArrowUp"));
    await session.step(56, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(60, "When user presses ArrowRight", () => pressKey(page, "ArrowRight"));
    await session.step(61, "Then Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(62, "And Files tree node inside browse tree should not be selected", () => shouldNotBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(63, "When user presses ArrowLeft", () => pressKey(page, "ArrowLeft"));
    await session.step(64, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(65, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(66, "When user presses ArrowLeft", () => pressKey(page, "ArrowLeft"));
    await session.step(67, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(68, "And no errors should have been logged", () => noErrors(page));
    await session.step(69, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The tree keeps what it had open across a visit to another view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(72, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(73, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(74, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(75, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(78, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(79, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(80, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(81, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(82, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(83, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(84, "And no errors should have been logged", () => noErrors(page));
    await session.step(85, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Hiding and showing the panel leaves a closed child closed", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(88, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(89, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(90, "And Files---Demo---chem tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "visible"));
    await session.step(91, "When user collapses Files---Demo tree node inside browse tree", () => collapse(page, el("Files---Demo tree node inside browse tree")));
    await session.step(92, "Then Files---Demo---chem tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "hidden"));
    await session.step(93, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(94, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(95, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(96, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(97, "And the browse tree should have finished loading", () => browseTreeLoaded(page));
    await session.step(98, "And Files tree node inside browse tree should be expanded", () => shouldBe(page, el("Files tree node inside browse tree"), "expanded"));
    await session.step(99, "And Files---Demo tree node inside browse tree should be collapsed", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "collapsed"));
    await session.step(100, "And Files---Demo---chem tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "hidden"));
    await session.step(101, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(102, "And no errors should have been logged", () => noErrors(page));
    await session.step(103, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("An arrow key opens the file it moves to", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(106, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(107, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(108, "When user clicks on Files---Demo---demog-1000.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(109, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(110, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(111, "Then Files---Demo---demog.csv tree node inside browse tree should be selected", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "selected"));
    await session.step(112, "And the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(113, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(114, "And no errors should have been logged", () => noErrors(page));
  });
  test("Enter keeps the previewed file open through the next single click", {tag: ["@browse", "@realizes:views.browse", "@known-failure:125"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(28, "Given user is logged in", () => loggedIn(page));
    await session.step(29, "And the browse panel is open", () => browsePanelOpen(page));
    await knownFailure(async () => {
      await session.step(118, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
      await session.step(119, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
      await session.step(120, "When user clicks on Files---Demo---demog-1000.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
      await session.step(121, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
      await session.step(122, "When user presses Enter", () => pressKey(page, "Enter"));
      await session.step(123, "And user clicks on Dashboards tree node inside browse tree", () => clickOn(page, el("Dashboards tree node inside browse tree")));
      await session.step(124, "Then Projects view should be visible", () => shouldBe(page, el("Projects view"), "visible"));
      await session.step(125, "And demog-1000 view should be present", () => shouldBe(page, el("demog-1000 view"), "present"));
    }, 125);
  });
});
