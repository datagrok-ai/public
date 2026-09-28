/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-tree.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/connections.js';
import '../../bindings/grid.js';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, isExpanded, pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browseTreeLoaded} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Working with the nodes of the Browse tree", () => {
  const session = feature(test, "features/browse/browse-tree.feature", import.meta.url);
  test("A node opens and closes on its own twistie", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(32, "And Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(33, "When user expands Files tree node inside browse tree", () => expand(page, el("Files tree node inside browse tree")));
    await session.step(34, "Then Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(35, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(36, "When user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(37, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(38, "And Files---App-Data tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "hidden"));
    await session.step(39, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
    await session.step(41, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The arrow keys walk the tree and open a node", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(46, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(47, "When user clicks on Files tree node inside browse tree", () => clickOn(page, el("Files tree node inside browse tree")));
    await session.step(48, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(49, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(50, "Then Dashboards tree node inside browse tree should be selected", () => shouldBe(page, el("Dashboards tree node inside browse tree"), "selected"));
    await session.step(51, "And Files tree node inside browse tree should not be selected", () => shouldNotBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(52, "When user presses ArrowUp", () => pressKey(page, "ArrowUp"));
    await session.step(53, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(57, "When user presses ArrowRight", () => pressKey(page, "ArrowRight"));
    await session.step(58, "Then Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(59, "And Files tree node inside browse tree should not be selected", () => shouldNotBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(60, "When user presses ArrowLeft", () => pressKey(page, "ArrowLeft"));
    await session.step(61, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(62, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(63, "When user presses ArrowLeft", () => pressKey(page, "ArrowLeft"));
    await session.step(64, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
    await session.step(66, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The tree keeps what it had open across a visit to another view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(69, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(70, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(71, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(72, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(75, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(76, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(77, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(78, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(79, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(80, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(81, "And no errors should have been logged", () => noErrors(page));
    await session.step(82, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Hiding and showing the panel leaves a closed child closed", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(85, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(86, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(87, "And Files---Demo---chem tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "visible"));
    await session.step(88, "When user collapses Files---Demo tree node inside browse tree", () => collapse(page, el("Files---Demo tree node inside browse tree")));
    await session.step(89, "Then Files---Demo---chem tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "hidden"));
    await session.step(90, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(91, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(92, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(93, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(94, "And the browse tree should have finished loading", () => browseTreeLoaded(page));
    await session.step(95, "And Files tree node inside browse tree should be expanded", () => shouldBe(page, el("Files tree node inside browse tree"), "expanded"));
    await session.step(96, "And Files---Demo tree node inside browse tree should be collapsed", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "collapsed"));
    await session.step(97, "And Files---Demo---chem tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---chem tree node inside browse tree"), "hidden"));
    await session.step(98, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(99, "And no errors should have been logged", () => noErrors(page));
    await session.step(100, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Enter on a selected item opens it", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(103, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(104, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(105, "When user clicks on Files---Demo---demog-1000.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog-1000.csv tree node inside browse tree")));
    await session.step(106, "Then the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(107, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(108, "Then Files---Demo---demog.csv tree node inside browse tree should be selected", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "selected"));
    await session.step(110, "And the \"demog-1000\" view should be current", () => viewIsCurrent(page, "demog-1000"));
    await session.step(111, "When user presses Enter", () => pressKey(page, "Enter"));
    await session.step(112, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
  });
});
