/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-tree.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [views.browse]
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, collapse, expand, isExpanded, pressKey, shouldBe, shouldNotBe} from '@datagrok-libraries/bdd/bindings/common/steps';
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
    await session.step(36, "And Files---My-files tree node inside browse tree should be visible", () => shouldBe(page, el("Files---My-files tree node inside browse tree"), "visible"));
    await session.step(37, "When user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(38, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(39, "And Files---App-Data tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "hidden"));
    await session.step(40, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The arrow keys walk the tree and open a node", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(25, "Given user is logged in", () => loggedIn(page));
    await session.step(26, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(47, "Given user collapses Files tree node inside browse tree", () => collapse(page, el("Files tree node inside browse tree")));
    await session.step(48, "When user clicks on Files tree node inside browse tree", () => clickOn(page, el("Files tree node inside browse tree")));
    await session.step(49, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(50, "When user presses ArrowDown", () => pressKey(page, "ArrowDown"));
    await session.step(51, "Then Dashboards tree node inside browse tree should be selected", () => shouldBe(page, el("Dashboards tree node inside browse tree"), "selected"));
    await session.step(52, "And Files tree node inside browse tree should not be selected", () => shouldNotBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(53, "When user presses ArrowUp", () => pressKey(page, "ArrowUp"));
    await session.step(54, "Then Files tree node inside browse tree should be selected", () => shouldBe(page, el("Files tree node inside browse tree"), "selected"));
    await session.step(57, "When user presses ArrowRight", () => pressKey(page, "ArrowRight"));
    await session.step(58, "Then Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(59, "And Files---My-files tree node inside browse tree should be selected", () => shouldBe(page, el("Files---My-files tree node inside browse tree"), "selected"));
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
});
