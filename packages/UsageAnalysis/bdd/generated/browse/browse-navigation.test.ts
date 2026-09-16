/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-navigation.feature
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
import {clickOn, followingShouldBe, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Browse panel and the icons of its toolbar", () => {
  const session = feature(test, "features/browse/browse-navigation.feature", import.meta.url);
  test("The tree opens with every top-level section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(30, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
    await session.step(39, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Browse tab toggles the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(42, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(43, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(44, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(45, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(46, "And Apps tree node inside browse tree should be visible", () => shouldBe(page, el("Apps tree node inside browse tree"), "visible"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
    await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Home icon returns to the Home view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(53, "Given user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(54, "And the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(55, "When user clicks on \"Home\" icon inside browse toolbar", () => clickOn(page, el("\"Home\" icon inside browse toolbar")));
    await session.step(56, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(57, "And no errors should have been logged", () => noErrors(page));
    await session.step(58, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse tree closes an open section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(61, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(62, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(63, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(64, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(65, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(66, "And no errors should have been logged", () => noErrors(page));
    await session.step(67, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Find path reveals the node of the object that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(70, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(71, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(72, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "And user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(74, "Then Tutorials tree node inside browse tree should be hidden", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "hidden"));
    await session.step(75, "When user clicks on \"Find path\" icon inside browse toolbar", () => clickOn(page, el("\"Find path\" icon inside browse toolbar")));
    await session.step(76, "Then Tutorials tree node inside browse tree should be visible", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "visible"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
    await session.step(78, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open local file imports a CSV into a table view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(81, "When user uploads \"fixtures/browse-import.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/browse-import.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(82, "Then the \"browse-import\" view should be current", () => viewIsCurrent(page, "browse-import"));
    await session.step(83, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(84, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(85, "And the table should have a column \"population\"", () => hasColumn(page, "population"));
    await session.step(86, "And no errors should have been logged", () => noErrors(page));
    await session.step(87, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open text opens the text import view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(26, "Given user is logged in", () => loggedIn(page));
    await session.step(27, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(90, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(91, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(92, "And no errors should have been logged", () => noErrors(page));
    await session.step(93, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
