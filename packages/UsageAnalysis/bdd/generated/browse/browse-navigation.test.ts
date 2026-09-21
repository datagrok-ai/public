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
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(31, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]), [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]);
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Browse tab toggles the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(43, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(44, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(45, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(46, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(47, "And Apps tree node inside browse tree should be visible", () => shouldBe(page, el("Apps tree node inside browse tree"), "visible"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
    await session.step(49, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Home icon returns to the Home view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(54, "Given user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(55, "And the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(56, "When user clicks on \"Home\" icon inside browse toolbar", () => clickOn(page, el("\"Home\" icon inside browse toolbar")));
    await session.step(57, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(58, "And no errors should have been logged", () => noErrors(page));
    await session.step(59, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse tree closes an open section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(62, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(63, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(64, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(65, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(66, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(67, "And no errors should have been logged", () => noErrors(page));
    await session.step(68, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Find path reveals the node of the object that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(71, "Given Apps tree node inside browse tree is expanded", () => isExpanded(page, el("Apps tree node inside browse tree")));
    await session.step(72, "When user clicks on Tutorials tree node inside browse tree", () => clickOn(page, el("Tutorials tree node inside browse tree")));
    await session.step(73, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(74, "And user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(75, "Then Tutorials tree node inside browse tree should be hidden", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "hidden"));
    await session.step(76, "When user clicks on \"Find path\" icon inside browse toolbar", () => clickOn(page, el("\"Find path\" icon inside browse toolbar")));
    await session.step(77, "Then Tutorials tree node inside browse tree should be visible", () => shouldBe(page, el("Tutorials tree node inside browse tree"), "visible"));
    await session.step(78, "And no errors should have been logged", () => noErrors(page));
    await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open local file imports a CSV into a table view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(82, "When user uploads \"fixtures/browse-import.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/browse-import.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(83, "Then the \"browse-import\" view should be current", () => viewIsCurrent(page, "browse-import"));
    await session.step(84, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(85, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(86, "And the table should have a column \"population\"", () => hasColumn(page, "population"));
    await session.step(87, "And no errors should have been logged", () => noErrors(page));
    await session.step(88, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open text opens the text import view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(91, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(92, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(93, "And no errors should have been logged", () => noErrors(page));
    await session.step(94, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
