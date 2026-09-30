/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/browse/browse-navigation.feature
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
import {clickOn, collapse, followingShouldBe, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, expandMyStuffBucket, listedInMyStuffBucket, noScriptOnServer, notListedInMyStuffBucket, refreshBrowse, scriptOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Browse panel and the icons of its toolbar", () => {
  const session = feature(test, "features/browse/browse-navigation.feature", import.meta.url);
  test("The tree opens with every top-level section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]), [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]);
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Browse tab toggles the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(45, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(46, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(47, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(48, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(49, "And Apps tree node inside browse tree should be visible", () => shouldBe(page, el("Apps tree node inside browse tree"), "visible"));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Home icon returns to the Home view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(56, "Given user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(57, "And the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(58, "When user clicks on \"Home\" icon inside browse toolbar", () => clickOn(page, el("\"Home\" icon inside browse toolbar")));
    await session.step(59, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(60, "And no errors should have been logged", () => noErrors(page));
    await session.step(61, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse tree closes an open section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(64, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(65, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(66, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(67, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(68, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
    await session.step(70, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Find path reveals and selects the node of the object that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(73, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(74, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(75, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(76, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(77, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(78, "Then Files---Demo---demog.csv tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "hidden"));
    await session.step(79, "When user clicks on \"Find path\" icon inside browse toolbar", () => clickOn(page, el("\"Find path\" icon inside browse toolbar")));
    await session.step(80, "Then Files---Demo---demog.csv tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "visible"));
    await session.step(81, "And Files---Demo---demog.csv tree node inside browse tree should be selected", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "selected"));
    await session.step(82, "And no errors should have been logged", () => noErrors(page));
    await session.step(83, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open local file imports a CSV into a table view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(86, "When user uploads \"fixtures/browse-import.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/browse-import.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(87, "Then the \"browse-import\" view should be current", () => viewIsCurrent(page, "browse-import"));
    await session.step(88, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(89, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(90, "And the table should have a column \"population\"", () => hasColumn(page, "population"));
    await session.step(91, "And no errors should have been logged", () => noErrors(page));
    await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open text opens the text import view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(95, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(96, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(97, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
    await session.step(98, "And no errors should have been logged", () => noErrors(page));
    await session.step(99, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Refresh brings in a script saved on the server meanwhile", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(102, "Given no script named \"BDD-Browse-Script-{run}\" is on the server", () => noScriptOnServer(page, session.text("BDD-Browse-Script-{run}")));
    await session.step(103, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(104, "And a script \"BDD-Browse-Script-{run}\" is on the server:", () => scriptOnServer(page, session.text("BDD-Browse-Script-{run}"), "//language: javascript\nlet x = 1;"));
    await session.step(109, "Then \"BDD-Browse-Script-{run}\" should not be listed in the \"Scripts\" bucket of My stuff", () => notListedInMyStuffBucket(page, session.text("BDD-Browse-Script-{run}"), "Scripts"));
    await session.step(110, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(111, "And user expands the \"Scripts\" bucket of My stuff", () => expandMyStuffBucket(page, "Scripts"));
    await session.step(112, "Then \"BDD-Browse-Script-{run}\" should be listed in the \"Scripts\" bucket of My stuff", () => listedInMyStuffBucket(page, session.text("BDD-Browse-Script-{run}"), "Scripts"));
    await session.step(113, "And no errors should have been logged", () => noErrors(page));
    await session.step(114, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(116, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
  });
  test("Refresh keeps an open section open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(121, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(122, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(123, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(124, "Then Files tree node inside browse tree should be expanded", () => shouldBe(page, el("Files tree node inside browse tree"), "expanded"));
    await session.step(125, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(126, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(127, "And no errors should have been logged", () => noErrors(page));
    await session.step(128, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Refresh keeps the preview that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(131, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(132, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(133, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(134, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(135, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(136, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(137, "And demog view should be visible", () => shouldBe(page, el("demog view"), "visible"));
    await session.step(138, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(139, "And no errors should have been logged", () => noErrors(page));
    await session.step(140, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The panel's close icon hides it and the Browse tab brings it back", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(29, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(143, "When user clicks on browse panel close icon", () => clickOn(page, el("browse panel close icon")));
    await session.step(144, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(145, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(146, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(147, "And Apps tree node inside browse tree should be visible", () => shouldBe(page, el("Apps tree node inside browse tree"), "visible"));
    await session.step(148, "And no errors should have been logged", () => noErrors(page));
    await session.step(149, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
