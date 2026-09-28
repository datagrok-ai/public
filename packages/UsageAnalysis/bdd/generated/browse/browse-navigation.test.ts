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
import {clickOn, collapse, expand, followingShouldBe, isExpanded, shouldBe, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {refreshBrowse} from '@datagrok-libraries/bdd/bindings/platform/browse';
import {columnCount, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, noScriptOnServer, scriptOnServer, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The Browse panel and the icons of its toolbar", () => {
  const session = feature(test, "features/browse/browse-navigation.feature", import.meta.url);
  test("The tree opens with every top-level section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(36, "Then the following elements should be visible:", () => followingShouldBe(page, "visible", [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]), [["My stuff tree node inside browse tree"],["Spaces tree node inside browse tree"],["Apps tree node inside browse tree"],["Files tree node inside browse tree"],["Dashboards tree node inside browse tree"],["Databases tree node inside browse tree"],["Platform tree node inside browse tree"]]);
    await session.step(44, "And no errors should have been logged", () => noErrors(page));
    await session.step(45, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Browse tab toggles the panel", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(48, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(49, "Then the browse tree should be hidden", () => shouldBe(page, el("the browse tree"), "hidden"));
    await session.step(50, "When user clicks on browse tab", () => clickOn(page, el("browse tab")));
    await session.step(51, "Then the browse tree should be visible", () => shouldBe(page, el("the browse tree"), "visible"));
    await session.step(52, "And Apps tree node inside browse tree should be visible", () => shouldBe(page, el("Apps tree node inside browse tree"), "visible"));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
    await session.step(54, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The Home icon returns to the Home view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(59, "Given user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(60, "And the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(61, "When user clicks on \"Home\" icon inside browse toolbar", () => clickOn(page, el("\"Home\" icon inside browse toolbar")));
    await session.step(62, "Then the \"Home\" view should be current", () => viewIsCurrent(page, "Home"));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
    await session.step(64, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Collapse tree closes an open section", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(67, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(68, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(69, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(70, "Then Files---Demo tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "hidden"));
    await session.step(71, "And Files tree node inside browse tree should be visible", () => shouldBe(page, el("Files tree node inside browse tree"), "visible"));
    await session.step(72, "And no errors should have been logged", () => noErrors(page));
    await session.step(73, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Find path reveals and selects the node of the object that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(76, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(77, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(78, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(79, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(80, "When user clicks on \"Collapse tree\" icon inside browse toolbar", () => clickOn(page, el("\"Collapse tree\" icon inside browse toolbar")));
    await session.step(81, "Then Files---Demo---demog.csv tree node inside browse tree should be hidden", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "hidden"));
    await session.step(82, "When user clicks on \"Find path\" icon inside browse toolbar", () => clickOn(page, el("\"Find path\" icon inside browse toolbar")));
    await session.step(83, "Then Files---Demo---demog.csv tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "visible"));
    await session.step(84, "And Files---Demo---demog.csv tree node inside browse tree should be selected", () => shouldBe(page, el("Files---Demo---demog.csv tree node inside browse tree"), "selected"));
    await session.step(85, "And no errors should have been logged", () => noErrors(page));
    await session.step(86, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open local file imports a CSV into a table view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(89, "When user uploads \"fixtures/browse-import.csv\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/browse-import.csv", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(90, "Then the \"browse-import\" view should be current", () => viewIsCurrent(page, "browse-import"));
    await session.step(91, "And the table should have 5 rows", () => rowCount(page, 5));
    await session.step(92, "And the table should have 3 columns", () => columnCount(page, 3));
    await session.step(93, "And the table should have a column \"population\"", () => hasColumn(page, "population"));
    await session.step(94, "And no errors should have been logged", () => noErrors(page));
    await session.step(95, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Open text opens the text import view", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(98, "When user clicks on \"Open text\" icon inside browse toolbar", () => clickOn(page, el("\"Open text\" icon inside browse toolbar")));
    await session.step(99, "Then the \"Import text\" view should be current", () => viewIsCurrent(page, "Import text"));
    await session.step(100, "And code editor should be visible", () => shouldBe(page, el("code editor"), "visible"));
    await session.step(101, "And no errors should have been logged", () => noErrors(page));
    await session.step(102, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Refresh brings in a script saved on the server meanwhile", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(105, "Given no script named \"BDD-Browse-Script-{run}\" is on the server", () => noScriptOnServer(page, session.text("BDD-Browse-Script-{run}")));
    await session.step(106, "And My stuff tree node inside browse tree is expanded", () => isExpanded(page, el("My stuff tree node inside browse tree")));
    await session.step(107, "And a script \"BDD-Browse-Script-{run}\" is on the server:", () => scriptOnServer(page, session.text("BDD-Browse-Script-{run}"), "//language: javascript\nlet x = 1;"));
    await session.step(112, "Then My-stuff---Scripts---BDD-Browse-Script-{run} tree node inside browse tree should be absent", () => shouldBe(page, el(session.text("My-stuff---Scripts---BDD-Browse-Script-{run} tree node inside browse tree")), "absent"));
    await session.step(113, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(114, "And user expands My-stuff---Scripts tree node inside browse tree", () => expand(page, el("My-stuff---Scripts tree node inside browse tree")));
    await session.step(115, "Then My-stuff---Scripts---BDD-Browse-Script-{run} tree node inside browse tree should be visible", () => shouldBe(page, el(session.text("My-stuff---Scripts---BDD-Browse-Script-{run} tree node inside browse tree")), "visible"));
    await session.step(116, "And no errors should have been logged", () => noErrors(page));
    await session.step(117, "And no error or warning balloon should have been shown", () => noBalloons(page));
    await session.step(119, "When user collapses My stuff tree node inside browse tree", () => collapse(page, el("My stuff tree node inside browse tree")));
  });
  test("Refresh keeps an open section open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(124, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(125, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(126, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(127, "Then Files tree node inside browse tree should be expanded", () => shouldBe(page, el("Files tree node inside browse tree"), "expanded"));
    await session.step(128, "And Files---Demo tree node inside browse tree should be visible", () => shouldBe(page, el("Files---Demo tree node inside browse tree"), "visible"));
    await session.step(129, "And Files---App-Data tree node inside browse tree should be visible", () => shouldBe(page, el("Files---App-Data tree node inside browse tree"), "visible"));
    await session.step(130, "And no errors should have been logged", () => noErrors(page));
    await session.step(131, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("Refresh keeps the preview that is open", {tag: ["@browse", "@realizes:views.browse"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(32, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(134, "Given Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(135, "And Files---Demo tree node inside browse tree is expanded", () => isExpanded(page, el("Files---Demo tree node inside browse tree")));
    await session.step(136, "When user clicks on Files---Demo---demog.csv tree node inside browse tree", () => clickOn(page, el("Files---Demo---demog.csv tree node inside browse tree")));
    await session.step(137, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(138, "When user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(139, "Then the \"demog\" view should be current", () => viewIsCurrent(page, "demog"));
    await session.step(140, "And demog view should be visible", () => shouldBe(page, el("demog view"), "visible"));
    await session.step(141, "And grid should show 5850 rows", () => showsRows(page, el("grid"), 5850));
    await session.step(142, "And no errors should have been logged", () => noErrors(page));
    await session.step(143, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
