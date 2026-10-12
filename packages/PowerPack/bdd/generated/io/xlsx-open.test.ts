/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/io/xlsx-open.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [powerpack.import.xlsx]
--- */
import {test} from '@playwright/test';
import '../../bindings/add-new-column.js';
import '../../bindings/enrichment.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelOpen, refreshBrowse, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {dropFile, fixtureInUsersFiles, menuFileChooser, noTableViewOpen, openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {noBalloons, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("An Excel workbook opens from every entry path", () => {
  const session = feature(test, "features/io/xlsx-open.feature", import.meta.url);
  test("The workbook opens from My files in the Browse tree", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(31, "And the \"fixtures/xlsx-open-test.xlsx\" file of the project is in the user's files as \"xlsx-open-{time}.xlsx\"", () => fixtureInUsersFiles(page, "fixtures/xlsx-open-test.xlsx", session.text("xlsx-open-{time}.xlsx")));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "And user refreshes the browse tree", () => refreshBrowse(page));
    await session.step(34, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree", () => doubleClickOn(page, el(session.text("Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree"))));
    await session.step(37, "Then the open table views should be exactly \"Customers, Orders, Products\"", () => openTableViewsExactly(page, "Customers, Orders, Products"));
    await session.step(38, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(39, "And no errors should have been logged", () => noErrors(page));
    await session.step(40, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens when dropped onto the window", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(43, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(44, "When user drops the \"fixtures/xlsx-open-test.xlsx\" file of the project onto home widgets panel", () => dropFile(page, "fixtures/xlsx-open-test.xlsx", el("home widgets panel")));
    await session.step(45, "Then the open table views should be exactly \"Customers, Orders, Products\"", () => openTableViewsExactly(page, "Customers, Orders, Products"));
    await session.step(46, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
    await session.step(48, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens through the Open local file icon and the file chooser", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(52, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(53, "When user uploads \"fixtures/xlsx-open-test.xlsx\" through \"Open local file\" icon in browse toolbar", () => uploadThrough(page, "fixtures/xlsx-open-test.xlsx", el("\"Open local file\" icon in browse toolbar")));
    await session.step(54, "Then the open table views should be exactly \"Customers, Orders, Products\"", () => openTableViewsExactly(page, "Customers, Orders, Products"));
    await session.step(55, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(56, "And no errors should have been logged", () => noErrors(page));
    await session.step(57, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens through File > Open > File... and the file chooser", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(60, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(61, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
    await session.step(62, "And user picks \"File > Open > File...\" from the open menu and chooses the \"fixtures/xlsx-open-test.xlsx\" file of the project", () => menuFileChooser(page, "File > Open > File...", "fixtures/xlsx-open-test.xlsx"));
    await session.step(63, "Then the open table views should be exactly \"Customers, Orders, Products\"", () => openTableViewsExactly(page, "Customers, Orders, Products"));
    await session.step(64, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(65, "And no errors should have been logged", () => noErrors(page));
    await session.step(66, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
