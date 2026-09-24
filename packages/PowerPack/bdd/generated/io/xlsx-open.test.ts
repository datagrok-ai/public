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
import '../../bindings/formula-lines.js';
import '../../bindings/home.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {cellOfTable, columnTypeOfTable, dropFile, fixtureInHome, menuFileChooser, noTableViewOpen, refreshBrowseTree, tableViewsOpen} from '../../bindings/io.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, doubleClickOn, isExpanded, uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {tableColumns, tableRows} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelOpen, viewIsCurrent} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noBalloons, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("An Excel workbook opens from every entry path", () => {
  const session = feature(test, "features/io/xlsx-open.feature", import.meta.url);
  test("The workbook opens from My files in the Browse tree", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(30, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(31, "And the \"fixtures/xlsx-open-test.xlsx\" file of the project is in the home folder as \"xlsx-open-{time}.xlsx\"", () => fixtureInHome(page, "fixtures/xlsx-open-test.xlsx", session.text("xlsx-open-{time}.xlsx")));
    await session.step(32, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(33, "And user refreshes the browse tree", () => refreshBrowseTree(page));
    await session.step(34, "And Files tree node inside browse tree is expanded", () => isExpanded(page, el("Files tree node inside browse tree")));
    await session.step(35, "And Files---My-files tree node inside browse tree is expanded", () => isExpanded(page, el("Files---My-files tree node inside browse tree")));
    await session.step(36, "When user double-clicks on Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree", () => doubleClickOn(page, el(session.text("Files---My-files---xlsx-open-{time}.xlsx tree node inside browse tree"))));
    await session.step(37, "Then the table views \"Customers, Orders, Products\" should be open", () => tableViewsOpen(page, "Customers, Orders, Products"));
    await session.step(38, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(39, "And table \"Customers\" should have 5 rows", () => tableRows(page, "Customers", 5));
    await session.step(40, "And table \"Customers\" should have columns \"CustomerID, Name, Country, Since\"", () => tableColumns(page, "Customers", "CustomerID, Name, Country, Since"));
    await session.step(41, "And the value of \"Name\" column in row 2 of table \"Customers\" should be \"Bottom Dollar\"", () => cellOfTable(page, "Name", 2, "Customers", "Bottom Dollar"));
    await session.step(42, "And \"Since\" column of table \"Customers\" should have type \"datetime\"", () => columnTypeOfTable(page, "Since", "Customers", "datetime"));
    await session.step(43, "And table \"Orders\" should have 6 rows", () => tableRows(page, "Orders", 6));
    await session.step(44, "And table \"Orders\" should have columns \"OrderID, CustomerID, Amount, Shipped\"", () => tableColumns(page, "Orders", "OrderID, CustomerID, Amount, Shipped"));
    await session.step(45, "And the value of \"Amount\" column in row 3 of table \"Orders\" should be \"1200.00\"", () => cellOfTable(page, "Amount", 3, "Orders", "1200.00"));
    await session.step(46, "And table \"Products\" should have 4 rows", () => tableRows(page, "Products", 4));
    await session.step(47, "And table \"Products\" should have columns \"ProductID, Product, Price, InStock\"", () => tableColumns(page, "Products", "ProductID, Product, Price, InStock"));
    await session.step(48, "And the value of \"Product\" column in row 1 of table \"Products\" should be \"Chai\"", () => cellOfTable(page, "Product", 1, "Products", "Chai"));
    await session.step(49, "And grid should show 4 rows", () => showsRows(page, el("grid"), 4));
    await session.step(50, "And no errors should have been logged", () => noErrors(page));
    await session.step(51, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens when dropped onto the window", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(54, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(55, "When user drops the \"fixtures/xlsx-open-test.xlsx\" file of the project onto home widgets panel", () => dropFile(page, "fixtures/xlsx-open-test.xlsx", el("home widgets panel")));
    await session.step(56, "Then the table views \"Customers, Orders, Products\" should be open", () => tableViewsOpen(page, "Customers, Orders, Products"));
    await session.step(57, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(58, "And table \"Customers\" should have 5 rows", () => tableRows(page, "Customers", 5));
    await session.step(59, "And the value of \"Country\" column in row 3 of table \"Customers\" should be \"Switzerland\"", () => cellOfTable(page, "Country", 3, "Customers", "Switzerland"));
    await session.step(60, "And table \"Orders\" should have 6 rows", () => tableRows(page, "Orders", 6));
    await session.step(61, "And table \"Products\" should have 4 rows", () => tableRows(page, "Products", 4));
    await session.step(62, "And the value of \"Price\" column in row 4 of table \"Products\" should be \"23.25\"", () => cellOfTable(page, "Price", 4, "Products", "23.25"));
    await session.step(63, "And grid should show 4 rows", () => showsRows(page, el("grid"), 4));
    await session.step(64, "And no errors should have been logged", () => noErrors(page));
    await session.step(65, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens through the Open local file icon and the file chooser", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(68, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(69, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(70, "When user uploads \"fixtures/xlsx-open-test.xlsx\" through \"Open local file\" icon in browse toolbar", () => uploadThrough(page, "fixtures/xlsx-open-test.xlsx", el("\"Open local file\" icon in browse toolbar")));
    await session.step(71, "Then the table views \"Customers, Orders, Products\" should be open", () => tableViewsOpen(page, "Customers, Orders, Products"));
    await session.step(72, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(73, "And table \"Customers\" should have 5 rows", () => tableRows(page, "Customers", 5));
    await session.step(74, "And table \"Orders\" should have 6 rows", () => tableRows(page, "Orders", 6));
    await session.step(75, "And the value of \"Shipped\" column in row 2 of table \"Orders\" should be \"no\"", () => cellOfTable(page, "Shipped", 2, "Orders", "no"));
    await session.step(76, "And table \"Products\" should have 4 rows", () => tableRows(page, "Products", 4));
    await session.step(77, "And grid should show 4 rows", () => showsRows(page, el("grid"), 4));
    await session.step(78, "And no errors should have been logged", () => noErrors(page));
    await session.step(79, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
  test("The workbook opens through File > Open > File... and the file chooser", {tag: ["@realizes:powerpack.import.xlsx"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(27, "Given user is logged in", () => loggedIn(page));
    await session.step(82, "Given no table view is open", () => noTableViewOpen(page));
    await session.step(83, "When user clicks on Menu tab", () => clickOn(page, el("Menu tab")));
    await session.step(84, "And user picks \"File > Open > File...\" from the open menu and chooses the \"fixtures/xlsx-open-test.xlsx\" file of the project", () => menuFileChooser(page, "File > Open > File...", "fixtures/xlsx-open-test.xlsx"));
    await session.step(85, "Then the table views \"Customers, Orders, Products\" should be open", () => tableViewsOpen(page, "Customers, Orders, Products"));
    await session.step(86, "And the \"Products\" view should be current", () => viewIsCurrent(page, "Products"));
    await session.step(87, "And table \"Customers\" should have 5 rows", () => tableRows(page, "Customers", 5));
    await session.step(88, "And table \"Orders\" should have 6 rows", () => tableRows(page, "Orders", 6));
    await session.step(89, "And table \"Products\" should have 4 rows", () => tableRows(page, "Products", 4));
    await session.step(90, "And the value of \"Product\" column in row 2 of table \"Products\" should be \"Chang\"", () => cellOfTable(page, "Product", 2, "Products", "Chang"));
    await session.step(91, "And no errors should have been logged", () => noErrors(page));
    await session.step(92, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
