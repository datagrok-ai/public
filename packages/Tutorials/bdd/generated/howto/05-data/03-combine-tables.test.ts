/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/03-combine-tables.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '../../../bindings/howto-basics.js';
import '../../../bindings/howto-compute.js';
import '../../../bindings/howto-curves.js';
import '../../../bindings/howto-data-2.js';
import '../../../bindings/howto-data.js';
import '../../../bindings/howto-databases.js';
import '../../../bindings/howto-filters.js';
import '../../../bindings/howto-ml.js';
import '../../../bindings/howto-molecules.js';
import '../../../bindings/howto-more-viewers.js';
import '../../../bindings/howto-scripting.js';
import '../../../bindings/howto-sequences.js';
import '../../../bindings/howto-viewers.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnComplete, columnsExactly, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, clickViewTab, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {hoverArea, noErrors, wheelOverAreaTimes} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/03-combine-tables.feature", import.meta.url);
  test("Combine tables by joining or stacking them", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(20, "Given user is logged in", () => loggedIn(page));
    await session.step(21, "And simple mode is off", () => simpleModeOff(page));
    await session.step(22, "And user opens plate-96-1 dataset", () => openDataset(page, ds("plate-96-1")));
    await session.step(23, "And user opens plate-96-2 dataset", () => openDataset(page, ds("plate-96-2")));
    await session.step(24, "And user opens northwind-order-details dataset", () => openDataset(page, ds("northwind-order-details")));
    await session.step(25, "And user opens northwind-products dataset", () => openDataset(page, ds("northwind-products")));
    await session.step(26, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(27, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(28, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(29, "When the video says \"products: the name of each product ID\"", () => videoSays(page, "products: the name of each product ID"));
    await session.step(30, "And user clicks on the tab of the \"order-details\" view", () => clickViewTab(page, "order-details"));
    await session.step(31, "And the video says \"Order lines refer to products only by ID\"", () => videoSays(page, "Order lines refer to products only by ID"));
    await session.step(33, "And user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."), undefined, "Join Tables: add columns from a table that shares a key");
    await session.step(34, "Then Join Tables dialog should be visible", () => shouldBe(page, el("Join Tables dialog"), "visible"));
    await session.step(36, "When user selects \"order-details\" in \"Tables\" input in Join Tables dialog", () => selectIn(page, "order-details", el("\"Tables\" input in Join Tables dialog")), undefined, "With more tables open, pick the two to join");
    await session.step(37, "And user selects \"products\" in 2nd choice input in Join Tables dialog", () => selectIn(page, "products", el("2nd choice input in Join Tables dialog")));
    await session.step(38, "Then \"Tables\" input in Join Tables dialog should have value \"order-details\"", () => shouldHaveValue(page, el("\"Tables\" input in Join Tables dialog"), "order-details"));
    await session.step(39, "And 2nd choice input in Join Tables dialog should have value \"products\"", () => shouldHaveValue(page, el("2nd choice input in Join Tables dialog"), "products"));
    await session.step(41, "When user selects \"ProductID\" in Key Columns input in Join Tables dialog", () => selectIn(page, "ProductID", el("Key Columns input in Join Tables dialog")), undefined, "Check the key on both sides: ProductID");
    await session.step(42, "Then Join Tables dialog should contain text \"All key rows are matched\"", () => shouldContainText(page, el("Join Tables dialog"), "All key rows are matched"));
    await session.step(44, "When user clicks on 2nd columns input in Join Tables dialog", () => clickOn(page, el("2nd columns input in Join Tables dialog")), undefined, "Choose the columns to bring over from each table");
    await session.step(45, "And user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(46, "And user toggles the \"ProductName\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "ProductName", el("\"Select columns...\" dialog")));
    await session.step(47, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(49, "And user selects \"left\" in \"Join Type\" input in Join Tables dialog", () => selectIn(page, "left", el("\"Join Type\" input in Join Tables dialog")), undefined, "Inner keeps matched rows only; left keeps every order line");
    await session.step(50, "And user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(51, "Then Join Tables dialog should be hidden", () => shouldBe(page, el("Join Tables dialog"), "hidden"));
    await session.step(52, "And the open table views should be exactly \"plate-96-1, plate-96-2, order-details, products, result\"", () => openTableViewsExactly(page, "plate-96-1, plate-96-2, order-details, products, result"));
    await session.step(53, "And the table should have 2155 rows", () => rowCount(page, 2155));
    await session.step(54, "And the table should have the columns \"OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName\"", () => columnsExactly(page, "OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName"));
    await session.step(55, "And \"ProductName\" column should have no missing values", () => columnComplete(page, "ProductName"));
    await session.step(57, "And the value of \"ProductName\" column in row 1 should be \"Chai\"", () => valueInRow(page, "ProductName", 1, "Chai"));
    await session.step(58, "When the video says \"A new table: every order line with its product's name\" at the bottom", () => videoSaysAt(page, "A new table: every order line with its product's name", "bottom"));
    await session.step(59, "And user hovers over the \"cell 1 of ProductName\" area of grid", () => hoverArea(page, "cell 1 of ProductName", el("grid")));
    await session.step(60, "When user clicks on the tab of the \"plate-96-1\" view", () => clickViewTab(page, "plate-96-1"));
    await session.step(62, "And user picks \"Data > Append Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Append Tables..."), undefined, "Two plates with the same columns: Append Tables stacks them");
    await session.step(63, "Then Append Tables dialog should be visible", () => shouldBe(page, el("Append Tables dialog"), "visible"));
    await session.step(65, "When user checks \"plate-96-1\" input in Append Tables dialog", () => check(page, el("\"plate-96-1\" input in Append Tables dialog")), undefined, "Tick the tables to stack");
    await session.step(66, "And user checks \"plate-96-2\" input in Append Tables dialog", () => check(page, el("\"plate-96-2\" input in Append Tables dialog")));
    await session.step(67, "And user clicks on OK button in Append Tables dialog", () => clickOn(page, el("OK button in Append Tables dialog")));
    await session.step(68, "Then Append Tables dialog should be hidden", () => shouldBe(page, el("Append Tables dialog"), "hidden"));
    await session.step(69, "And the open table views should be exactly \"plate-96-1, plate-96-2, order-details, products, result, result (2)\"", () => openTableViewsExactly(page, "plate-96-1, plate-96-2, order-details, products, result, result (2)"));
    await session.step(70, "And the table should have 192 rows", () => rowCount(page, 192));
    await session.step(71, "And the table should have the columns \"barcode, row, col, pos, volume, concentration, activity, role\"", () => columnsExactly(page, "barcode, row, col, pos, volume, concentration, activity, role"));
    await session.step(72, "When user scrolls the mouse wheel down 8 times over the \"cell 15 of pos\" area of grid", () => wheelOverAreaTimes(page, "down", 8, "cell 15 of pos", el("grid")));
    await session.step(74, "And user hovers over the \"cell 97 of barcode\" area of grid", () => hoverArea(page, "cell 97 of barcode", el("grid")), undefined, "192 rows in one table; plate 2 starts at row 97");
    await session.step(75, "Then the value of \"barcode\" column in row 96 should be \"1737398\"", () => valueInRow(page, "barcode", 96, "1737398"));
    await session.step(76, "And the value of \"barcode\" column in row 97 should be \"5899277\"", () => valueInRow(page, "barcode", 97, "5899277"));
    await session.step(77, "And no errors should have been logged", () => noErrors(page));
  });
});
