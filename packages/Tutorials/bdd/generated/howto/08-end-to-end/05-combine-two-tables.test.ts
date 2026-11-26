/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/08-end-to-end/05-combine-two-tables.feature
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
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {pickAggregation} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, selectIn, shouldBe, shouldContainText, shouldHaveValue, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, clickViewTab, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {hasArea, hoverArea, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector, toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · End to end", () => {
  const session = feature(test, "features/howto/08-end-to-end/05-combine-two-tables.feature", import.meta.url);
  test("Combine two tables to answer a question", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens northwind-order-details dataset", () => openDataset(page, ds("northwind-order-details")));
    await session.step(18, "And user opens northwind-products dataset", () => openDataset(page, ds("northwind-products")));
    await session.step(19, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(20, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(21, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(22, "When the video says \"The task: which products sell the most\"", () => videoSays(page, "The task: which products sell the most"));
    await session.step(24, "And user clicks on the tab of the \"order-details\" view", () => clickViewTab(page, "order-details"), undefined, "The order lines have the quantities, but name a product only by its ID");
    await session.step(25, "And user hovers over the \"cell 1 of ProductID\" area of grid", () => hoverArea(page, "cell 1 of ProductID", el("grid")));
    await session.step(27, "And user clicks on the tab of the \"products\" view", () => clickViewTab(page, "products"), undefined, "The products table has the names");
    await session.step(28, "And user hovers over the \"cell 1 of ProductName\" area of grid", () => hoverArea(page, "cell 1 of ProductName", el("grid")));
    await session.step(30, "And user clicks on the tab of the \"order-details\" view", () => clickViewTab(page, "order-details"), undefined, "Join Tables brings the names over, matched on ProductID");
    await session.step(31, "And user picks \"Data > Join Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Join Tables..."));
    await session.step(32, "Then Join Tables dialog should be visible", () => shouldBe(page, el("Join Tables dialog"), "visible"));
    await session.step(33, "And \"Tables\" input in Join Tables dialog should have value \"order-details\"", () => shouldHaveValue(page, el("\"Tables\" input in Join Tables dialog"), "order-details"));
    await session.step(34, "When user selects \"ProductID\" in Key Columns input in Join Tables dialog", () => selectIn(page, "ProductID", el("Key Columns input in Join Tables dialog")));
    await session.step(35, "Then Join Tables dialog should contain text \"All key rows are matched\"", () => shouldContainText(page, el("Join Tables dialog"), "All key rows are matched"));
    await session.step(37, "When user clicks on 2nd columns input in Join Tables dialog", () => clickOn(page, el("2nd columns input in Join Tables dialog")), undefined, "From products, take just the name");
    await session.step(38, "And user clicks on \"None\" link in \"Select columns...\" dialog", () => clickOn(page, el("\"None\" link in \"Select columns...\" dialog")));
    await session.step(39, "And user toggles the \"ProductName\" column in the column list of \"Select columns...\" dialog", () => toggleInColumnList(page, "ProductName", el("\"Select columns...\" dialog")));
    await session.step(40, "And user clicks on OK button in \"Select columns...\" dialog", () => clickOn(page, el("OK button in \"Select columns...\" dialog")));
    await session.step(42, "And user selects \"left\" in \"Join Type\" input in Join Tables dialog", () => selectIn(page, "left", el("\"Join Type\" input in Join Tables dialog")), undefined, "Left join keeps every order line");
    await session.step(43, "And user clicks on OK button in Join Tables dialog", () => clickOn(page, el("OK button in Join Tables dialog")));
    await session.step(44, "Then Join Tables dialog should be hidden", () => shouldBe(page, el("Join Tables dialog"), "hidden"));
    await session.step(45, "And the open table views should be exactly \"order-details, products, result\"", () => openTableViewsExactly(page, "order-details, products, result"));
    await session.step(46, "And the table should have 2155 rows", () => rowCount(page, 2155));
    await session.step(47, "And the table should have the columns \"OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName\"", () => columnsExactly(page, "OrderID, ProductID, UnitPrice, Quantity, Discount, ProductName"));
    await session.step(48, "When the video says \"The result: every order line with its product name\"", () => videoSays(page, "The result: every order line with its product name"));
    await session.step(49, "And user hovers over the \"cell 1 of ProductName\" area of grid", () => hoverArea(page, "cell 1 of ProductName", el("grid")));
    await session.step(51, "And user clicks on \"bar chart\" icon in toolbox", () => clickOn(page, el("\"bar chart\" icon in toolbox")), undefined, "Now the answer: a bar chart of quantity per product");
    await session.step(52, "Then bar chart viewer should be visible", () => shouldBe(page, el("bar chart viewer"), "visible"));
    await session.step(53, "When user picks \"ProductName\" in the \"split\" column selector of bar chart viewer", () => pickInColumnSelector(page, "ProductName", "split", el("bar chart viewer")));
    await session.step(54, "And user picks \"Quantity\" in the \"value\" column selector of bar chart viewer", () => pickInColumnSelector(page, "Quantity", "value", el("bar chart viewer")));
    await session.step(56, "And user picks the \"sum\" aggregation in the \"value\" column selector of bar chart viewer", () => pickAggregation(page, "sum", "value", el("bar chart viewer")), undefined, "Sum the quantities");
    await session.step(57, "Then bar chart viewer should have a \"bar Camembert Pierrot\" area", () => hasArea(page, el("bar chart viewer"), "bar Camembert Pierrot"));
    await session.step(58, "When the video says \"Camembert Pierrot sells the most: 1,577 units\"", () => videoSays(page, "Camembert Pierrot sells the most: 1,577 units"));
    await session.step(59, "And user hovers over the \"bar Camembert Pierrot\" area of bar chart viewer", () => hoverArea(page, "bar Camembert Pierrot", el("bar chart viewer")));
    await session.step(60, "Then no errors should have been logged", () => noErrors(page));
  });
});
