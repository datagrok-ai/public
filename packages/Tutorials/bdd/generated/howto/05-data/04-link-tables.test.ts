/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/04-link-tables.feature
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
import {clickOn, selectIn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {pickFromTopMenu} from '@datagrok-libraries/bdd/bindings/platform/commands';
import {selectedRowCount, setLinkTables, tableFilterCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, clickViewTab, contextPanelClosed, dialogCloses, dockViewBeside, openDataset, simpleModeOff, toolboxPaneHidden, viewDockedBeside} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, clickAreaHolding, noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/04-link-tables.feature", import.meta.url);
  test("Link tables", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens northwind-order-details dataset", () => openDataset(page, ds("northwind-order-details")));
    await session.step(17, "And user opens northwind-products dataset", () => openDataset(page, ds("northwind-products")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user drags the tab of the \"order-details\" view to the right side of the \"products\" view", () => dockViewBeside(page, "order-details", "right", "products"));
    await session.step(23, "And user clicks on the tab of the \"products\" view", () => clickViewTab(page, "products"));
    await session.step(24, "And the \"order-details\" view should be docked right-of the \"products\" view", () => viewDockedBeside(page, "order-details", "right-of", "products"));
    await session.step(25, "When the video says \"products and its order lines, side by side\"", () => videoSays(page, "products and its order lines, side by side"));
    await session.step(27, "And user picks \"Data > Link Tables...\" from the top menu", () => pickFromTopMenu(page, "Data > Link Tables..."), undefined, "Link Tables: what you do in one table shows in the other");
    await session.step(28, "Then \"Link Tables\" dialog should be visible", () => shouldBe(page, el("\"Link Tables\" dialog"), "visible"));
    await session.step(30, "When user sets the tables of the Link Tables dialog to \"products\" and \"order-details\"", () => setLinkTables(page, "products", "order-details"), undefined, "products drives order-details, matched on ProductID");
    await session.step(31, "And user selects \"ProductID\" in 2nd column input in \"Link Tables\" dialog", () => selectIn(page, "ProductID", el("2nd column input in \"Link Tables\" dialog")));
    await session.step(32, "Then Key Columns input in \"Link Tables\" dialog should contain text \"ProductID\"", () => shouldContainText(page, el("Key Columns input in \"Link Tables\" dialog"), "ProductID"));
    await session.step(34, "When user selects \"selection to filter\" in \"Link Type\" input in \"Link Tables\" dialog", () => selectIn(page, "selection to filter", el("\"Link Type\" input in \"Link Tables\" dialog")), undefined, "Link type: what travels, and what it does in the other table");
    await session.step(35, "And user clicks on LINK button in \"Link Tables\" dialog", () => clickOn(page, el("LINK button in \"Link Tables\" dialog")));
    await session.step(36, "Then \"products -> order-details\" text in \"Link Tables\" dialog should be visible", () => shouldBe(page, el("\"products -> order-details\" text in \"Link Tables\" dialog"), "visible"));
    await session.step(37, "When user clicks on CLOSE button in \"Link Tables\" dialog", () => clickOn(page, el("CLOSE button in \"Link Tables\" dialog")));
    await session.step(38, "Then the \"Link Tables\" dialog should close", () => dialogCloses(page, "Link Tables"));
    await session.step(40, "When user clicks on the \"row header 1\" area of grid", () => clickArea(page, "row header 1", el("grid")), undefined, "Select products: the order lines filter to theirs");
    await session.step(41, "And user clicks on the \"row header 3\" area of grid holding Shift", () => clickAreaHolding(page, "row header 3", el("grid"), "Shift"));
    await session.step(42, "Then 3 rows should be selected", () => selectedRowCount(page, 3));
    await session.step(43, "And 94 rows of table \"order-details\" should pass the filter", () => tableFilterCount(page, 94, "order-details"));
    await session.step(45, "When user clicks on the \"row header 5\" area of grid holding Control", () => clickAreaHolding(page, "row header 5", el("grid"), "Control"), undefined, "Ctrl+click adds a product, and its order lines");
    await session.step(46, "Then 4 rows should be selected", () => selectedRowCount(page, 4));
    await session.step(47, "And 104 rows of table \"order-details\" should pass the filter", () => tableFilterCount(page, 104, "order-details"));
    await session.step(48, "And no errors should have been logged", () => noErrors(page));
  });
});
