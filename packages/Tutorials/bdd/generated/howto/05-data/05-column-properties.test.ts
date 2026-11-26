/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/05-data/05-column-properties.feature
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
import {clickOn, enterInto, expand, hoverText, pressKey, selectIn, shouldBe, shouldContainText, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnType, columnUnits, displayedInRow, hasColumn} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, dialogCloses, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, pickFromAreaContextMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {gridShowsColumns} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Data", () => {
  const session = feature(test, "features/howto/05-data/05-column-properties.feature", import.meta.url);
  test("Column properties", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(17, "And simple mode is off", () => simpleModeOff(page));
    await session.step(18, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(19, "And the grid shows only the columns \"model, diesel, turbo, two.doors, horsepower, city.mpg, price\"", () => gridShowsColumns(page, "model, diesel, turbo, two.doors, horsepower, city.mpg, price"));
    await session.step(20, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(21, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(22, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(24, "When user clicks on the \"cell 1 of price\" area of grid", () => clickArea(page, "cell 1 of price", el("grid")), undefined, "F2 opens the properties of the current column");
    await session.step(25, "And user presses F2", () => pressKey(page, "F2"));
    await session.step(26, "Then \"price\" dialog should be visible", () => shouldBe(page, el("\"price\" dialog"), "visible"));
    await session.step(28, "When user enters \"Price\" into \"New name:\" input in \"price\" dialog", () => enterInto(page, "Price", el("\"New name:\" input in \"price\" dialog")), undefined, "Rename it");
    await session.step(30, "And user clicks on first \"Add row\" button in \"price\" dialog", () => clickOn(page, el("first \"Add row\" button in \"price\" dialog")), undefined, "Tags: a format, units, anything to keep with the column");
    await session.step(31, "And user enters \"format\" into \"Key 1\" text input in \"price\" dialog", () => enterInto(page, "format", el("\"Key 1\" text input in \"price\" dialog")));
    await session.step(32, "And user enters \"#,##0\" into \"Value 1\" text input in \"price\" dialog", () => enterInto(page, "#,##0", el("\"Value 1\" text input in \"price\" dialog")));
    await session.step(33, "And user clicks on first \"Add row\" button in \"price\" dialog", () => clickOn(page, el("first \"Add row\" button in \"price\" dialog")));
    await session.step(34, "And user enters \"units\" into \"Key 1\" text input in \"price\" dialog", () => enterInto(page, "units", el("\"Key 1\" text input in \"price\" dialog")));
    await session.step(35, "And user enters \"USD\" into \"Value 1\" text input in \"price\" dialog", () => enterInto(page, "USD", el("\"Value 1\" text input in \"price\" dialog")));
    await session.step(36, "And user clicks on OK button in \"price\" dialog", () => clickOn(page, el("OK button in \"price\" dialog")));
    await session.step(37, "Then the \"price\" dialog should close", () => dialogCloses(page, "price"));
    await session.step(38, "And the table should have a column \"Price\"", () => hasColumn(page, "Price"));
    await session.step(39, "And \"Price\" column should have units \"USD\"", () => columnUnits(page, "Price", "USD"));
    await session.step(40, "And the \"Price\" cell of row 1 should be displayed as \"16,500\"", () => displayedInRow(page, "Price", 1, "16,500"));
    await session.step(42, "When user clicks on \"Context Panel\" status bar toggle", () => clickOn(page, el("\"Context Panel\" status bar toggle")), undefined, "The tags are listed in the column's Details too");
    await session.step(43, "And user clicks on the \"header Price\" area of grid", () => clickArea(page, "header Price", el("grid")));
    await session.step(44, "Then the context panel should show \"Price\"", () => contextPanelShows(page, "Price"));
    await session.step(45, "When user expands Details pane in context panel", () => expand(page, el("Details pane in context panel")));
    await session.step(46, "Then Details pane in context panel should contain text \"USD\"", () => shouldContainText(page, el("Details pane in context panel"), "USD"));
    await session.step(47, "When user hovers over the text \"USD\" in context panel", () => hoverText(page, "USD", el("context panel")));
    await session.step(49, "When user picks \"Change Type...\" from the context menu of the \"header horsepower\" area of grid", () => pickFromAreaContextMenu(page, "Change Type...", "header horsepower", el("grid")), undefined, "Change the type from the column's menu");
    await session.step(50, "Then \"Change column type\" dialog should be visible", () => shouldBe(page, el("\"Change column type\" dialog"), "visible"));
    await session.step(51, "When user selects \"double\" in \"New type\" input in \"Change column type\" dialog", () => selectIn(page, "double", el("\"New type\" input in \"Change column type\" dialog")));
    await session.step(52, "And user clicks on OK button in \"Change column type\" dialog", () => clickOn(page, el("OK button in \"Change column type\" dialog")));
    await session.step(53, "Then the \"Change column type\" dialog should close", () => dialogCloses(page, "Change column type"));
    await session.step(54, "And \"horsepower\" column should have type \"double\"", () => columnType(page, "horsepower", "double"));
    await session.step(55, "And the \"horsepower\" cell of row 1 should be displayed as \"154.00\"", () => displayedInRow(page, "horsepower", 1, "154.00"));
    await session.step(56, "When the video says \"horsepower is a double now: it can hold fractions\"", () => videoSays(page, "horsepower is a double now: it can hold fractions"));
    await session.step(57, "And user hovers over the \"left edge of cell 1 of horsepower\" area of grid", () => hoverArea(page, "left edge of cell 1 of horsepower", el("grid")));
    await session.step(58, "Then no errors should have been logged", () => noErrors(page));
  });
});
