/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/12-heat-map.feature
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
import {clickOn, expand, hoverOver, shouldBe, videoSays, videoSaysAt} from '@datagrok-libraries/bdd/bindings/common/steps';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {doubleClickArea, hoverArea, noErrors, readingIs, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {toggleInColumnList} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/12-heat-map.feature", import.meta.url);
  test("See a whole table at once with a heat map", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And simple mode is off", () => simpleModeOff(page));
    await session.step(13, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(14, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(15, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(16, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(17, "When the video says \"A heat map colors every value: high red, low blue\"", () => videoSays(page, "A heat map colors every value: high red, low blue"));
    await session.step(18, "And user clicks on heat map icon on toolbox", () => clickOn(page, el("heat map icon on toolbox")));
    await session.step(19, "Then heat map viewer should be added to the open tableview", () => viewerAdded(page, "heat map"));
    await session.step(21, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the map room");
    await session.step(22, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(24, "When user double-clicks on the \"header price\" area of heat map viewer", () => doubleClickArea(page, "header price", el("heat map viewer")), undefined, "Double-click a header to sort the rows by it");
    await session.step(25, "Then the \"sort column\" reading of heat map viewer should be \"price\"", () => readingReads(page, "sort column", el("heat map viewer"), "price"));
    await session.step(26, "When the video says \"Sorted by price: the big, fast cars line up at the top\" at the top", () => videoSaysAt(page, "Sorted by price: the big, fast cars line up at the top", "top"));
    await session.step(27, "And user hovers over the \"header horsepower\" area of heat map viewer", () => hoverArea(page, "header horsepower", el("heat map viewer")));
    await session.step(29, "And user hovers over heat map viewer", () => hoverOver(page, el("heat map viewer")), undefined, "The rest is in the settings");
    await session.step(30, "And user clicks on settings icon of heat map viewer", () => clickOn(page, el("settings icon of heat map viewer")));
    await session.step(32, "And user expands \"Columns\" category in context panel", () => expand(page, el("\"Columns\" category in context panel")), undefined, "Columns: hide the ones you don't need");
    await session.step(33, "And user clicks on first button in \"Columns\" property in context panel", () => clickOn(page, el("first button in \"Columns\" property in context panel")));
    await session.step(34, "And user toggles the \"diesel\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "diesel", el("Order or Hide Columns dialog")));
    await session.step(35, "And user toggles the \"turbo\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "turbo", el("Order or Hide Columns dialog")));
    await session.step(36, "And user toggles the \"two.doors\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "two.doors", el("Order or Hide Columns dialog")));
    await session.step(37, "And user toggles the \"hatchback\" column in the column list of Order or Hide Columns dialog", () => toggleInColumnList(page, "hatchback", el("Order or Hide Columns dialog")));
    await session.step(38, "And user clicks on CLOSE button in Order or Hide Columns dialog", () => clickOn(page, el("CLOSE button in Order or Hide Columns dialog")));
    await session.step(39, "Then the \"columns shown\" reading of heat map viewer should be 13", () => readingIs(page, "columns shown", el("heat map viewer"), 13));
    await session.step(40, "And no errors should have been logged", () => noErrors(page));
  });
});
