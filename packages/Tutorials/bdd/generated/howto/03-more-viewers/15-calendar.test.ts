/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/15-calendar.feature
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
import {clickOn, hoverOver, selectIn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPasses, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, noErrors, readingReads, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/15-calendar.feature", import.meta.url);
  test("See dates on a calendar", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(20, "When the video says \"A calendar shows how many rows fall on each day\"", () => videoSays(page, "A calendar shows how many rows fall on each day"));
    await session.step(21, "And user clicks on calendar icon on toolbox", () => clickOn(page, el("calendar icon on toolbox")));
    await session.step(22, "Then calendar viewer should be added to the open tableview", () => viewerAdded(page, "calendar"));
    await session.step(23, "And the \"date column\" reading of calendar viewer should be \"STARTED\"", () => readingReads(page, "date column", el("calendar viewer"), "STARTED"));
    await session.step(25, "When user clicks on the \"month 1990-01\" area of calendar viewer", () => clickArea(page, "month 1990-01", el("calendar viewer")), undefined, "Click a month to select its rows");
    await session.step(26, "Then 218 rows should be selected", () => selectedRowCount(page, 218));
    await session.step(28, "When user clicks on the \"weekday Sunday\" area of calendar viewer", () => clickArea(page, "weekday Sunday", el("calendar viewer")), undefined, "or a weekday to select all of those days");
    await session.step(29, "Then 837 rows should be selected", () => selectedRowCount(page, 837));
    await session.step(31, "When user hovers over calendar viewer", () => hoverOver(page, el("calendar viewer")), undefined, "The settings hold the rest");
    await session.step(32, "And user clicks on settings icon of calendar viewer", () => clickOn(page, el("settings icon of calendar viewer")));
    await session.step(33, "Then the context panel should show \"Calendar\"", () => contextPanelShows(page, "Calendar"));
    await session.step(35, "When user selects \"Filter\" in \"On Click\" property in context panel", () => selectIn(page, "Filter", el("\"On Click\" property in context panel")), undefined, "Data: make a click filter the table instead");
    await session.step(36, "And user clicks on the \"month 1989-12\" area of calendar viewer", () => clickArea(page, "month 1989-12", el("calendar viewer")));
    await session.step(37, "Then 264 rows should pass the filter", () => filterPasses(page, 264));
    await session.step(38, "And no errors should have been logged", () => noErrors(page));
  });
});
