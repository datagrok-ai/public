/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/03-more-viewers/02-pie-chart.feature
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
import {check, clickOn, expand, hoverOver, selectIn, shouldBe, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneShown} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {clickArea, hoverArea, noErrors, propertyShouldBe, readingIs, viewerAdded} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInColumnSelector} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · More viewers", () => {
  const session = feature(test, "features/howto/03-more-viewers/02-pie-chart.feature", import.meta.url);
  test("Show shares with a pie chart", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(13, "And simple mode is off", () => simpleModeOff(page));
    await session.step(14, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(15, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(16, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(17, "And the toolbox pane is shown", () => toolboxPaneShown(page));
    await session.step(18, "When the video says \"A pie chart shows each category's share of the rows\"", () => videoSays(page, "A pie chart shows each category's share of the rows"));
    await session.step(19, "And user clicks on pie chart icon on toolbox", () => clickOn(page, el("pie chart icon on toolbox")));
    await session.step(20, "Then pie chart viewer should be added to the open tableview", () => viewerAdded(page, "pie chart"));
    await session.step(22, "When user clicks on toolbox tab", () => clickOn(page, el("toolbox tab")), undefined, "Hide the toolbox to give the chart room");
    await session.step(23, "Then toolbox should be hidden", () => shouldBe(page, el("toolbox"), "hidden"));
    await session.step(25, "When user picks \"SEVERITY\" in the \"category\" column selector of pie chart viewer", () => pickInColumnSelector(page, "SEVERITY", "category", el("pie chart viewer")), undefined, "Pick the category right on the chart");
    await session.step(26, "Then the \"slices\" reading of pie chart viewer should be 5", () => readingIs(page, "slices", el("pie chart viewer"), 5));
    await session.step(28, "When user hovers over the \"slice High\" area of pie chart viewer", () => hoverArea(page, "slice High", el("pie chart viewer")), undefined, "Hover a slice for its row count");
    await session.step(30, "And user clicks on the \"slice High\" area of pie chart viewer", () => clickArea(page, "slice High", el("pie chart viewer")), undefined, "Click a slice to select its rows");
    await session.step(31, "Then 650 rows should be selected", () => selectedRowCount(page, 650));
    await session.step(33, "When user hovers over pie chart viewer", () => hoverOver(page, el("pie chart viewer")), undefined, "The rest is in the settings");
    await session.step(34, "And user clicks on settings icon of pie chart viewer", () => clickOn(page, el("settings icon of pie chart viewer")));
    await session.step(35, "Then the context panel should show \"Pie chart\"", () => contextPanelShows(page, "Pie chart"));
    await session.step(37, "When user selects \"desc\" in \"Pie Sort Order\" property in context panel", () => selectIn(page, "desc", el("\"Pie Sort Order\" property in context panel")), undefined, "Data: sort the slices, largest first");
    await session.step(38, "Then the \"start angle of None\" reading of pie chart viewer should be 0", () => readingIs(page, "start angle of None", el("pie chart viewer"), 0));
    await session.step(40, "When user expands \"Style\" category in context panel", () => expand(page, el("\"Style\" category in context panel")), undefined, "Style: turn the pie into a donut");
    await session.step(41, "And user selects \"Donut\" in \"Mode\" property in context panel", () => selectIn(page, "Donut", el("\"Mode\" property in context panel")));
    await session.step(42, "Then \"Mode\" property of pie chart viewer should be \"Donut\"", () => propertyShouldBe(page, "Mode", el("pie chart viewer"), "Donut"));
    await session.step(44, "When user expands \"Misc\" category in context panel", () => expand(page, el("\"Misc\" category in context panel")), undefined, "Misc: show the row count on each slice");
    await session.step(45, "And user checks \"Show Value\" property in context panel", () => check(page, el("\"Show Value\" property in context panel")));
    await session.step(46, "Then no errors should have been logged", () => noErrors(page));
  });
});
