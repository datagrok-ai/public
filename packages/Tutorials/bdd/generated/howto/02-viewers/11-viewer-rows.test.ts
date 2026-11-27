/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/11-viewer-rows.feature
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
import {viewerShare} from '../../../bindings/howto-viewers.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, enterInto, selectIn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterPassesAll, selectRowsByHeaders, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, contextPanelShows, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {hoverEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/11-viewer-rows.feature", import.meta.url);
  test("Show only some rows in a viewer", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And simple mode is off", () => simpleModeOff(page));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(19, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(20, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(21, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","WEIGHT"],["Y","HEIGHT"]]), [["X","WEIGHT"],["Y","HEIGHT"]]);
    await session.step(24, "And scatter plot viewer takes 60% of the space it shares with the grid", () => viewerShare(page, el("scatter plot viewer"), 60));
    await session.step(25, "Then scatter plot viewer should show 872 rows", () => showsRows(page, el("scatter plot viewer"), 872));
    await session.step(26, "When the video says \"A viewer can show fewer rows than its table\"", () => videoSays(page, "A viewer can show fewer rows than its table"));
    await session.step(28, "And user hovers over empty plot space of scatter plot viewer", () => hoverEmptySpace(page, el("scatter plot viewer")), undefined, "The gear opens the viewer's settings");
    await session.step(29, "And user clicks on settings icon of scatter plot viewer", () => clickOn(page, el("settings icon of scatter plot viewer")));
    await session.step(30, "Then the context panel should show \"Scatter plot\"", () => contextPanelShows(page, "Scatter plot"));
    await session.step(32, "When user enters \"${AGE} > 44\" into \"Filter\" property in context panel", () => enterInto(page, "${AGE} > 44", el("\"Filter\" property in context panel")), undefined, "Data > Filter: a formula picks the rows this viewer shows");
    await session.step(33, "Then scatter plot viewer should show 466 rows", () => showsRows(page, el("scatter plot viewer"), 466));
    await session.step(35, "And all rows should pass the filter", () => filterPassesAll(page), undefined, "The grid keeps all 1,000 rows");
    await session.step(37, "When user selects \"Selected\" in \"Row Source\" property in context panel", () => selectIn(page, "Selected", el("\"Row Source\" property in context panel")), undefined, "Row Source: here, only the selected rows");
    await session.step(38, "Then scatter plot viewer should show 0 rows", () => showsRows(page, el("scatter plot viewer"), 0));
    await session.step(40, "When user selects rows 1 to 25 by their row headers", () => selectRowsByHeaders(page, 1, 25), undefined, "Select rows in the grid: the plot shows those over 44");
    await session.step(41, "Then 25 rows should be selected", () => selectedRowCount(page, 25));
    await session.step(42, "And scatter plot viewer should show 8 rows", () => showsRows(page, el("scatter plot viewer"), 8));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
