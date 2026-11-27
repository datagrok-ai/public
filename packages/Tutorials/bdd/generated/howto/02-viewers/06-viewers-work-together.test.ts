/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/02-viewers/06-viewers-work-together.feature
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
import {clickOn, videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {filterPasses, onlyOfSelected, selectedRowCount, someSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, addViewerWith, clickArea, dragSelectionBetweenAreas, hasArea, hoverArea, noErrors, readingIs, repainted, showsRows, someHighlight, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Viewers", () => {
  const session = feature(test, "features/howto/02-viewers/06-viewers-work-together.feature", import.meta.url);
  test("Viewers work together", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And simple mode is off", () => simpleModeOff(page));
    await session.step(15, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(16, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(17, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(18, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(19, "And user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(20, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"]]), [["Value","AGE"]]);
    await session.step(22, "And user adds a bar chart viewer", () => addViewer(page, "bar chart"));
    await session.step(23, "When the video says \"The grid and all viewers of a table show the same rows\"", () => videoSays(page, "The grid and all viewers of a table show the same rows"));
    await session.step(25, "And user takes a snapshot of scatter plot viewer", () => takeSnapshot(page, el("scatter plot viewer")), undefined, "Point at a bin: its rows light up on the other viewers");
    await session.step(26, "And user hovers over the \"bin 8\" area of histogram viewer", () => hoverArea(page, "bin 8", el("histogram viewer")));
    await session.step(27, "Then scatter plot viewer should have repainted", () => repainted(page, el("scatter plot viewer")));
    await session.step(29, "When user clicks on the \"cell 10 of AGE\" area of grid", () => clickArea(page, "cell 10 of AGE", el("grid")), undefined, "Click a row to make it current everywhere");
    await session.step(30, "Then row 10 should be current", () => currentRowIs(page, 10));
    await session.step(31, "And histogram viewer should have a \"current row marker\" area", () => hasArea(page, el("histogram viewer"), "current row marker"));
    await session.step(33, "When user drags a selection box from the \"bin 4\" area to the \"bin 7\" area of histogram viewer", () => dragSelectionBetweenAreas(page, "bin 4", "bin 7", el("histogram viewer")), undefined, "Shift+drag across bins selects their rows, shown everywhere");
    await session.step(34, "Then some rows should be selected", () => someSelected(page));
    await session.step(35, "And scatter plot viewer should show a selection highlight", () => someHighlight(page, el("scatter plot viewer")));
    await session.step(36, "And bar chart viewer should show a selection highlight", () => someHighlight(page, el("bar chart viewer")));
    await session.step(38, "When user clicks on the \"bar RA\" area of bar chart viewer", () => clickArea(page, "bar RA", el("bar chart viewer")), undefined, "Click a bar to select its rows");
    await session.step(39, "Then only rows where \"DIS_POP\" is \"RA\" should be selected", () => onlyOfSelected(page, "DIS_POP", "RA"));
    await session.step(40, "And 2550 rows should be selected", () => selectedRowCount(page, 2550));
    await session.step(42, "When user clicks on filter icon in toolbar", () => clickOn(page, el("filter icon in toolbar")), undefined, "A filter applies to the grid and every viewer");
    await session.step(43, "And user clicks on the \"category Psoriasis of DIS_POP\" area of filter panel", () => clickArea(page, "category Psoriasis of DIS_POP", el("filter panel")));
    await session.step(44, "Then 1209 rows should pass the filter", () => filterPasses(page, 1209));
    await session.step(45, "And histogram viewer should show 1209 rows", () => showsRows(page, el("histogram viewer"), 1209));
    await session.step(46, "And the \"bars\" reading of bar chart viewer should be 1", () => readingIs(page, "bars", el("bar chart viewer"), 1));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
});
