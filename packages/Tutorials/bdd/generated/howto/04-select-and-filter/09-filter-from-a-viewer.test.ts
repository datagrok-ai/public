/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/howto/04-select-and-filter/09-filter-from-a-viewer.feature
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
import {videoSays} from '@datagrok-libraries/bdd/bindings/common/steps';
import {filterIsExactlyCategory, filterPasses, filterPassesAll, filterPassesFewer} from '@datagrok-libraries/bdd/bindings/platform/data';
import {browsePanelClosed, contextPanelClosed, openDataset, simpleModeOff, toolboxPaneHidden} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, closeContextMenu, noErrors, pickFromContextMenu, pointerAway, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {clickEmptySpace, dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("How-to · Select and filter", () => {
  const session = feature(test, "features/howto/04-select-and-filter/09-filter-from-a-viewer.feature", import.meta.url);
  test("Filter from a viewer", {tag: ["@guide", "@tutorial"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And simple mode is off", () => simpleModeOff(page));
    await session.step(16, "And user opens demog dataset", () => openDataset(page, ds("demog")));
    await session.step(17, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","DIS_POP"]]), [["Split","DIS_POP"]]);
    await session.step(19, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(22, "And the toolbox pane is hidden", () => toolboxPaneHidden(page));
    await session.step(23, "And the context panel is closed", () => contextPanelClosed(page));
    await session.step(24, "And the browse panel is closed", () => browsePanelClosed(page));
    await session.step(25, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(27, "When user picks \"On Click > Filter\" from the context menu of bar chart viewer", () => pickFromContextMenu(page, "On Click > Filter", el("bar chart viewer")), undefined, "A bar chart can filter on click instead of selecting");
    await session.step(28, "And user closes the context menu", () => closeContextMenu(page));
    await session.step(29, "Then \"On Click\" property of bar chart viewer should be \"Filter\"", () => propertyShouldBe(page, "On Click", el("bar chart viewer"), "Filter"));
    await session.step(31, "When user clicks on the \"bar RA\" area of bar chart viewer", () => clickArea(page, "bar RA", el("bar chart viewer")), undefined, "Click a bar to keep only its rows");
    await session.step(32, "Then 2550 rows should pass the filter", () => filterPasses(page, 2550));
    await session.step(33, "And the filter should pass exactly the rows where \"DIS_POP\" is \"RA\"", () => filterIsExactlyCategory(page, "DIS_POP", "RA"));
    await session.step(34, "When user clicks on the \"bar UC\" area of bar chart viewer", () => clickArea(page, "bar UC", el("bar chart viewer")));
    await session.step(35, "Then 728 rows should pass the filter", () => filterPasses(page, 728));
    await session.step(37, "When user clicks on empty plot space of bar chart viewer", () => clickEmptySpace(page, el("bar chart viewer")), undefined, "Click empty space to let go");
    await session.step(38, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(40, "When user drags the max handle of the \"x\" range slider of scatter plot viewer by -350 pixels", () => dragRangeHandle(page, "max", "x", el("scatter plot viewer"), -350), undefined, "Drag the end of an axis slider to zoom in");
    await session.step(41, "Then fewer than 4000 rows should pass the filter", () => filterPassesFewer(page, 4000));
    await session.step(42, "When user moves the pointer away from scatter plot viewer", () => pointerAway(page, el("scatter plot viewer")));
    await session.step(43, "And the video says \"Zooming a scatter plot filters the table to the points in view\"", () => videoSays(page, "Zooming a scatter plot filters the table to the points in view"));
    await session.step(45, "And user picks \"Reset View\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Reset View", el("scatter plot viewer")), undefined, "Reset View brings every row back");
    await session.step(46, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(47, "And no errors should have been logged", () => noErrors(page));
  });
});
