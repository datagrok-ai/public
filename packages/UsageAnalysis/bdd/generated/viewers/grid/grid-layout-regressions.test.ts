/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-layout-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, public#1625, public#679]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowsRangeSelected, selectedRowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, clickAreaHolding, loadLayout, noErrors, pickFromAreaContextMenu, saveLayoutToServer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {areaLiesBeside, gridPins, stackViewer, switchTab, viewersTabbed} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Grid regressions — what the grid keeps after a layout is applied", () => {
  const session = feature(test, "features/viewers/grid/grid-layout-regressions.feature", import.meta.url);
  test("Pinned columns stay side by side after a layout with a viewer stacked over the grid (public#1625)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#1625"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "When user picks \"Pin > Pin Column\" from the context menu of the \"header USUBJID\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header USUBJID", el("grid")));
    await session.step(18, "And user picks \"Pin > Pin Column\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Pin > Pin Column", "header AGE", el("grid")));
    await session.step(19, "Then the grid should pin the columns \"USUBJID, AGE\"", () => gridPins(page, "USUBJID, AGE"));
    await session.step(20, "And the \"header AGE\" area of grid should lie to the right of the \"header USUBJID\" area", () => areaLiesBeside(page, "header AGE", el("grid"), "right", "header USUBJID"));
    await session.step(21, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(22, "When user stacks scatter plot viewer onto grid as a tab", () => stackViewer(page, el("scatter plot viewer"), el("grid")));
    await session.step(23, "Then the viewers \"demog-1000, Scatter plot\" should be the tabs of one panel", () => viewersTabbed(page, "demog-1000, Scatter plot"));
    await session.step(24, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(25, "And user loads the saved layout", () => loadLayout(page));
    await session.step(26, "And user switches the tabbed panel to the \"demog-1000\" tab", () => switchTab(page, "demog-1000"));
    await session.step(27, "Then the grid should pin the columns \"USUBJID, AGE\"", () => gridPins(page, "USUBJID, AGE"));
    await session.step(28, "And the \"header AGE\" area of grid should lie to the right of the \"header USUBJID\" area", () => areaLiesBeside(page, "header AGE", el("grid"), "right", "header USUBJID"));
    await session.step(29, "And no errors should have been logged", () => noErrors(page));
  });
  test("A row-number click still makes the row current after a layout is applied (public#679)", {tag: ["@viewers", "@realizes:viewers.grid", "@realizes:public#679"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(34, "Given user adds a scatter plot viewer", () => addViewer(page, "scatter plot"));
    await session.step(35, "When user saves the layout of the current table view to the server", () => saveLayoutToServer(page));
    await session.step(36, "And user loads the saved layout", () => loadLayout(page));
    await session.step(37, "Then scatter plot viewer should be visible", () => shouldBe(page, el("scatter plot viewer"), "visible"));
    await session.step(38, "When user clicks on the \"row header 5\" area of grid", () => clickArea(page, "row header 5", el("grid")));
    await session.step(39, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(40, "When user clicks on the \"row header 7\" area of grid holding Shift", () => clickAreaHolding(page, "row header 7", el("grid"), "Shift"));
    await session.step(41, "Then rows 5 to 7 should be selected", () => rowsRangeSelected(page, 5, 7));
    await session.step(42, "And 3 rows should be selected", () => selectedRowCount(page, 3));
    await session.step(43, "And no errors should have been logged", () => noErrors(page));
  });
});
