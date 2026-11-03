/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/grid/grid-extract-styles.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.grid, core#231]
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
import {clickOn, isExpanded} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnsExactly, setGridColumnBackColor} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {columnsSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {closeAllViews, openDataset, switchTableView} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {openTableViewsExactly} from '@datagrok-libraries/bdd/bindings/platform/workspace';
import {areaColor, clickAreaHolding, dragAreaBy, noErrors, pickFromAreaContextMenu, readingsDiffer} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Extract keeps how the columns look in the grid", () => {
  const session = feature(test, "features/viewers/grid/grid-extract-styles.feature", import.meta.url);
  test("Extract keeps how the columns look in the grid", {tag: ["@journey", "@viewers", "@realizes:viewers.grid", "@realizes:core#231"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await run.scenario("Extract keeps the look of several selected columns", async () => {
      await session.step(15, "When user sets the back color of the cells of \"HEIGHT\" column to \"#F5F5A3\"", () => setGridColumnBackColor(page, "cells", "HEIGHT", "#F5F5A3"));
      await session.step(16, "And user sets the back color of the header of \"HEIGHT\" column to \"#904040\"", () => setGridColumnBackColor(page, "header", "HEIGHT", "#904040"));
      await session.step(17, "And user sets the back color of the cells of \"WEIGHT\" column to \"#C8E6C9\"", () => setGridColumnBackColor(page, "cells", "WEIGHT", "#C8E6C9"));
      await session.step(18, "And user drags the \"column resizer HEIGHT\" area of grid by 60 pixels to the right", () => dragAreaBy(page, "column resizer HEIGHT", el("grid"), 60, "right"));
      await session.step(19, "And user clicks on the \"header HEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header HEIGHT", el("grid"), "Control"));
      await session.step(20, "And user clicks on the \"header WEIGHT\" area of grid holding Control", () => clickAreaHolding(page, "header WEIGHT", el("grid"), "Control"));
      await session.step(21, "Then columns \"HEIGHT, WEIGHT\" should be selected", () => columnsSelected(page, "HEIGHT, WEIGHT"));
      await session.step(22, "Given \"Actions\" section in context panel is expanded", () => isExpanded(page, el("\"Actions\" section in context panel")));
      await session.step(23, "When user clicks on \"Extract\" label in context panel", () => clickOn(page, el("\"Extract\" label in context panel")));
      await session.step(24, "Then the open table views should be exactly \"demog-1000, result\"", () => openTableViewsExactly(page, "demog-1000, result"));
      await session.step(25, "When user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(26, "Then the table should have the columns \"HEIGHT, WEIGHT\"", () => columnsExactly(page, "HEIGHT, WEIGHT"));
      await session.step(27, "And the \"cell 1 of HEIGHT\" area of grid should contain the color \"#F5F5A3\"", () => areaColor(page, "cell 1 of HEIGHT", el("grid"), "#F5F5A3"));
      await session.step(28, "And the \"header HEIGHT\" area of grid should contain the color \"#904040\"", () => areaColor(page, "header HEIGHT", el("grid"), "#904040"));
      await session.step(29, "And the \"cell 1 of WEIGHT\" area of grid should contain the color \"#C8E6C9\"", () => areaColor(page, "cell 1 of WEIGHT", el("grid"), "#C8E6C9"));
      await session.step(30, "And the \"column width of HEIGHT\" and \"column width of WEIGHT\" readings of grid should differ", () => readingsDiffer(page, "column width of HEIGHT", "column width of WEIGHT", el("grid")));
      await session.step(31, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Extract keeps the look of a single column", async () => {
      await session.step(35, "When user closes all views", () => closeAllViews(page));
      await session.step(36, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(37, "And user sets the back color of the cells of \"AGE\" column to \"#F5F5A3\"", () => setGridColumnBackColor(page, "cells", "AGE", "#F5F5A3"));
      await session.step(38, "And user picks \"Extract\" from the context menu of the \"header AGE\" area of grid", () => pickFromAreaContextMenu(page, "Extract", "header AGE", el("grid")));
      await session.step(39, "Then the open table views should be exactly \"demog-1000, result\"", () => openTableViewsExactly(page, "demog-1000, result"));
      await session.step(40, "When user switches to the \"result\" table view", () => switchTableView(page, "result"));
      await session.step(41, "Then the table should have the columns \"AGE\"", () => columnsExactly(page, "AGE"));
      await session.step(42, "And the \"cell 1 of AGE\" area of grid should contain the color \"#F5F5A3\"", () => areaColor(page, "cell 1 of AGE", el("grid"), "#F5F5A3"));
      await session.step(43, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
