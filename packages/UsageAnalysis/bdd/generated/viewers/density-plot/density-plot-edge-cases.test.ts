/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/density-plot/density-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.density-plot, core#1048]
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
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, renameColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragZoomOverArea, hasArea, hoverArea, readingAtLeast, readingIs, readingReads, readingsEqual, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Density plot edge cases", () => {
  const session = feature(test, "features/viewers/density-plot/density-plot-edge-cases.feature", import.meta.url);
  test("A click on a rectangular bin selects exactly the rows counted in it", {tag: ["@viewers", "@realizes:viewers.density-plot", "@edge-cases", "@realizes:core#1048"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["binShape","rectangle"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["binShape","rectangle"]]);
    await session.step(23, "Then the \"bin shape\" reading of density plot viewer should be \"rectangle\"", () => readingReads(page, "bin shape", el("density plot viewer"), "rectangle"));
    await session.step(24, "And density plot viewer should have a \"densest bin\" area", () => hasArea(page, el("density plot viewer"), "densest bin"));
    await session.step(25, "When user clears the row selection", () => clearSelection(page));
    await session.step(26, "And user clicks on the \"densest bin\" area of density plot viewer", () => clickArea(page, "densest bin", el("density plot viewer")));
    await session.step(27, "Then the \"rows in densest bin\" reading of density plot viewer should be at least 1", () => readingAtLeast(page, "rows in densest bin", el("density plot viewer"), 1));
    await session.step(28, "And the \"rows selected\" and \"rows in densest bin\" readings of density plot viewer should be the same", () => readingsEqual(page, "rows selected", "rows in densest bin", el("density plot viewer")));
  });
  test("Changing Bins after a zoom bins the same rows while Bin To Range is off", {tag: ["@viewers", "@realizes:viewers.density-plot", "@edge-cases", "@realizes:core#1048"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(34, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(35, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","AGE"],["yColumnName","WEIGHT"],["binShape","rectangle"]]), [["xColumnName","AGE"],["yColumnName","WEIGHT"],["binShape","rectangle"]]);
    await session.step(39, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
    await session.step(40, "And the \"bin to range\" reading of density plot viewer should be \"false\"", () => readingReads(page, "bin to range", el("density plot viewer"), "false"));
    await session.step(41, "When user drags a zoom box over the \"view\" area of density plot viewer", () => dragZoomOverArea(page, "view", el("density plot viewer")));
    await session.step(42, "And user sets \"bins\" property of density plot viewer to \"30\"", () => setProperty(page, "bins", el("density plot viewer"), "30"));
    await session.step(43, "Then the \"rows shown\" reading of density plot viewer should be 1000", () => readingIs(page, "rows shown", el("density plot viewer"), 1000));
  });
  test("A logarithmic X axis bins only the positive values", {tag: ["@viewers", "@realizes:viewers.density-plot", "@edge-cases", "@realizes:core#1048"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(50, "Given user opens a table \"log-x\" with:", () => openTableOf(page, "log-x", [["X","Y"],["-5","1"],["-1","2"],["0","3"],["2","4"],["10","5"],["100","6"]]), [["X","Y"],["-5","1"],["-1","2"],["0","3"],["2","4"],["10","5"],["100","6"]]);
    await session.step(58, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","X"],["yColumnName","Y"]]), [["xColumnName","X"],["yColumnName","Y"]]);
    await session.step(61, "Then the \"rows shown\" reading of density plot viewer should be 6", () => readingIs(page, "rows shown", el("density plot viewer"), 6));
    await session.step(62, "When user sets \"xAxisType\" property of density plot viewer to \"logarithmic\"", () => setProperty(page, "xAxisType", el("density plot viewer"), "logarithmic"));
    await session.step(63, "Then the \"rows shown\" reading of density plot viewer should be 3", () => readingIs(page, "rows shown", el("density plot viewer"), 3));
  });
  test("A column name with markup characters is shown in the bin tooltip as typed", {tag: ["@viewers", "@realizes:viewers.density-plot", "@edge-cases", "@realizes:core#1048"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(12, "Given user is logged in", () => loggedIn(page));
    await session.step(69, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(70, "And user renames \"AGE\" column to \"<b>X</b>\"", () => renameColumn(page, "AGE", "<b>X</b>"));
    await session.step(71, "And user adds a density plot viewer with:", () => addViewerWith(page, "density plot", [["xColumnName","<b>X</b>"],["yColumnName","WEIGHT"]]), [["xColumnName","<b>X</b>"],["yColumnName","WEIGHT"]]);
    await session.step(74, "When user hovers over the \"densest bin\" area of density plot viewer", () => hoverArea(page, "densest bin", el("density plot viewer")));
    await session.step(75, "Then tooltip should contain text \"WEIGHT\"", () => shouldContainText(page, el("tooltip"), "WEIGHT"));
    await session.step(76, "And tooltip should contain text \"<b>X</b>\"", () => shouldContainText(page, el("tooltip"), "<b>X</b>"));
  });
});
