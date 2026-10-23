/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/line-chart/line-chart-mouse-over-and-overview.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.line-chart]
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
import {mouseOverRowIs} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, areaShorter, hasArea, hasNoArea, hoverArea, moreInk, noErrors, notRepainted, pickFromContextMenu, pointerAway, readingIs, reportsNoError, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Line chart mouse-over row line and the stacked bar overview", () => {
  const session = feature(test, "features/viewers/line-chart/line-chart-mouse-over-and-overview.feature", import.meta.url);
  test("Line chart mouse-over row line and the stacked bar overview", {tag: ["@journey", "@viewers", "@realizes:viewers.line-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(13, "And user adds a line chart viewer with:", () => addViewerWith(page, "line chart", [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]), [["xColumnName","CAST Idea ID"],["yColumnNames","Chemical Space X"]]);
    await session.step(16, "Then the \"rows shown\" reading of line chart viewer should be 100", () => readingIs(page, "rows shown", el("line chart viewer"), 100));
    await session.step(17, "And line chart viewer should report no error", () => reportsNoError(page, el("line chart viewer")));
    await run.scenario("Overview Stacked Bar Chart draws the overview strip under the plot", async () => {
      await session.step(20, "Then line chart viewer should not have an \"overview\" area", () => hasNoArea(page, el("line chart viewer"), "overview"));
      await session.step(21, "When user picks \"Overview > Stacked Bar Chart\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Overview > Stacked Bar Chart", el("line chart viewer")));
      await session.step(22, "Then line chart viewer should have an \"overview\" area", () => hasArea(page, el("line chart viewer"), "overview"));
      await session.step(23, "And the \"overview\" area of line chart viewer should be painted", () => areaPainted(page, "overview", el("line chart viewer")));
      await session.step(24, "And the \"plot\" area of line chart viewer should be shorter than before", () => areaShorter(page, "plot", el("line chart viewer")));
      await session.step(25, "When user picks \"Overview > None\" from the context menu of line chart viewer", () => pickFromContextMenu(page, "Overview > None", el("line chart viewer")));
      await session.step(26, "Then line chart viewer should not have an \"overview\" area", () => hasNoArea(page, el("line chart viewer"), "overview"));
      await session.step(27, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A grid hover draws the mouse-over row line, and Show Mouse-over Row Line switches it off", async () => {
      await session.step(30, "When user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(31, "And user hovers over the \"cell 3 of Id\" area of grid", () => hoverArea(page, "cell 3 of Id", el("grid")));
      await session.step(32, "Then the mouse-over row of the table should be 3", () => mouseOverRowIs(page, 3));
      await session.step(33, "And line chart viewer should have more ink than before", () => moreInk(page, el("line chart viewer")));
      await session.step(34, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(35, "And user sets \"showMouseOverRowLine\" property of line chart viewer to \"false\"", () => setProperty(page, "showMouseOverRowLine", el("line chart viewer"), "false"));
      await session.step(36, "And user takes a snapshot of line chart viewer", () => takeSnapshot(page, el("line chart viewer")));
      await session.step(37, "And user hovers over the \"cell 3 of Id\" area of grid", () => hoverArea(page, "cell 3 of Id", el("grid")));
      await session.step(38, "Then the mouse-over row of the table should be 3", () => mouseOverRowIs(page, 3));
      await session.step(39, "And line chart viewer should not have repainted", () => notRepainted(page, el("line chart viewer")));
      await session.step(40, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(41, "And user sets \"showMouseOverRowLine\" property of line chart viewer to \"true\"", () => setProperty(page, "showMouseOverRowLine", el("line chart viewer"), "true"));
      await session.step(42, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
