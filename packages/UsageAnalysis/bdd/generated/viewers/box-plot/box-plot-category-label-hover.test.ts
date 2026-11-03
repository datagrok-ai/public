/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-category-label-hover.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, hoverArea, noErrors, notRepainted, pointerAway, repaintedBy, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot — hovering a category label highlights its row group", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-category-label-hover.feature", import.meta.url);
  test("Box plot — hovering a category label highlights its row group", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(12, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]), [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]);
    await run.scenario("Hovering a category's axis label highlights its row group, and Show Mouse Over Row Group switches it off", async () => {
      await session.step(18, "When user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(19, "And user hovers over the \"category Caucasian\" area of box plot viewer", () => hoverArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(20, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(21, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(22, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "false"));
      await session.step(23, "And user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(24, "And user hovers over the \"category Caucasian\" area of box plot viewer", () => hoverArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(25, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(26, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(27, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "true"));
      await session.step(28, "Then no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
