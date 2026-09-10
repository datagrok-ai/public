/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-colors.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
--- */
import {test} from '@playwright/test';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, areaPainted, areaRepainted, hasArea, noErrors, readingBetween, readingIs, readingReads, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map colouring", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-colors.feature", import.meta.url);
  test("Heat map colouring", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(14, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(15, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
    await session.step(16, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
    await session.step(17, "And the \"column AGE\" area of heat map viewer should be painted", () => areaPainted(page, "column AGE", el("heat map viewer")));
    await run.scenario("Global Color Scaling recolours every numerical column", async () => {
      await session.step(20, "Then the \"global color scaling\" reading of heat map viewer should be \"false\"", () => readingReads(page, "global color scaling", el("heat map viewer"), "false"));
      await session.step(21, "When user sets \"globalColorScaling\" property of heat map viewer to \"true\"", () => setProperty(page, "globalColorScaling", el("heat map viewer"), "true"));
      await session.step(22, "Then the \"global color scaling\" reading of heat map viewer should be \"true\"", () => readingReads(page, "global color scaling", el("heat map viewer"), "true"));
      await session.step(23, "And the \"column AGE\" area of heat map viewer should have repainted", () => areaRepainted(page, "column AGE", el("heat map viewer")));
      await session.step(24, "And the \"column WEIGHT\" area of heat map viewer should have repainted", () => areaRepainted(page, "column WEIGHT", el("heat map viewer")));
      await session.step(25, "And the \"rows shown\" reading of heat map viewer should be 1000", () => readingIs(page, "rows shown", el("heat map viewer"), 1000));
      await session.step(26, "When user sets \"globalColorScaling\" property of heat map viewer to \"false\"", () => setProperty(page, "globalColorScaling", el("heat map viewer"), "false"));
      await session.step(27, "Then the \"global color scaling\" reading of heat map viewer should be \"false\"", () => readingReads(page, "global color scaling", el("heat map viewer"), "false"));
      await session.step(28, "And the \"column AGE\" area of heat map viewer should have repainted", () => areaRepainted(page, "column AGE", el("heat map viewer")));
      await session.step(29, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Heatmap Colors off stops filling the cells with colour (GROK-20619)", async () => {
      await session.step(39, "Then the \"heatmap colors\" reading of heat map viewer should be \"true\"", () => readingReads(page, "heatmap colors", el("heat map viewer"), "true"));
      await session.step(40, "When user sets \"heatmapColors\" property of heat map viewer to \"false\"", () => setProperty(page, "heatmapColors", el("heat map viewer"), "false"));
      await session.step(41, "Then the \"heatmap colors\" reading of heat map viewer should be \"false\"", () => readingReads(page, "heatmap colors", el("heat map viewer"), "false"));
      await session.step(42, "And the \"column AGE\" area of heat map viewer should have repainted", () => areaRepainted(page, "column AGE", el("heat map viewer")));
    }, {knownFailure: true});
    run.finish();
  });
});
