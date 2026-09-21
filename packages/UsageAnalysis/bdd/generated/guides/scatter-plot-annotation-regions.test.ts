/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/guides/scatter-plot-annotation-regions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/spaces.js';
import '../../bindings/tile-viewer.js';
import '../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaPainted, dragAcrossArea, hasArea, pickFromContextMenu, propertyShouldContain, readingIs, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Draw a custom annotation region on a scatter plot", () => {
  const session = feature(test, "features/guides/scatter-plot-annotation-regions.feature", import.meta.url);
  test("Draw an annotation region on a scatter plot", {tag: ["@guide", "@help:visualize/viewers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a scatter plot viewer with:", () => addViewerWith(page, "scatter plot", [["X","AGE"],["Y","WEIGHT"]]), [["X","AGE"],["Y","WEIGHT"]]);
    await session.step(16, "When user picks \"Tools > Draw Annotation Region\" from the context menu of scatter plot viewer", () => pickFromContextMenu(page, "Tools > Draw Annotation Region", el("scatter plot viewer")));
    await session.step(17, "Then the \"region drawing mode\" reading of scatter plot viewer should be \"true\"", () => readingReads(page, "region drawing mode", el("scatter plot viewer"), "true"));
    await session.step(18, "When user drags across the \"view\" area of scatter plot viewer", () => dragAcrossArea(page, "view", el("scatter plot viewer")));
    await session.step(19, "Then scatter plot viewer should have a \"region 1\" area", () => hasArea(page, el("scatter plot viewer"), "region 1"));
    await session.step(20, "And the \"region 1\" area of scatter plot viewer should be painted", () => areaPainted(page, "region 1", el("scatter plot viewer")));
    await session.step(21, "When user clicks OK button in \"Formula Lines\" dialog", () => clickOn(page, el("OK button in \"Formula Lines\" dialog")));
    await session.step(22, "Then the \"regions shown\" reading of scatter plot viewer should be 1", () => readingIs(page, "regions shown", el("scatter plot viewer"), 1));
    await session.step(23, "And \"annotationRegions\" property of scatter plot viewer should contain \"area\"", () => propertyShouldContain(page, "annotationRegions", el("scatter plot viewer"), "area"));
  });
});
