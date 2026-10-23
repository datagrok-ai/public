/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/statistics/statistics-back-color.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.stats-viewer]
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
import {addViewer, areaColor, areaNotColor, noErrors, readingIs, readingReads, resizeTo, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Statistics — Back Color paints the backdrop of the statistics table", () => {
  const session = feature(test, "features/viewers/statistics/statistics-back-color.feature", import.meta.url);
  test("Back Color repaints the backdrop, and white takes it away again", {tag: ["@viewers", "@realizes:viewers.stats-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a statistics viewer", () => addViewer(page, "statistics"));
    await session.step(14, "And user resizes statistics viewer to 900 by 400", () => resizeTo(page, el("statistics viewer"), 900, 400));
    await session.step(15, "Then the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(18, "Then the \"view\" area of statistics viewer should not contain the color \"#9900FF\"", () => areaNotColor(page, "view", el("statistics viewer"), "#9900FF"));
    await session.step(19, "When user sets \"backColor\" property of statistics viewer to \"#9900FF\"", () => setProperty(page, "backColor", el("statistics viewer"), "#9900FF"));
    await session.step(20, "Then the \"view\" area of statistics viewer should contain the color \"#9900FF\"", () => areaColor(page, "view", el("statistics viewer"), "#9900FF"));
    await session.step(21, "And the \"nulls of HEIGHT\" reading of statistics viewer should be \"128\"", () => readingReads(page, "nulls of HEIGHT", el("statistics viewer"), "128"));
    await session.step(22, "And the \"columns shown\" reading of statistics viewer should be 11", () => readingIs(page, "columns shown", el("statistics viewer"), 11));
    await session.step(23, "When user sets \"backColor\" property of statistics viewer to \"#FFFFFF\"", () => setProperty(page, "backColor", el("statistics viewer"), "#FFFFFF"));
    await session.step(24, "Then the \"view\" area of statistics viewer should not contain the color \"#9900FF\"", () => areaNotColor(page, "view", el("statistics viewer"), "#9900FF"));
    await session.step(25, "And no errors should have been logged", () => noErrors(page));
  });
});
