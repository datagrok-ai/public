/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-p-value-filter.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot, core#1774]
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
import {addViewerWith, hasArea, noErrors, readingNotAsRemembered, rememberReading, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot p-value under the viewer's own filter", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-p-value-filter.feature", import.meta.url);
  test("The p-value follows the viewer's Filter", {tag: ["@viewers", "@realizes:viewers.box-plot", "@realizes:core#1774"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(9, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"]]), [["Value","AGE"],["Category 1","RACE"]]);
    await session.step(12, "Then box plot viewer should have a \"p value\" area", () => hasArea(page, el("box plot viewer"), "p value"));
    await session.step(13, "When user remembers the \"p value\" reading of box plot viewer", () => rememberReading(page, "p value", el("box plot viewer")));
    await session.step(14, "And user sets \"Filter\" property of box plot viewer to \"${SEX} == \\\"F\\\"\"", () => setProperty(page, "Filter", el("box plot viewer"), "${SEX} == \"F\""));
    await session.step(15, "Then the \"p value\" reading of box plot viewer should not be as remembered", () => readingNotAsRemembered(page, "p value", el("box plot viewer")));
    await session.step(16, "When user sets \"Show Group Comparison\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Group Comparison", el("box plot viewer"), "true"));
    await session.step(17, "Then box plot viewer should have a \"group comparison\" area", () => hasArea(page, el("box plot viewer"), "group comparison"));
    await session.step(18, "When user remembers the \"p value\" reading of box plot viewer", () => rememberReading(page, "p value", el("box plot viewer")));
    await session.step(19, "And user sets \"Filter\" property of box plot viewer to \"${SEX} == \\\"M\\\"\"", () => setProperty(page, "Filter", el("box plot viewer"), "${SEX} == \"M\""));
    await session.step(20, "Then the \"p value\" reading of box plot viewer should not be as remembered", () => readingNotAsRemembered(page, "p value", el("box plot viewer")));
    await session.step(21, "When user remembers the \"p value\" reading of box plot viewer", () => rememberReading(page, "p value", el("box plot viewer")));
    await session.step(22, "And user sets \"Filter\" property of box plot viewer to \"\"", () => setProperty(page, "Filter", el("box plot viewer"), ""));
    await session.step(23, "Then the \"p value\" reading of box plot viewer should not be as remembered", () => readingNotAsRemembered(page, "p value", el("box plot viewer")));
    await session.step(24, "And no errors should have been logged", () => noErrors(page));
  });
});
