/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/histogram/histogram-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.histogram, core#998]
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
import {addCategoricalFilter, filterIsExactly, filterPasses} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, enterIntoArea, hasNoArea, noErrors, readingIs, readingReads, resizeTo, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Histogram edge cases", () => {
  const session = feature(test, "features/viewers/histogram/histogram-edge-cases.feature", import.meta.url);
  test("Bins set to 0 is reported instead of drawn", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases", "@realizes:core#998"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(21, "When user sets \"Bins\" property of histogram viewer to \"0\"", () => setProperty(page, "Bins", el("histogram viewer"), "0"));
    await session.step(22, "Then the \"error\" reading of histogram viewer should be \"bins should be a positive integer\"", () => readingReads(page, "error", el("histogram viewer"), "bins should be a positive integer"));
    await session.step(23, "And histogram viewer should not have a \"bin 1\" area", () => hasNoArea(page, el("histogram viewer"), "bin 1"));
    await session.step(24, "When user sets \"Bins\" property of histogram viewer to \"20\"", () => setProperty(page, "Bins", el("histogram viewer"), "20"));
    await session.step(25, "Then the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(26, "And the \"bins shown\" reading of histogram viewer should be 20", () => readingIs(page, "bins shown", el("histogram viewer"), 20));
    await session.step(27, "And no errors should have been logged", () => noErrors(page));
  });
  test("A range normalized to the filter logs no error when every row is filtered out", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases", "@realizes:core#998"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(31, "When user sets \"Normalize To Filter\" property of histogram viewer to \"true\"", () => setProperty(page, "Normalize To Filter", el("histogram viewer"), "true"));
    await session.step(32, "And user enters \"60\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "60", "range min input", el("histogram viewer")));
    await session.step(33, "Then the filter should pass exactly the rows where \"AGE\" is between 60 and 89", () => filterIsExactly(page, "AGE", 60, 89));
    await session.step(34, "When user adds a categorical filter on \"SEX\" keeping \"\"", () => addCategoricalFilter(page, "SEX", ""));
    await session.step(35, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(36, "And histogram viewer should show 0 rows", () => showsRows(page, el("histogram viewer"), 0));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
