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
import {rangeInputReads} from '../../../bindings/histogram.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterIsExactly, filterPasses, filterPassesAll, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, enterIntoArea, hasNoArea, loadLayout, noErrors, readingIs, readingReads, resizeTo, saveLayout, setProperty, showsRows} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
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
  test("On a date column a maximum below the minimum empties the filter and is reported, and a bound past the dates is clamped", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(41, "When user sets \"Value\" property of histogram viewer to \"STARTED\"", () => setProperty(page, "Value", el("histogram viewer"), "STARTED"));
    await session.step(42, "And user enters \"1/1/1991\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "1/1/1991", "range min input", el("histogram viewer")));
    await session.step(43, "Then 483 rows should pass the filter", () => filterPasses(page, 483));
    await session.step(44, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(45, "When user enters \"1/1/1990\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "1/1/1990", "range max input", el("histogram viewer")));
    await session.step(46, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(47, "And the \"error\" reading of histogram viewer should be \"max should be greater than min\"", () => readingReads(page, "error", el("histogram viewer"), "max should be greater than min"));
    await session.step(48, "When user enters \"1/1/1992\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "1/1/1992", "range max input", el("histogram viewer")));
    await session.step(49, "Then 483 rows should pass the filter", () => filterPasses(page, 483));
    await session.step(50, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(51, "And the range max input of histogram viewer should read \"1/1/1992\"", () => rangeInputReads(page, "max", "1/1/1992"));
    await session.step(52, "And no errors should have been logged", () => noErrors(page));
  });
  test("On a date column a range normalized to the filter logs no error when every row is filtered out", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(56, "When user sets \"Value\" property of histogram viewer to \"STARTED\"", () => setProperty(page, "Value", el("histogram viewer"), "STARTED"));
    await session.step(57, "And user sets \"Normalize To Filter\" property of histogram viewer to \"true\"", () => setProperty(page, "Normalize To Filter", el("histogram viewer"), "true"));
    await session.step(58, "And user enters \"1/1/1991\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "1/1/1991", "range min input", el("histogram viewer")));
    await session.step(59, "Then 483 rows should pass the filter", () => filterPasses(page, 483));
    await session.step(60, "When user adds a categorical filter on \"SEX\" keeping \"\"", () => addCategoricalFilter(page, "SEX", ""));
    await session.step(61, "Then 0 rows should pass the filter", () => filterPasses(page, 0));
    await session.step(62, "And histogram viewer should show 0 rows", () => showsRows(page, el("histogram viewer"), 0));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
  test("Text typed into the range min input is not taken and the filter stays as it was", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(67, "When user enters \"abc\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "abc", "range min input", el("histogram viewer")));
    await session.step(68, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(69, "And histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(70, "And no errors should have been logged", () => noErrors(page));
  });
  test("Switching Value while a range filters releases the range", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(74, "When user enters \"30\" into the \"range min input\" area of histogram viewer", () => enterIntoArea(page, "30", "range min input", el("histogram viewer")));
    await session.step(75, "And user enters \"50\" into the \"range max input\" area of histogram viewer", () => enterIntoArea(page, "50", "range max input", el("histogram viewer")));
    await session.step(76, "Then the filter should pass exactly the rows where \"AGE\" is between 30 and 50", () => filterIsExactly(page, "AGE", 30, 50));
    await session.step(77, "When user sets \"Value\" property of histogram viewer to \"WEIGHT\"", () => setProperty(page, "Value", el("histogram viewer"), "WEIGHT"));
    await session.step(78, "Then all rows should pass the filter", () => filterPassesAll(page));
    await session.step(79, "And no errors should have been logged", () => noErrors(page));
  });
  test("A layout whose Value column the table no longer has leaves the histogram empty with a message", {tag: ["@viewers", "@realizes:viewers.histogram", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(11, "And user adds a histogram viewer with:", () => addViewerWith(page, "histogram", [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]), [["Value","AGE"],["Show Range Inputs","true"],["Filtering Enabled","true"]]);
    await session.step(15, "And user resizes histogram viewer to 500 by 400", () => resizeTo(page, el("histogram viewer"), 500, 400));
    await session.step(16, "Then histogram viewer should show 1000 rows", () => showsRows(page, el("histogram viewer"), 1000));
    await session.step(17, "And the \"error\" reading of histogram viewer should be \"\"", () => readingReads(page, "error", el("histogram viewer"), ""));
    await session.step(83, "When user sets \"Value\" property of histogram viewer to \"HEIGHT\"", () => setProperty(page, "Value", el("histogram viewer"), "HEIGHT"));
    await session.step(84, "And user saves the layout of the current table view", () => saveLayout(page));
    await session.step(85, "And user removes \"HEIGHT\" column", () => removeColumn(page, "HEIGHT"));
    await session.step(86, "And user loads the saved layout", () => loadLayout(page));
    await session.step(87, "Then histogram viewer should contain the text \"does not exist\"", () => shouldContainText(page, el("histogram viewer"), "does not exist"));
    await session.step(88, "And histogram viewer should not have a \"bin 1\" area", () => hasNoArea(page, el("histogram viewer"), "bin 1"));
    await session.step(89, "And no errors should have been logged", () => noErrors(page));
  });
});
