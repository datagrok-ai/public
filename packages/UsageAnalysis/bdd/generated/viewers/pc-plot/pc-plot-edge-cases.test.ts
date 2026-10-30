/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/pc-plot/pc-plot-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.pc-plot, core#967, core#1203]
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
import {axesShouldBe} from '../../../bindings/pc-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {hoverOver} from '@datagrok-libraries/bdd/bindings/common/steps';
import {deleteSelected, filterPassesAll, filterPassesFewer, rowCount, selectAllRows, setCell} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, applyNamedLayout, hoverArea, noErrors, painted, pointerAway, propertyShouldBe, readingDoesNotRead, readingIs, readingReads, saveNamedLayout, setProperties, setProperty, showsRows, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {dragRangeHandle} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("PC plot edge cases — a layout with line width 0, an unparsable font, one column, an empty source", () => {
  const session = feature(test, "features/viewers/pc-plot/pc-plot-edge-cases.feature", import.meta.url);
  test("A layout saved with Line Width 0 and density on is applied and hovered without an error (core#967)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(22, "When user sets properties of pc plot viewer:", () => setProperties(page, el("pc plot viewer"), [["Show Density","true"],["Line Width","0"]]), [["Show Density","true"],["Line Width","0"]]);
    await session.step(25, "And user saves the layout of the current table view as \"pc-width-0\"", () => saveNamedLayout(page, "pc-width-0"));
    await session.step(26, "And user applies the layout \"pc-width-0\" to the current table view", () => applyNamedLayout(page, "pc-width-0"));
    await session.step(27, "Then \"Line Width\" property of pc plot viewer should be \"0\"", () => propertyShouldBe(page, "Line Width", el("pc plot viewer"), "0"));
    await session.step(28, "And the \"density style\" reading of pc plot viewer should be \"circles\"", () => readingReads(page, "density style", el("pc plot viewer"), "circles"));
    await session.step(29, "And pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
    await session.step(30, "When user hovers over the \"view\" area of pc plot viewer", () => hoverArea(page, "view", el("pc plot viewer")));
    await session.step(31, "And user moves the pointer away from pc plot viewer", () => pointerAway(page, el("pc plot viewer")));
    await session.step(32, "Then the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(33, "And no errors should have been logged", () => noErrors(page));
  });
  test("An Axis Font the browser cannot parse leaves the plot painting without an error (core#967)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(36, "When user sets \"Axis Font\" property of pc plot viewer to \"bold\"", () => setProperty(page, "Axis Font", el("pc plot viewer"), "bold"));
    await session.step(37, "And user takes a snapshot of pc plot viewer", () => takeSnapshot(page, el("pc plot viewer")));
    await session.step(38, "And user hovers over the \"view\" area of pc plot viewer", () => hoverArea(page, "view", el("pc plot viewer")));
    await session.step(39, "Then pc plot viewer should be painted", () => painted(page, el("pc plot viewer")));
    await session.step(40, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("With a single column the plot asks for at least two, and a second column brings it back", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(44, "When user sets \"Column Names\" property of pc plot viewer to \"AGE\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE"));
    await session.step(45, "Then the \"error\" reading of pc plot viewer should be \"please select at least 2 columns\"", () => readingReads(page, "error", el("pc plot viewer"), "please select at least 2 columns"));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
    await session.step(47, "When user sets \"Column Names\" property of pc plot viewer to \"AGE, WEIGHT\"", () => setProperty(page, "Column Names", el("pc plot viewer"), "AGE, WEIGHT"));
    await session.step(48, "Then the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(49, "And the axes of pc plot viewer should be \"AGE, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, WEIGHT"));
    await session.step(50, "And the \"lines drawn\" reading of pc plot viewer should be 1000", () => readingIs(page, "lines drawn", el("pc plot viewer"), 1000));
    await session.step(51, "And no errors should have been logged", () => noErrors(page));
  });
  test("Once every row of the table is deleted the plot draws no line and logs no error", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(54, "When user selects all rows", () => selectAllRows(page));
    await session.step(55, "And user deletes the selected rows", () => deleteSelected(page));
    await session.step(56, "Then the table should have 0 rows", () => rowCount(page, 0));
    await session.step(57, "And pc plot viewer should show 0 rows", () => showsRows(page, el("pc plot viewer"), 0));
    await session.step(58, "And the \"lines drawn\" reading of pc plot viewer should be 0", () => readingIs(page, "lines drawn", el("pc plot viewer"), 0));
    await session.step(59, "And no errors should have been logged", () => noErrors(page));
  });
  test("A transformation turned on releases the filter of the axis sliders it replaces (core#967)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot", "@realizes:core#967"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await session.step(63, "When user hovers over pc plot viewer", () => hoverOver(page, el("pc plot viewer")));
    await session.step(64, "And user drags the max handle of the \"AGE\" range slider of pc plot viewer by 120 pixels", () => dragRangeHandle(page, "max", "AGE", el("pc plot viewer"), 120));
    await session.step(65, "Then fewer than 1000 rows should pass the filter", () => filterPassesFewer(page, 1000));
    await session.step(66, "When user sets \"Transformation\" property of pc plot viewer to '[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"AGE\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]'", () => setProperty(page, "Transformation", el("pc plot viewer"), "[{\"#type\":\"GroupAggregation\",\"aggType\":\"key\",\"colName\":\"SEX\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"AGE\"},{\"#type\":\"GroupAggregation\",\"aggType\":\"avg\",\"colName\":\"WEIGHT\"}]"));
    await session.step(67, "Then the \"transformed\" reading of pc plot viewer should be \"true\"", () => readingReads(page, "transformed", el("pc plot viewer"), "true"));
    await session.step(68, "And all rows should pass the filter", () => filterPassesAll(page));
    await session.step(69, "And no errors should have been logged", () => noErrors(page));
  });
  test("Log Columns on a column holding zero says the axis cannot be logarithmic (core#1203)", {tag: ["@viewers", "@edge-cases", "@realizes:viewers.pc-plot", "@known-failure:76", "@realizes:core#1203"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(13, "Given user is logged in", () => loggedIn(page));
    await session.step(14, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(15, "And user adds a pc plot viewer with:", () => addViewerWith(page, "pc plot", [["Column Names","AGE, HEIGHT, WEIGHT"]]), [["Column Names","AGE, HEIGHT, WEIGHT"]]);
    await session.step(17, "Then the axes of pc plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => axesShouldBe(page, el("pc plot viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(18, "And pc plot viewer should show 1000 rows", () => showsRows(page, el("pc plot viewer"), 1000));
    await session.step(19, "And the \"error\" reading of pc plot viewer should be \"\"", () => readingReads(page, "error", el("pc plot viewer"), ""));
    await knownFailure(async () => {
      await session.step(73, "When user sets \"WEIGHT\" column in row 1 to \"0\"", () => setCell(page, "WEIGHT", 1, "0"));
      await session.step(74, "And user sets \"Log Columns\" property of pc plot viewer to \"WEIGHT\"", () => setProperty(page, "Log Columns", el("pc plot viewer"), "WEIGHT"));
      await session.step(75, "Then \"Log Columns\" property of pc plot viewer should be \"WEIGHT\"", () => propertyShouldBe(page, "Log Columns", el("pc plot viewer"), "WEIGHT"));
      await session.step(76, "And the \"error\" reading of pc plot viewer should not be \"\"", () => readingDoesNotRead(page, "error", el("pc plot viewer"), ""));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    }, 76);
  });
});
