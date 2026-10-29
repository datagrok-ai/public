/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/form/form-edge-cases.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.form, core#1138, core#1192]
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
import {uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs, makeLastRowCurrent, valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, enterIntoArea, errorBalloonText, hasNoArea, noErrors, readingIs, reportsNoError, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Form viewer edge cases", () => {
  const session = feature(test, "features/viewers/form/form-edge-cases.feature", import.meta.url);
  test("The row arrows stop at the first and the last row", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(22, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
    await session.step(23, "Then row 2 should be current", () => currentRowIs(page, 2));
    await session.step(24, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
    await session.step(25, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(26, "When user clicks on the \"previous row\" area of form viewer", () => clickArea(page, "previous row", el("form viewer")));
    await session.step(27, "Then row 1 should be current", () => currentRowIs(page, 1));
    await session.step(28, "And the \"row\" reading of form viewer should be 1", () => readingIs(page, "row", el("form viewer"), 1));
    await session.step(29, "When user makes the last row current", () => makeLastRowCurrent(page));
    await session.step(30, "Then the \"row\" reading of form viewer should be 1000", () => readingIs(page, "row", el("form viewer"), 1000));
    await session.step(31, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
    await session.step(32, "Then row 1000 should be current", () => currentRowIs(page, 1000));
    await session.step(33, "And the \"row\" reading of form viewer should be 1000", () => readingIs(page, "row", el("form viewer"), 1000));
    await session.step(34, "And no errors should have been logged", () => noErrors(page));
  });
  test("A file that is not a form is refused with a message and the form keeps its fields [file=fixtures/cars-small.csv, name=cars-small.csv]", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@realizes:core#1138"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(38, "When user uploads \"fixtures/cars-small.csv\" through folder open icon in form viewer", () => uploadThrough(page, "fixtures/cars-small.csv", el("folder open icon in form viewer")));
    await session.step(39, "Then an error balloon containing \"cars-small.csv is not a valid form file\" should have been shown", () => errorBalloonText(page, "cars-small.csv is not a valid form file"));
    await session.step(40, "And the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("A file that is not a form is refused with a message and the form keeps its fields [file=fixtures/nx/spgi-legend-position-old.layout, name=spgi-legend-position-old.layout]", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@realizes:core#1138"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(38, "When user uploads \"fixtures/nx/spgi-legend-position-old.layout\" through folder open icon in form viewer", () => uploadThrough(page, "fixtures/nx/spgi-legend-position-old.layout", el("folder open icon in form viewer")));
    await session.step(39, "Then an error balloon containing \"spgi-legend-position-old.layout is not a valid form file\" should have been shown", () => errorBalloonText(page, "spgi-legend-position-old.layout is not a valid form file"));
    await session.step(40, "And the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
  });
  test("Text typed into a numeric field of an editable form leaves the cell as it was", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@known-failure:54", "@realizes:core#1192"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await knownFailure(async () => {
      await session.step(50, "When user clicks on the \"edit\" area of form viewer", () => clickArea(page, "edit", el("form viewer")));
      await session.step(51, "And user enters \"31\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "31", "field AGE", el("form viewer")));
      await session.step(52, "Then the value of \"AGE\" column in row 1 should be \"31\"", () => valueInRow(page, "AGE", 1, "31"));
      await session.step(53, "When user enters \"abc\" into the \"field AGE\" area of form viewer", () => enterIntoArea(page, "abc", "field AGE", el("form viewer")));
      await session.step(54, "Then the value of \"AGE\" column in row 1 should be \"31\"", () => valueInRow(page, "AGE", 1, "31"));
      await session.step(55, "And no errors should have been logged", () => noErrors(page));
    }, 54);
  });
  test("A column removed from the table takes its field off the form", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await session.step(59, "When user removes \"AGE\" column", () => removeColumn(page, "AGE"));
    await session.step(60, "Then the \"fields shown\" reading of form viewer should be 10", () => readingIs(page, "fields shown", el("form viewer"), 10));
    await session.step(61, "And form viewer should not have a \"field AGE\" area", () => hasNoArea(page, el("form viewer"), "field AGE"));
    await session.step(62, "And form viewer should report no error", () => reportsNoError(page, el("form viewer")));
    await session.step(63, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column removed from the table takes its field off a form built from chosen columns", {tag: ["@viewers", "@realizes:viewers.form", "@edge-cases", "@known-failure:70", "@ticket-pending"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(18, "Then the \"fields shown\" reading of form viewer should be 11", () => readingIs(page, "fields shown", el("form viewer"), 11));
    await knownFailure(async () => {
      await session.step(67, "Given user sets \"columnNames\" property of form viewer to \"AGE, HEIGHT, WEIGHT\"", () => setProperty(page, "columnNames", el("form viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(68, "Then the \"fields shown\" reading of form viewer should be 3", () => readingIs(page, "fields shown", el("form viewer"), 3));
      await session.step(69, "When user removes \"WEIGHT\" column", () => removeColumn(page, "WEIGHT"));
      await session.step(70, "Then the \"fields shown\" reading of form viewer should be 2", () => readingIs(page, "fields shown", el("form viewer"), 2));
      await session.step(71, "And form viewer should not have a \"field WEIGHT\" area", () => hasNoArea(page, el("form viewer"), "field WEIGHT"));
      await session.step(72, "And no errors should have been logged", () => noErrors(page));
    }, 70);
  });
});
