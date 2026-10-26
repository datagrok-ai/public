/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/form/form-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.form, public#3561, public#940, public#2644]
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
import {columnType, currentRowIs, makeRowCurrent} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {colorLinear} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, clickArea, doubleClickArea, noErrors, readingAsRemembered, readingDoesNotRead, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Form viewer regressions from closed issues", () => {
  const session = feature(test, "features/viewers/form/form-regressions.feature", import.meta.url);
  test("An empty integer cell shows an empty field on the form viewer (public#3561)", {tag: ["@viewers", "@realizes:viewers.form", "@realizes:public#3561"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "Given user opens a table \"ints\" with:", () => openTableOf(page, "ints", [["id","n"],["1","5"],["2",""],["3","7"]]), [["id","n"],["1","5"],["2",""],["3","7"]]);
    await session.step(21, "Then \"n\" column should have type \"int\"", () => columnType(page, "n", "int"));
    await session.step(22, "When user adds a form viewer", () => addViewer(page, "form"));
    await session.step(23, "And user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(24, "Then the \"n\" reading of form viewer should be \"5\"", () => readingReads(page, "n", el("form viewer"), "5"));
    await session.step(25, "When user makes row 2 current", () => makeRowCurrent(page, 2));
    await session.step(26, "Then the \"row\" reading of form viewer should be 2", () => readingIs(page, "row", el("form viewer"), 2));
    await session.step(27, "And the \"n\" reading of form viewer should be \"\"", () => readingReads(page, "n", el("form viewer"), ""));
    await session.step(28, "And no errors should have been logged", () => noErrors(page));
  });
  test("A field keeps its own row's colour after the grid is sorted (public#940)", {tag: ["@viewers", "@realizes:viewers.form", "@realizes:public#940"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(33, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(34, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(35, "When user colors \"AGE\" column linearly from \"#FF0000\" to \"#0000FF\"", () => colorLinear(page, "AGE", "#FF0000", "#0000FF"));
    await session.step(36, "And user makes row 5 current", () => makeRowCurrent(page, 5));
    await session.step(37, "Then the \"background of AGE\" reading of form viewer should not be \"\"", () => readingDoesNotRead(page, "background of AGE", el("form viewer"), ""));
    await session.step(38, "When user remembers the \"background of AGE\" reading of form viewer", () => rememberReading(page, "background of AGE", el("form viewer")));
    await session.step(39, "And user makes row 7 current", () => makeRowCurrent(page, 7));
    await session.step(40, "Then the \"background of AGE\" reading of form viewer should not be as remembered", () => readingNotAsRemembered(page, "background of AGE", el("form viewer")));
    await session.step(41, "When user double-clicks on the \"header AGE\" area of grid", () => doubleClickArea(page, "header AGE", el("grid")));
    await session.step(42, "Then the \"sort column\" reading of grid should be \"AGE\"", () => readingReads(page, "sort column", el("grid"), "AGE"));
    await session.step(43, "When user makes row 5 current", () => makeRowCurrent(page, 5));
    await session.step(44, "Then the \"row\" reading of form viewer should be 5", () => readingIs(page, "row", el("form viewer"), 5));
    await session.step(45, "And the \"background of AGE\" reading of form viewer should be as remembered", () => readingAsRemembered(page, "background of AGE", el("form viewer")));
    await session.step(46, "And no errors should have been logged", () => noErrors(page));
  });
  test("After the grid is sorted the next-row arrow takes the next row of the grid (public#2644)", {tag: ["@viewers", "@realizes:viewers.form", "@realizes:public#2644"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(51, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(52, "And user adds a form viewer", () => addViewer(page, "form"));
    await session.step(53, "When user double-clicks on the \"header HEIGHT\" area of grid", () => doubleClickArea(page, "header HEIGHT", el("grid")));
    await session.step(54, "Then the \"sort column\" reading of grid should be \"HEIGHT\"", () => readingReads(page, "sort column", el("grid"), "HEIGHT"));
    await session.step(55, "And the \"sort direction\" reading of grid should be \"descending\"", () => readingReads(page, "sort direction", el("grid"), "descending"));
    await session.step(56, "When user makes row 1 current", () => makeRowCurrent(page, 1));
    await session.step(57, "Then the \"HEIGHT\" reading of form viewer should be \"174.705\"", () => readingReads(page, "HEIGHT", el("form viewer"), "174.705"));
    await session.step(58, "When user clicks on the \"next row\" area of form viewer", () => clickArea(page, "next row", el("form viewer")));
    await session.step(59, "Then row 808 should be current", () => currentRowIs(page, 808));
    await session.step(60, "And the \"row\" reading of form viewer should be 808", () => readingIs(page, "row", el("form viewer"), 808));
    await session.step(61, "And the \"HEIGHT\" reading of form viewer should be \"174.652\"", () => readingReads(page, "HEIGHT", el("form viewer"), "174.652"));
    await session.step(62, "And no errors should have been logged", () => noErrors(page));
  });
});
