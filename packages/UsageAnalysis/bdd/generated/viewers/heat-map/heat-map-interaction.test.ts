/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/heat-map/heat-map-interaction.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.heat-map]
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
import {pressKey} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, noneSelected, rowsRangeSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, doubleClickArea, dragAreaBy, dragSelectionOverArea, hasArea, hasNoArea, noErrors, painted, readingBetween, readingHigher, readingIs, readingLower, readingReads, repainted, setProperty, someHighlight} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Heat map scroll sliders, rectangle selection and the Table property", () => {
  const session = feature(test, "features/viewers/heat-map/heat-map-interaction.feature", import.meta.url);
  test("Heat map scroll sliders, rectangle selection and the Table property", {tag: ["@journey", "@viewers", "@realizes:viewers.heat-map"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 4, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens cars dataset", () => openDataset(page, ds("cars")));
    await session.step(17, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(18, "And user adds a heat map viewer", () => addViewer(page, "heat map"));
    await session.step(19, "Then the \"is heatmap\" reading of heat map viewer should be \"true\"", () => readingReads(page, "is heatmap", el("heat map viewer"), "true"));
    await session.step(20, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
    await session.step(21, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
    await session.step(22, "And heat map viewer should be painted", () => painted(page, el("heat map viewer")));
    await run.scenario("Dragging the vertical slider's handle narrows the rows, and a double-click resets them", async () => {
      await session.step(25, "When user drags the \"y scroll max handle\" area of heat map viewer by 300 pixels to the up", () => dragAreaBy(page, "y scroll max handle", el("heat map viewer"), 300, "up"));
      await session.step(26, "Then the \"y scroll span\" reading of heat map viewer should be between 0.5 and 0.9", () => readingBetween(page, "y scroll span", el("heat map viewer"), 0.5, 0.9));
      await session.step(27, "And the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(28, "And the \"row height\" reading of heat map viewer should be higher than before", () => readingHigher(page, "row height", el("heat map viewer")));
      await session.step(29, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(30, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(31, "When user double-clicks on the \"y scroll slider\" area of heat map viewer", () => doubleClickArea(page, "y scroll slider", el("heat map viewer")));
      await session.step(32, "Then the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(33, "And the \"row height\" reading of heat map viewer should be lower than before", () => readingLower(page, "row height", el("heat map viewer")));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Dragging the horizontal slider's handle narrows the columns, and a double-click resets them", async () => {
      await session.step(37, "When user drags the \"x scroll min handle\" area of heat map viewer by 200 pixels to the right", () => dragAreaBy(page, "x scroll min handle", el("heat map viewer"), 200, "right"));
      await session.step(38, "Then the \"x scroll span\" reading of heat map viewer should be between 0.4 and 0.9", () => readingBetween(page, "x scroll span", el("heat map viewer"), 0.4, 0.9));
      await session.step(39, "And the \"y scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "y scroll span", el("heat map viewer"), 1));
      await session.step(40, "And the \"columns shown\" reading of heat map viewer should be lower than before", () => readingLower(page, "columns shown", el("heat map viewer")));
      await session.step(41, "And heat map viewer should not have a \"header USUBJID\" area", () => hasNoArea(page, el("heat map viewer"), "header USUBJID"));
      await session.step(42, "And heat map viewer should have a \"header SEVERITY\" area", () => hasArea(page, el("heat map viewer"), "header SEVERITY"));
      await session.step(43, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(44, "When user double-clicks on the \"x scroll slider\" area of heat map viewer", () => doubleClickArea(page, "x scroll slider", el("heat map viewer")));
      await session.step(45, "Then the \"x scroll span\" reading of heat map viewer should be 1", () => readingIs(page, "x scroll span", el("heat map viewer"), 1));
      await session.step(46, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(47, "And heat map viewer should have a \"header USUBJID\" area", () => hasArea(page, el("heat map viewer"), "header USUBJID"));
      await session.step(48, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A Shift-drag selects the rows it spans and Escape clears them", async () => {
      await session.step(51, "Given user clears the row selection", () => clearSelection(page));
      await session.step(52, "When user drags a selection box over the \"column AGE\" area of heat map viewer", () => dragSelectionOverArea(page, "column AGE", el("heat map viewer")));
      await session.step(53, "Then rows 100 to 892 should be selected", () => rowsRangeSelected(page, 100, 892));
      await session.step(54, "And heat map viewer should show a selection highlight", () => someHighlight(page, el("heat map viewer")));
      await session.step(55, "When user presses Escape", () => pressKey(page, "Escape"));
      await session.step(56, "Then no rows should be selected", () => noneSelected(page));
      await session.step(57, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Table property binds the map to another table and back", async () => {
      await session.step(60, "When user sets \"table\" property of heat map viewer to \"cars\"", () => setProperty(page, "table", el("heat map viewer"), "cars"));
      await session.step(61, "Then heat map viewer should be bound to table \"cars\"", () => boundTable(page, el("heat map viewer"), "cars"));
      await session.step(62, "And the \"rows\" reading of heat map viewer should be 30", () => readingIs(page, "rows", el("heat map viewer"), 30));
      await session.step(63, "And the \"columns shown\" reading of heat map viewer should be 17", () => readingIs(page, "columns shown", el("heat map viewer"), 17));
      await session.step(64, "And heat map viewer should have a \"cell 30 of model\" area", () => hasArea(page, el("heat map viewer"), "cell 30 of model"));
      await session.step(65, "And the \"text of cell 1 of model\" reading of heat map viewer should be \"alfaromeo\"", () => readingReads(page, "text of cell 1 of model", el("heat map viewer"), "alfaromeo"));
      await session.step(66, "And heat map viewer should not have a \"column AGE\" area", () => hasNoArea(page, el("heat map viewer"), "column AGE"));
      await session.step(67, "And heat map viewer should have repainted", () => repainted(page, el("heat map viewer")));
      await session.step(68, "When user sets \"table\" property of heat map viewer to \"demog-1000\"", () => setProperty(page, "table", el("heat map viewer"), "demog-1000"));
      await session.step(69, "Then heat map viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("heat map viewer"), "demog-1000"));
      await session.step(70, "And the \"rows\" reading of heat map viewer should be 1000", () => readingIs(page, "rows", el("heat map viewer"), 1000));
      await session.step(71, "And the \"columns shown\" reading of heat map viewer should be 11", () => readingIs(page, "columns shown", el("heat map viewer"), 11));
      await session.step(72, "And the \"column order\" reading of heat map viewer should be \"USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY\"", () => readingReads(page, "column order", el("heat map viewer"), "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"));
      await session.step(73, "And the \"row height\" reading of heat map viewer should be between 0 and 8", () => readingBetween(page, "row height", el("heat map viewer"), 0, 8));
      await session.step(74, "And heat map viewer should have a \"column AGE\" area", () => hasArea(page, el("heat map viewer"), "column AGE"));
      await session.step(75, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
