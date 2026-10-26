/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer, public#3855, public#1593]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/file-formats.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/queries.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {addTileViewer} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe, shouldContainText} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {hoverArea, noErrors, pickFromAreaContextMenu, pointerAway, readingReads} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {scrollGridTo} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer regressions from closed issues", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-regressions.feature", import.meta.url);
  test("A date format picked for the column reaches the card (public#3855)", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@realizes:public#3855"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(17, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(18, "Then the \"STARTED of row 1\" reading of tile viewer should be \"8/2/1990\"", () => readingReads(page, "STARTED of row 1", el("tile viewer"), "8/2/1990"));
    await session.step(19, "When user scrolls the grid to the \"STARTED\" column", () => scrollGridTo(page, "STARTED"));
    await session.step(20, "And user picks \"Format > dd.MM.yyyy (02.08.1990)\" from the context menu of the \"header STARTED\" area of grid", () => pickFromAreaContextMenu(page, "Format > dd.MM.yyyy (02.08.1990)", "header STARTED", el("grid")));
    await session.step(21, "Then the \"text of cell 1 of STARTED\" reading of grid should be \"02.08.1990\"", () => readingReads(page, "text of cell 1 of STARTED", el("grid"), "02.08.1990"));
    await session.step(22, "And the \"STARTED of row 1\" reading of tile viewer should be \"02.08.1990\"", () => readingReads(page, "STARTED of row 1", el("tile viewer"), "02.08.1990"));
    await session.step(23, "And no errors should have been logged", () => noErrors(page));
  });
  test("A column name the card cuts off shows whole in a tooltip over its caption (public#1593)", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@realizes:public#1593"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(28, "Given user opens a table \"long-names\" with:", () => openTableOf(page, "long-names", [["id","measured_concentration_of_compound_in_plasma"],["1","5"],["2","6"]]), [["id","measured_concentration_of_compound_in_plasma"],["1","5"],["2","6"]]);
    await session.step(32, "And user adds a tile viewer", () => addTileViewer(page));
    await session.step(33, "When user hovers over the \"label id of row 1\" area of tile viewer", () => hoverArea(page, "label id of row 1", el("tile viewer")));
    await session.step(34, "Then tooltip should be hidden", () => shouldBe(page, el("tooltip"), "hidden"));
    await session.step(35, "When user hovers over the \"label measured_concentration_of_compound_in_plasma of row 1\" area of tile viewer", () => hoverArea(page, "label measured_concentration_of_compound_in_plasma of row 1", el("tile viewer")));
    await session.step(36, "Then tooltip should be visible", () => shouldBe(page, el("tooltip"), "visible"));
    await session.step(37, "And tooltip should contain text \"measured_concentration_of_compound_in_plasma\"", () => shouldContainText(page, el("tooltip"), "measured_concentration_of_compound_in_plasma"));
    await session.step(38, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
    await session.step(39, "Then no errors should have been logged", () => noErrors(page));
  });
});
