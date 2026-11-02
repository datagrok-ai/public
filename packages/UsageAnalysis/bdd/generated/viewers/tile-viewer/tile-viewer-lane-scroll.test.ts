/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-lane-scroll.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer, core#1217]
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
import {addTileViewerWith, dragCardIntoLane} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {valueInRow} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {openDataset, openTableOf} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, noErrors, noViewerError, pointerAway, readingAsRemembered, readingHigher, readingIs, readingReads, rememberReading, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer — a scrolled lane keeps its place while viewers dock and close beside it", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-lane-scroll.feature", import.meta.url);
  test("The lane keeps its scroll when the docked histogram closes", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@known-failure:37"]}, async ({browser}) => {
    const page = await session.page(browser);
    await knownFailure(async () => {
      await session.step(19, "Given user is logged in", () => loggedIn(page));
      await session.step(20, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(21, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
      await session.step(23, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
      await session.step(24, "And the \"scroll of lane Caucasian\" reading of tile viewer should be 0", () => readingIs(page, "scroll of lane Caucasian", el("tile viewer"), 0));
      await session.step(25, "When user scrolls the mouse wheel down over the \"lane content Caucasian\" area of tile viewer", () => wheelOverArea(page, "down", "lane content Caucasian", el("tile viewer")));
      await session.step(26, "Then the \"scroll of lane Caucasian\" reading of tile viewer should be higher than before", () => readingHigher(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(27, "When user remembers the \"scroll of lane Caucasian\" reading of tile viewer", () => rememberReading(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(28, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(29, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
      await session.step(30, "And no viewer of the current view should report an error", () => noViewerError(page));
      await session.step(31, "And the \"scroll of lane Caucasian\" reading of tile viewer should be as remembered", () => readingAsRemembered(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(32, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(33, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(34, "And the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(35, "And the \"lanes\" reading of tile viewer should be 4", () => readingIs(page, "lanes", el("tile viewer"), 4));
      await session.step(36, "And no viewer of the current view should report an error", () => noViewerError(page));
      await session.step(37, "And the \"scroll of lane Caucasian\" reading of tile viewer should be as remembered", () => readingAsRemembered(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(38, "And no errors should have been logged", () => noErrors(page));
    }, 37);
  });
  test("A card dropped into the (empty) lane gets an empty value", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@known-failure:51", "@realizes:core#1217"]}, async ({browser}) => {
    const page = await session.page(browser);
    await knownFailure(async () => {
      await session.step(42, "Given user is logged in", () => loggedIn(page));
      await session.step(43, "And user opens a table \"with-empty\" with:", () => openTableOf(page, "with-empty", [["site","arm"],["A","x"],["A",""],["B","y"]]), [["site","arm"],["A","x"],["A",""],["B","y"]]);
      await session.step(48, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","arm"]]), [["Lanes Column Name","arm"]]);
      await session.step(50, "Then the \"lane of row 1\" reading of tile viewer should be \"x\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "x"));
      await session.step(51, "And tile viewer should have a \"lane (empty)\" area", () => hasArea(page, el("tile viewer"), "lane (empty)"));
      await session.step(52, "When user moves the pointer away from tile viewer", () => pointerAway(page, el("tile viewer")));
      await session.step(53, "And user drags the card of row 1 of tile viewer into lane \"(empty)\"", () => dragCardIntoLane(page, 1, el("tile viewer"), "(empty)"));
      await session.step(54, "Then the value of \"arm\" column in row 1 should be \"\"", () => valueInRow(page, "arm", 1, ""));
      await session.step(55, "And the \"lane of row 1\" reading of tile viewer should be \"(empty)\"", () => readingReads(page, "lane of row 1", el("tile viewer"), "(empty)"));
    }, 51);
  });
});
