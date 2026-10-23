/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/tile-viewer/tile-viewer-lane-scroll.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.tile-viewer]
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
import {addTileViewerWith} from '../../../bindings/tile-viewer.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, noErrors, noViewerError, readingAsRemembered, readingHigher, readingIs, rememberReading, viewerCount, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, knownFailure} from '@datagrok-libraries/bdd/runtime';

test.describe("Tile viewer — a scrolled lane keeps its place while viewers dock and close beside it", () => {
  const session = feature(test, "features/viewers/tile-viewer/tile-viewer-lane-scroll.feature", import.meta.url);
  test("The lane keeps its scroll when the docked histogram closes", {tag: ["@viewers", "@realizes:viewers.tile-viewer", "@known-failure:33"]}, async ({browser}) => {
    const page = await session.page(browser);
    await knownFailure(async () => {
      await session.step(15, "Given user is logged in", () => loggedIn(page));
      await session.step(16, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
      await session.step(17, "And user adds a tile viewer with:", () => addTileViewerWith(page, [["Lanes Column Name","RACE"]]), [["Lanes Column Name","RACE"]]);
      await session.step(19, "Then tile viewer should be visible", () => shouldBe(page, el("tile viewer"), "visible"));
      await session.step(20, "And the \"scroll of lane Caucasian\" reading of tile viewer should be 0", () => readingIs(page, "scroll of lane Caucasian", el("tile viewer"), 0));
      await session.step(21, "When user scrolls the mouse wheel down over the \"lane content Caucasian\" area of tile viewer", () => wheelOverArea(page, "down", "lane content Caucasian", el("tile viewer")));
      await session.step(22, "Then the \"scroll of lane Caucasian\" reading of tile viewer should be higher than before", () => readingHigher(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(23, "When user remembers the \"scroll of lane Caucasian\" reading of tile viewer", () => rememberReading(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(24, "And user adds a histogram viewer", () => addViewer(page, "histogram"));
      await session.step(25, "Then histogram viewer should be visible", () => shouldBe(page, el("histogram viewer"), "visible"));
      await session.step(26, "And no viewer of the current view should report an error", () => noViewerError(page));
      await session.step(27, "And the \"scroll of lane Caucasian\" reading of tile viewer should be as remembered", () => readingAsRemembered(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(28, "When user clicks on close icon of histogram viewer", () => clickOn(page, el("close icon of histogram viewer")));
      await session.step(29, "Then histogram viewer should be absent", () => shouldBe(page, el("histogram viewer"), "absent"));
      await session.step(30, "And the open tableview should have 0 histogram viewers", () => viewerCount(page, 0, "histogram"));
      await session.step(31, "And the \"lanes\" reading of tile viewer should be 4", () => readingIs(page, "lanes", el("tile viewer"), 4));
      await session.step(32, "And no viewer of the current view should report an error", () => noViewerError(page));
      await session.step(33, "And the \"scroll of lane Caucasian\" reading of tile viewer should be as remembered", () => readingAsRemembered(page, "scroll of lane Caucasian", el("tile viewer")));
      await session.step(34, "And no errors should have been logged", () => noErrors(page));
    }, 33);
  });
});
