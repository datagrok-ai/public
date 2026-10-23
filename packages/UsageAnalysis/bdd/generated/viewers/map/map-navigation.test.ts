/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/map/map-navigation.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.map-viewer]
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
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {clearSelection, noneSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, dragAreaBy, noErrors, readingAsRemembered, readingHigher, readingHigherThanRemembered, readingIs, readingLower, readingNotAsRemembered, rememberReading, wheelOverArea} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Map viewer — the mouse wheel zooms, a drag pans", () => {
  const session = feature(test, "features/viewers/map/map-navigation.feature", import.meta.url);
  test("Map viewer — the mouse wheel zooms, a drag pans", {tag: ["@journey", "@viewers", "@realizes:viewers.map-viewer"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens earthquakes dataset", () => openDataset(page, ds("earthquakes")));
    await session.step(13, "And user adds a map viewer", () => addViewer(page, "map"));
    await session.step(14, "Then map viewer should be visible", () => shouldBe(page, el("map viewer"), "visible"));
    await session.step(15, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
    await run.scenario("The mouse wheel zooms the map in, and back out to the zoom it had", async () => {
      await session.step(18, "Given user remembers the \"zoom\" reading of map viewer", () => rememberReading(page, "zoom", el("map viewer")));
      await session.step(19, "When user scrolls the mouse wheel up over the \"view\" area of map viewer", () => wheelOverArea(page, "up", "view", el("map viewer")));
      await session.step(20, "Then the \"zoom\" reading of map viewer should be higher than before", () => readingHigher(page, "zoom", el("map viewer")));
      await session.step(21, "And the \"zoom\" reading of map viewer should be higher than remembered", () => readingHigherThanRemembered(page, "zoom", el("map viewer")));
      await session.step(22, "When user scrolls the mouse wheel down over the \"view\" area of map viewer", () => wheelOverArea(page, "down", "view", el("map viewer")));
      await session.step(23, "Then the \"zoom\" reading of map viewer should be lower than before", () => readingLower(page, "zoom", el("map viewer")));
      await session.step(24, "And the \"zoom\" reading of map viewer should be as remembered", () => readingAsRemembered(page, "zoom", el("map viewer")));
      await session.step(25, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
      await session.step(26, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A drag on the map pans it and keeps the zoom and the selection", async () => {
      await session.step(29, "Given user clears the row selection", () => clearSelection(page));
      await session.step(30, "And user remembers the \"zoom\" reading of map viewer", () => rememberReading(page, "zoom", el("map viewer")));
      await session.step(31, "And user remembers the \"centre\" reading of map viewer", () => rememberReading(page, "centre", el("map viewer")));
      await session.step(32, "When user drags the \"view\" area of map viewer by 150 pixels to the left", () => dragAreaBy(page, "view", el("map viewer"), 150, "left"));
      await session.step(33, "Then the \"centre\" reading of map viewer should not be as remembered", () => readingNotAsRemembered(page, "centre", el("map viewer")));
      await session.step(34, "And the \"zoom\" reading of map viewer should be as remembered", () => readingAsRemembered(page, "zoom", el("map viewer")));
      await session.step(35, "And no rows should be selected", () => noneSelected(page));
      await session.step(36, "And the \"markers\" reading of map viewer should be 2426", () => readingIs(page, "markers", el("map viewer"), 2426));
      await session.step(37, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
