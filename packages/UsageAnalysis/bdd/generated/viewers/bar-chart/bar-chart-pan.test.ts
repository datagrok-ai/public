/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/bar-chart/bar-chart-pan.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.bar-chart]
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
import {zoomCategories} from '../../../bindings/bar-chart.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {noneSelected} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, dragAreaBy, hasArea, hasNoArea, noErrors, readingIs} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {doubleClickEmptySpace} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Bar chart panning through zoomed categories", () => {
  const session = feature(test, "features/viewers/bar-chart/bar-chart-pan.feature", import.meta.url);
  test("Bar chart panning through zoomed categories", {tag: ["@journey", "@viewers", "@realizes:viewers.bar-chart"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 1, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(16, "And user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"]]), [["Split","Primary Series Name"],["Value","CAST Idea ID"],["Value Aggr Type","count"]]);
    await session.step(20, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
    await run.scenario("A plain drag pans the zoomed categories, and a double click resets the view", async () => {
      await session.step(23, "When user zooms into the categories from the \"bar Pyrrolidines\" area to the \"bar Triazoles\" area of bar chart viewer", () => zoomCategories(page, "bar Pyrrolidines", "bar Triazoles", el("bar chart viewer")));
      await session.step(24, "Then bar chart viewer should have a \"bar Triazoles\" area", () => hasArea(page, el("bar chart viewer"), "bar Triazoles"));
      await session.step(25, "And the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(26, "And bar chart viewer should not have a \"bar Aminopiperidines\" area", () => hasNoArea(page, el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(27, "When user drags the \"view\" area of bar chart viewer by 200 pixels to the up", () => dragAreaBy(page, "view", el("bar chart viewer"), 200, "up"));
      await session.step(28, "Then bar chart viewer should have a \"bar Aminopiperidines\" area", () => hasArea(page, el("bar chart viewer"), "bar Aminopiperidines"));
      await session.step(29, "And bar chart viewer should not have a \"bar Triazoles\" area", () => hasNoArea(page, el("bar chart viewer"), "bar Triazoles"));
      await session.step(30, "And no rows should be selected", () => noneSelected(page));
      await session.step(31, "When user drags the \"view\" area of bar chart viewer by 200 pixels to the up", () => dragAreaBy(page, "view", el("bar chart viewer"), 200, "up"));
      await session.step(32, "Then bar chart viewer should have a \"bar Diazabicyclooctane\" area", () => hasArea(page, el("bar chart viewer"), "bar Diazabicyclooctane"));
      await session.step(33, "And the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(34, "And bar chart viewer should not have a \"bar Pyrrolidines\" area", () => hasNoArea(page, el("bar chart viewer"), "bar Pyrrolidines"));
      await session.step(35, "When user drags the \"view\" area of bar chart viewer by 400 pixels to the down", () => dragAreaBy(page, "view", el("bar chart viewer"), 400, "down"));
      await session.step(36, "Then bar chart viewer should have a \"bar Triazoles\" area", () => hasArea(page, el("bar chart viewer"), "bar Triazoles"));
      await session.step(37, "And the \"bars\" reading of bar chart viewer should be 2", () => readingIs(page, "bars", el("bar chart viewer"), 2));
      await session.step(38, "And bar chart viewer should not have a \"bar Diazabicyclooctane\" area", () => hasNoArea(page, el("bar chart viewer"), "bar Diazabicyclooctane"));
      await session.step(39, "When user double-clicks on empty plot space of bar chart viewer", () => doubleClickEmptySpace(page, el("bar chart viewer")));
      await session.step(40, "Then the \"bars\" reading of bar chart viewer should be 5", () => readingIs(page, "bars", el("bar chart viewer"), 5));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
