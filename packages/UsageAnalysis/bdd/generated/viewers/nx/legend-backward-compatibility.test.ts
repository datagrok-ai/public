/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/nx/legend-backward-compatibility.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [GROK-17281]
--- */
import {test} from '@playwright/test';
import '../../../bindings/connections.js';
import '../../../bindings/grid.js';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {noViewerError} from '../../../bindings/nx.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {uploadThrough} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, browsePanelOpen, openDataset, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {legendSide, noBalloons, noErrors, propertyShouldBe, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("A layout saved before the legend position existed puts every legend on the right", () => {
  const session = feature(test, "features/viewers/nx/legend-backward-compatibility.feature", import.meta.url);
  test("Every viewer of the old layout has its legend on the right", {tag: ["@viewers", "@realizes:GROK-17281"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(18, "Given user is logged in", () => loggedIn(page));
    await session.step(19, "And simple mode is off", () => simpleModeOff(page));
    await session.step(20, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(21, "And user opens spgi-3624 dataset", () => openDataset(page, ds("spgi-3624")));
    await session.step(22, "And the browse panel is open", () => browsePanelOpen(page));
    await session.step(23, "When user uploads \"fixtures/nx/spgi-legend-position-old.layout\" through \"Open local file\" icon inside browse toolbar", () => uploadThrough(page, "fixtures/nx/spgi-legend-position-old.layout", el("\"Open local file\" icon inside browse toolbar")));
    await session.step(24, "Then the open tableview should have 1 pie chart viewer", () => viewerCount(page, 1, "pie chart"));
    await session.step(25, "And the open tableview should have 1 scatter plot viewer", () => viewerCount(page, 1, "scatter plot"));
    await session.step(26, "And the open tableview should have 1 bar chart viewer", () => viewerCount(page, 1, "bar chart"));
    await session.step(27, "And the open tableview should have 1 line chart viewer", () => viewerCount(page, 1, "line chart"));
    await session.step(28, "And the open tableview should have 1 box plot viewer", () => viewerCount(page, 1, "box plot"));
    await session.step(29, "And \"legendPosition\" property of pie chart viewer should be \"Right\"", () => propertyShouldBe(page, "legendPosition", el("pie chart viewer"), "Right"));
    await session.step(30, "And \"legendPosition\" property of scatter plot viewer should be \"Right\"", () => propertyShouldBe(page, "legendPosition", el("scatter plot viewer"), "Right"));
    await session.step(31, "And \"legendPosition\" property of bar chart viewer should be \"Right\"", () => propertyShouldBe(page, "legendPosition", el("bar chart viewer"), "Right"));
    await session.step(32, "And \"legendPosition\" property of line chart viewer should be \"Right\"", () => propertyShouldBe(page, "legendPosition", el("line chart viewer"), "Right"));
    await session.step(33, "And \"legendPosition\" property of box plot viewer should be \"Right\"", () => propertyShouldBe(page, "legendPosition", el("box plot viewer"), "Right"));
    await session.step(34, "And the legend of scatter plot viewer should be on the right", () => legendSide(page, el("scatter plot viewer"), "right"));
    await session.step(35, "And the legend of bar chart viewer should be on the right", () => legendSide(page, el("bar chart viewer"), "right"));
    await session.step(36, "And the legend of line chart viewer should be on the right", () => legendSide(page, el("line chart viewer"), "right"));
    await session.step(37, "And the legend of box plot viewer should be on the right", () => legendSide(page, el("box plot viewer"), "right"));
    await session.step(38, "When user sets \"legendVisibility\" property of pie chart viewer to \"Always\"", () => setProperty(page, "legendVisibility", el("pie chart viewer"), "Always"));
    await session.step(39, "Then the legend of pie chart viewer should be on the right", () => legendSide(page, el("pie chart viewer"), "right"));
    await session.step(40, "And no viewer of the current view should report an error", () => noViewerError(page));
    await session.step(41, "And no errors should have been logged", () => noErrors(page));
    await session.step(42, "And no error or warning balloon should have been shown", () => noBalloons(page));
  });
});
