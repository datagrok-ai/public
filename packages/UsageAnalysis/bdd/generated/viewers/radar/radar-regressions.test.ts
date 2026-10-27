/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/radar/radar-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.radar, public#3415, public#1472]
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
import {clickOn, shouldHaveValue} from '@datagrok-libraries/bdd/bindings/common/steps';
import {autostartsCompleted, contextPanelOpen, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, boundTable, noErrors, propertyShouldBe, setProperty, viewerCount} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickFromViewerMenu} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Radar regressions — the Table property and Pick Up / Apply on a JS viewer", () => {
  const session = feature(test, "features/viewers/radar/radar-regressions.feature", import.meta.url);
  test("The Table property of a JS viewer names its table (public#3415)", {tag: ["@viewers", "@realizes:charts.viewer.radar", "@realizes:public#3415"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(18, "Given user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(19, "Then radar viewer should be bound to table \"demog-1000\"", () => boundTable(page, el("radar viewer"), "demog-1000"));
    await session.step(20, "When user clicks on settings icon of radar viewer", () => clickOn(page, el("settings icon of radar viewer")));
    await session.step(21, "Then \"Table\" property in context panel should have value \"demog-1000\"", () => shouldHaveValue(page, el("\"Table\" property in context panel"), "demog-1000"));
    await session.step(22, "And no errors should have been logged", () => noErrors(page));
  });
  test("Pick Up / Apply copies a setting between two Radar viewers, and the two stay apart (public#1472)", {tag: ["@viewers", "@realizes:charts.viewer.radar", "@realizes:public#1472"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(10, "Given user is logged in", () => loggedIn(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And the context panel is open", () => contextPanelOpen(page));
    await session.step(27, "Given user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(28, "And user adds a radar viewer", () => addViewer(page, "radar"));
    await session.step(29, "Then the open tableview should have 2 radar viewers", () => viewerCount(page, 2, "radar"));
    await session.step(30, "When user sets \"normalization\" property of first radar viewer to \"Global\"", () => setProperty(page, "normalization", el("first radar viewer"), "Global"));
    await session.step(31, "Then \"normalization\" property of second radar viewer should be \"Column\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Column"));
    await session.step(32, "When user picks \"Pick Up / Apply > Pick Up\" from the viewer menu of first radar viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Pick Up", el("first radar viewer")));
    await session.step(33, "And user picks \"Pick Up / Apply > Apply\" from the viewer menu of second radar viewer", () => pickFromViewerMenu(page, "Pick Up / Apply > Apply", el("second radar viewer")));
    await session.step(34, "Then \"normalization\" property of second radar viewer should be \"Global\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Global"));
    await session.step(35, "When user sets \"normalization\" property of first radar viewer to \"Column\"", () => setProperty(page, "normalization", el("first radar viewer"), "Column"));
    await session.step(36, "Then \"normalization\" property of second radar viewer should be \"Global\"", () => propertyShouldBe(page, "normalization", el("second radar viewer"), "Global"));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
});
