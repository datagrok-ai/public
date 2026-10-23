/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/box-plot/box-plot-hover-and-date-categories.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.box-plot]
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
import {clickOn, shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, areaColor, hasArea, hasNoArea, hoverArea, noErrors, notRepainted, pointerAway, repaintedBy, setProperties, setProperty, takeSnapshot} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Box plot hover highlights and date categories", () => {
  const session = feature(test, "features/viewers/box-plot/box-plot-hover-and-date-categories.feature", import.meta.url);
  test("Box plot hover highlights and date categories", {tag: ["@journey", "@viewers", "@realizes:viewers.box-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 3, page);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(13, "And user adds a box plot viewer with:", () => addViewerWith(page, "box plot", [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]), [["Value","AGE"],["Category 1","RACE"],["Marker Size","10"]]);
    await run.scenario("Hovering a category's plot column highlights its row group, and Show Mouse Over Row Group switches it off", async () => {
      await session.step(19, "When user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(20, "And user hovers over the \"category Caucasian\" area of box plot viewer", () => hoverArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(21, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(22, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(23, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"false\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "false"));
      await session.step(24, "And user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(25, "And user hovers over the \"category Caucasian\" area of box plot viewer", () => hoverArea(page, "category Caucasian", el("box plot viewer")));
      await session.step(26, "Then box plot viewer should not have repainted", () => notRepainted(page, el("box plot viewer")));
      await session.step(27, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(28, "And user sets \"Show Mouse Over Row Group\" property of box plot viewer to \"true\"", () => setProperty(page, "Show Mouse Over Row Group", el("box plot viewer"), "true"));
      await session.step(29, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A hover on another viewer's bar keeps the hovered group's own marker colours (github-3065)", async () => {
      await session.step(32, "When user adds a bar chart viewer with:", () => addViewerWith(page, "bar chart", [["Split","RACE"]]), [["Split","RACE"]]);
      await session.step(34, "And user sets \"Marker Color Column\" property of box plot viewer to \"RACE\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), "RACE"));
      await session.step(35, "Then the \"Caucasian values\" area of box plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "Caucasian values", el("box plot viewer"), "#2CA02C"));
      await session.step(36, "When user takes a snapshot of box plot viewer", () => takeSnapshot(page, el("box plot viewer")));
      await session.step(37, "And user hovers over the \"bar Caucasian\" area of bar chart viewer", () => hoverArea(page, "bar Caucasian", el("bar chart viewer")));
      await session.step(38, "Then box plot viewer should have repainted by at least 2000 pixels", () => repaintedBy(page, el("box plot viewer"), 2000));
      await session.step(39, "And the \"Caucasian values\" area of box plot viewer should contain the color \"#2CA02C\"", () => areaColor(page, "Caucasian values", el("box plot viewer"), "#2CA02C"));
      await session.step(40, "When user moves the pointer away from grid", () => pointerAway(page, el("grid")));
      await session.step(41, "And user clicks on close icon of bar chart viewer", () => clickOn(page, el("close icon of bar chart viewer")));
      await session.step(42, "Then bar chart viewer should be absent", () => shouldBe(page, el("bar chart viewer"), "absent"));
      await session.step(43, "When user sets \"Marker Color Column\" property of box plot viewer to \"\"", () => setProperty(page, "Marker Color Column", el("box plot viewer"), ""));
      await session.step(44, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A date category mapped to month and to quarter splits into months and quarters", async () => {
      await session.step(47, "When user sets \"Category 1\" property of box plot viewer to \"STARTED\"", () => setProperty(page, "Category 1", el("box plot viewer"), "STARTED"));
      await session.step(48, "And user sets \"Category 1 Map\" property of box plot viewer to \"month\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "month"));
      await session.step(49, "Then box plot viewer should have a \"category January\" area", () => hasArea(page, el("box plot viewer"), "category January"));
      await session.step(50, "And box plot viewer should have a \"category December\" area", () => hasArea(page, el("box plot viewer"), "category December"));
      await session.step(51, "And box plot viewer should have a \"June values\" area", () => hasArea(page, el("box plot viewer"), "June values"));
      await session.step(52, "And box plot viewer should not have a \"category Q1\" area", () => hasNoArea(page, el("box plot viewer"), "category Q1"));
      await session.step(53, "When user sets \"Category 1 Map\" property of box plot viewer to \"quarter\"", () => setProperty(page, "Category 1 Map", el("box plot viewer"), "quarter"));
      await session.step(54, "Then box plot viewer should have a \"category Q1\" area", () => hasArea(page, el("box plot viewer"), "category Q1"));
      await session.step(55, "And box plot viewer should have a \"category Q4\" area", () => hasArea(page, el("box plot viewer"), "category Q4"));
      await session.step(56, "And box plot viewer should have a \"Q3 values\" area", () => hasArea(page, el("box plot viewer"), "Q3 values"));
      await session.step(57, "And box plot viewer should not have a \"category January\" area", () => hasNoArea(page, el("box plot viewer"), "category January"));
      await session.step(58, "When user sets properties of box plot viewer:", () => setProperties(page, el("box plot viewer"), [["Category 1 Map",""],["Category 1","RACE"]]), [["Category 1 Map",""],["Category 1","RACE"]]);
      await session.step(61, "Then box plot viewer should have a \"category Caucasian\" area", () => hasArea(page, el("box plot viewer"), "category Caucasian"));
      await session.step(62, "And box plot viewer should not have a \"category Q1\" area", () => hasNoArea(page, el("box plot viewer"), "category Q1"));
      await session.step(63, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
