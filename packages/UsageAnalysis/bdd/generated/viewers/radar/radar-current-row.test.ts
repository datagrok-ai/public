/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/radar/radar-current-row.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [charts.viewer.radar]
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
import {check, clickOn, isExpanded, shouldBe, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {currentRowIs} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {autostartsCompleted, openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {areaColor, areaNotColor, clickArea, noErrors, propertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Radar highlights the current row's line with Show Current Row", () => {
  const session = feature(test, "features/viewers/radar/radar-current-row.feature", import.meta.url);
  test("Show Current Row draws the current row's line in green, and off it is gone", {tag: ["@viewers", "@realizes:charts.viewer.radar"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(11, "Given user is logged in", () => loggedIn(page));
    await session.step(12, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(13, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(14, "When user clicks on \"Add viewer\" icon", () => clickOn(page, el("\"Add viewer\" icon")));
    await session.step(15, "Then \"Add Viewer\" dialog should be visible", () => shouldBe(page, el("\"Add Viewer\" dialog"), "visible"));
    await session.step(16, "When user clicks on first \"Radar\" card in \"Add Viewer\" dialog", () => clickOn(page, el("first \"Radar\" card in \"Add Viewer\" dialog")));
    await session.step(17, "Then radar viewer should be visible", () => shouldBe(page, el("radar viewer"), "visible"));
    await session.step(18, "And \"Values\" property of radar viewer should be \"AGE, HEIGHT, WEIGHT\"", () => propertyShouldBe(page, "Values", el("radar viewer"), "AGE, HEIGHT, WEIGHT"));
    await session.step(21, "When user clicks on the \"cell 5 of AGE\" area of grid", () => clickArea(page, "cell 5 of AGE", el("grid")));
    await session.step(22, "Then row 5 should be current", () => currentRowIs(page, 5));
    await session.step(23, "And the \"view\" area of radar viewer should contain the color \"#00FF00\"", () => areaColor(page, "view", el("radar viewer"), "#00FF00"));
    await session.step(24, "When user clicks on settings icon of radar viewer", () => clickOn(page, el("settings icon of radar viewer")));
    await session.step(25, "Given \"Selection\" category in context panel is expanded", () => isExpanded(page, el("\"Selection\" category in context panel")));
    await session.step(26, "When user unchecks \"Show Current Row\" property in context panel", () => uncheck(page, el("\"Show Current Row\" property in context panel")));
    await session.step(27, "Then \"Show Current Row\" property of radar viewer should be \"false\"", () => propertyShouldBe(page, "Show Current Row", el("radar viewer"), "false"));
    await session.step(28, "And the \"view\" area of radar viewer should not contain the color \"#00FF00\"", () => areaNotColor(page, "view", el("radar viewer"), "#00FF00"));
    await session.step(29, "When user checks \"Show Current Row\" property in context panel", () => check(page, el("\"Show Current Row\" property in context panel")));
    await session.step(30, "Then the \"view\" area of radar viewer should contain the color \"#00FF00\"", () => areaColor(page, "view", el("radar viewer"), "#00FF00"));
    await session.step(31, "And no errors should have been logged", () => noErrors(page));
  });
});
