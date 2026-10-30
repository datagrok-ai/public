/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/matrix-plot/matrix-plot-filter-panel-and-axes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.matrix-plot]
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
import {check, clickOn, hoverOver, isExpanded, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {addCategoricalFilter, filterPasses, openEmptyFilterPanel} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, propertyShouldBe, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Matrix plot — a filter-panel filter in every cell, and the Show X Axes checkbox a person clicks", () => {
  const session = feature(test, "features/viewers/matrix-plot/matrix-plot-filter-panel-and-axes.feature", import.meta.url);
  test("Matrix plot — a filter-panel filter in every cell, and the Show X Axes checkbox a person clicks", {tag: ["@journey", "@viewers", "@realizes:viewers.matrix-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(14, "Given user is logged in", () => loggedIn(page));
    await session.step(15, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(16, "And user adds a matrix plot viewer", () => addViewer(page, "matrix plot"));
    await session.step(17, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(18, "And the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
    await session.step(19, "And the \"cells\" reading of matrix plot viewer should be 16", () => readingIs(page, "cells", el("matrix plot viewer"), 16));
    await run.scenario("A SEX = M card narrows every cell to its own rows, and closing it restores them", async () => {
      await session.step(24, "When user opens an empty filter panel", () => openEmptyFilterPanel(page));
      await session.step(25, "Then the \"distinct cell signatures\" reading of matrix plot viewer should be 16", () => readingIs(page, "distinct cell signatures", el("matrix plot viewer"), 16));
      await session.step(26, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(27, "When user remembers the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer", () => rememberReading(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(28, "And user adds a categorical filter on \"SEX\" keeping \"M\"", () => addCategoricalFilter(page, "SEX", "M"));
      await session.step(29, "Then 447 rows should pass the filter", () => filterPasses(page, 447));
      await session.step(30, "And the \"rows shown\" reading of matrix plot viewer should be 447", () => readingIs(page, "rows shown", el("matrix plot viewer"), 447));
      await session.step(31, "And the \"cell rows shown of AGE x AGE\" reading of matrix plot viewer should be 447", () => readingIs(page, "cell rows shown of AGE x AGE", el("matrix plot viewer"), 447));
      await session.step(32, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 373", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 373));
      await session.step(33, "And the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(34, "When user hovers over \"SEX\" filter card", () => hoverOver(page, el("\"SEX\" filter card")));
      await session.step(35, "And user clicks on close of \"SEX\" filter card", () => clickOn(page, el("close of \"SEX\" filter card")));
      await session.step(36, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
      await session.step(37, "And the \"rows shown\" reading of matrix plot viewer should be 1000", () => readingIs(page, "rows shown", el("matrix plot viewer"), 1000));
      await session.step(38, "And the \"cell rows shown of HEIGHT x AGE\" reading of matrix plot viewer should be 872", () => readingIs(page, "cell rows shown of HEIGHT x AGE", el("matrix plot viewer"), 872));
      await session.step(39, "And the \"cell signature HEIGHT x AGE\" reading of matrix plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature HEIGHT x AGE", el("matrix plot viewer")));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The Show X Axes checkbox in the context panel takes the axis strip away and brings it back (GROK-19106)", async () => {
      await session.step(43, "Then the \"auto layout\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "auto layout", el("matrix plot viewer"), "true"));
      await session.step(44, "And the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(45, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(46, "When user clicks on settings icon of matrix plot viewer", () => clickOn(page, el("settings icon of matrix plot viewer")));
      await session.step(47, "Given \"Axes\" category in context panel is expanded", () => isExpanded(page, el("\"Axes\" category in context panel")));
      await session.step(48, "When user unchecks \"Show X Axes\" property in context panel", () => uncheck(page, el("\"Show X Axes\" property in context panel")));
      await session.step(49, "Then \"showXAxes\" property of matrix plot viewer should be \"false\"", () => propertyShouldBe(page, "showXAxes", el("matrix plot viewer"), "false"));
      await session.step(50, "And the \"x axis shown\" reading of matrix plot viewer should be \"false\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "false"));
      await session.step(51, "And matrix plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(52, "And the \"y axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "y axis shown", el("matrix plot viewer"), "true"));
      await session.step(53, "When user checks \"Show X Axes\" property in context panel", () => check(page, el("\"Show X Axes\" property in context panel")));
      await session.step(54, "Then the \"x axis shown\" reading of matrix plot viewer should be \"true\"", () => readingReads(page, "x axis shown", el("matrix plot viewer"), "true"));
      await session.step(55, "And matrix plot viewer should have an \"x axis\" area", () => hasArea(page, el("matrix plot viewer"), "x axis"));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
