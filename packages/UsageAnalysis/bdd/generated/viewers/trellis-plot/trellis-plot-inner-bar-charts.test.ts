/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-inner-bar-charts.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, public#3617, public#1547]
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
import {clickOn} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingAsRemembered, readingIs, readingNotAsRemembered, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {pickInnerViewer, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot — inner bar charts on spgi-100", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-inner-bar-charts.feature", import.meta.url);
  test("Trellis plot — inner bar charts on spgi-100", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot", "@realizes:public#3617", "@realizes:public#1547"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(15, "Given user is logged in", () => loggedIn(page));
    await session.step(16, "And user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await run.scenario("Include Nulls off redraws the inner bar charts, and on brings their picture back (public#3617)", async () => {
      await session.step(20, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["Y Column Names","Stereo Category"]]), [["Y Column Names","Stereo Category"]]);
      await session.step(22, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(23, "And user sets \"splitColumnName\" inner property of trellis plot viewer to \"Series\"", () => setInnerProperty(page, "splitColumnName", el("trellis plot viewer"), "Series"));
      await session.step(24, "Then the \"cells drawn\" reading of trellis plot viewer should be 15", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 15));
      await session.step(25, "When user remembers the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer", () => rememberReading(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
      await session.step(26, "And user sets \"includeNulls\" inner property of trellis plot viewer to \"false\"", () => setInnerProperty(page, "includeNulls", el("trellis plot viewer"), "false"));
      await session.step(27, "Then the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
      await session.step(28, "And the \"cells drawn\" reading of trellis plot viewer should be 15", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 15));
      await session.step(29, "When user sets \"includeNulls\" inner property of trellis plot viewer to \"true\"", () => setInnerProperty(page, "includeNulls", el("trellis plot viewer"), "true"));
      await session.step(30, "Then the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
      await session.step(31, "When user clicks on close icon of trellis plot viewer", () => clickOn(page, el("close icon of trellis plot viewer")));
      await session.step(32, "Then no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Inner bar charts draw complete on the first pass (public#1547)", async () => {
      await session.step(36, "Given user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","Scaffold Names"],["Y Column Names","Series"]]), [["X Column Names","Scaffold Names"],["Y Column Names","Series"]]);
      await session.step(39, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
      await session.step(40, "And user sets \"valueColumnName\" inner property of trellis plot viewer to \"CAST Idea ID\"", () => setInnerProperty(page, "valueColumnName", el("trellis plot viewer"), "CAST Idea ID"));
      await session.step(41, "And user sets \"splitColumnName\" inner property of trellis plot viewer to \"Stereo Category\"", () => setInnerProperty(page, "splitColumnName", el("trellis plot viewer"), "Stereo Category"));
      await session.step(42, "Then the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
      await session.step(43, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(44, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 6", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 6));
      await session.step(45, "When user remembers the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer", () => rememberReading(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
      await session.step(46, "And user sets \"barSortOrder\" inner property of trellis plot viewer to \"desc\"", () => setInnerProperty(page, "barSortOrder", el("trellis plot viewer"), "desc"));
      await session.step(47, "Then the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
      await session.step(48, "When user sets \"barSortOrder\" inner property of trellis plot viewer to \"asc\"", () => setInnerProperty(page, "barSortOrder", el("trellis plot viewer"), "asc"));
      await session.step(49, "Then the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
      await session.step(50, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
      await session.step(51, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
