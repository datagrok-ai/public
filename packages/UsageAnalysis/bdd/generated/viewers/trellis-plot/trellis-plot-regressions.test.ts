/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-regressions.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot, public#2144, public#3617, public#1547]
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
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, noErrors, readingAsRemembered, readingIs, readingNotAsRemembered, readingReads, rememberReading} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {innerPropertyShouldBe, pickInnerViewer, setInnerProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot regressions of closed issues", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-regressions.feature", import.meta.url);
  test("An inner bar chart keeps its own setting through a pie chart and back (public#2144)", {tag: ["@viewers", "@realizes:viewers.trellis-plot", "@realizes:public#2144"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(20, "Given user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(21, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"]]), [["X Column Names","SEX"],["Y Column Names","RACE"]]);
    await session.step(24, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
    await session.step(25, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
    await session.step(26, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(27, "And user sets \"valueAggrType\" inner property of trellis plot viewer to \"avg\"", () => setInnerProperty(page, "valueAggrType", el("trellis plot viewer"), "avg"));
    await session.step(28, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(29, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(30, "And user picks \"Pie chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Pie chart", el("trellis plot viewer")));
    await session.step(31, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Pie chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Pie chart"));
    await session.step(32, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(33, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
    await session.step(34, "Then the \"inner viewer type\" reading of trellis plot viewer should be \"Bar chart\"", () => readingReads(page, "inner viewer type", el("trellis plot viewer"), "Bar chart"));
    await session.step(35, "And \"valueAggrType\" inner property of trellis plot viewer should be \"avg\"", () => innerPropertyShouldBe(page, "valueAggrType", el("trellis plot viewer"), "avg"));
    await session.step(36, "And the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
    await session.step(37, "And no errors should have been logged", () => noErrors(page));
  });
  test("Include Nulls off redraws the inner bar charts, and on brings their picture back (public#3617)", {tag: ["@viewers", "@realizes:viewers.trellis-plot", "@realizes:public#3617"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(41, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(42, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["Y Column Names","Stereo Category"]]), [["Y Column Names","Stereo Category"]]);
    await session.step(44, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
    await session.step(45, "And user sets \"splitColumnName\" inner property of trellis plot viewer to \"Series\"", () => setInnerProperty(page, "splitColumnName", el("trellis plot viewer"), "Series"));
    await session.step(46, "Then the \"cells drawn\" reading of trellis plot viewer should be 15", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 15));
    await session.step(47, "When user remembers the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer", () => rememberReading(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
    await session.step(48, "And user sets \"includeNulls\" inner property of trellis plot viewer to \"false\"", () => setInnerProperty(page, "includeNulls", el("trellis plot viewer"), "false"));
    await session.step(49, "Then the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
    await session.step(50, "And the \"cells drawn\" reading of trellis plot viewer should be 15", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 15));
    await session.step(51, "When user sets \"includeNulls\" inner property of trellis plot viewer to \"true\"", () => setInnerProperty(page, "includeNulls", el("trellis plot viewer"), "true"));
    await session.step(52, "Then the \"cell signature Chemist 13 | R_ONE\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature Chemist 13 | R_ONE", el("trellis plot viewer")));
    await session.step(53, "And no errors should have been logged", () => noErrors(page));
  });
  test("Inner bar charts draw complete on the first pass (public#1547)", {tag: ["@viewers", "@realizes:viewers.trellis-plot", "@realizes:public#1547"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(16, "Given user is logged in", () => loggedIn(page));
    await session.step(57, "Given user opens spgi dataset", () => openDataset(page, ds("spgi")));
    await session.step(58, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","Scaffold Names"],["Y Column Names","Series"]]), [["X Column Names","Scaffold Names"],["Y Column Names","Series"]]);
    await session.step(61, "When user picks \"Bar chart\" in the viewer selector of trellis plot viewer", () => pickInnerViewer(page, "Bar chart", el("trellis plot viewer")));
    await session.step(62, "And user sets \"valueColumnName\" inner property of trellis plot viewer to \"CAST Idea ID\"", () => setInnerProperty(page, "valueColumnName", el("trellis plot viewer"), "CAST Idea ID"));
    await session.step(63, "And user sets \"splitColumnName\" inner property of trellis plot viewer to \"Stereo Category\"", () => setInnerProperty(page, "splitColumnName", el("trellis plot viewer"), "Stereo Category"));
    await session.step(64, "Then the \"cells drawn\" reading of trellis plot viewer should be 6", () => readingIs(page, "cells drawn", el("trellis plot viewer"), 6));
    await session.step(65, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
    await session.step(66, "And the \"distinct cell signatures\" reading of trellis plot viewer should be 6", () => readingIs(page, "distinct cell signatures", el("trellis plot viewer"), 6));
    await session.step(67, "When user remembers the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer", () => rememberReading(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
    await session.step(68, "And user sets \"barSortOrder\" inner property of trellis plot viewer to \"desc\"", () => setInnerProperty(page, "barSortOrder", el("trellis plot viewer"), "desc"));
    await session.step(69, "Then the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
    await session.step(70, "When user sets \"barSortOrder\" inner property of trellis plot viewer to \"asc\"", () => setInnerProperty(page, "barSortOrder", el("trellis plot viewer"), "asc"));
    await session.step(71, "Then the \"cell signature TRISUBSTITUTED | Triazoles\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature TRISUBSTITUTED | Triazoles", el("trellis plot viewer")));
    await session.step(72, "And the \"blank cells\" reading of trellis plot viewer should be 0", () => readingIs(page, "blank cells", el("trellis plot viewer"), 0));
    await session.step(73, "And no errors should have been logged", () => noErrors(page));
  });
});
