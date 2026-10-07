/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/trellis-plot/trellis-plot-scales-and-paging.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.trellis-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/biostructure.js';
import '../../../bindings/connections.js';
import '../../../bindings/flow.js';
import '../../../bindings/grid.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {shouldBe} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewerWith, clickArea, dragAreaBy, hasArea, hasNoArea, hoverArea, noErrors, pickFromAreaContextMenu, propertyShouldBe, readingAsRemembered, readingIs, readingNotAsRemembered, rememberReading, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {cellsWideTall, innerPropertyShouldBe} from '@datagrok-libraries/bdd/bindings/tiers/viewers/widgets';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Trellis plot axes without global scale, the Y slider, paging ends and packing", () => {
  const session = feature(test, "features/viewers/trellis-plot/trellis-plot-scales-and-paging.feature", import.meta.url);
  test("Trellis plot axes without global scale, the Y slider, paging ends and packing", {tag: ["@journey", "@viewers", "@realizes:viewers.trellis-plot"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 5, page);
    await session.step(17, "Given user is logged in", () => loggedIn(page));
    await session.step(18, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(19, "And user adds a trellis plot viewer with:", () => addViewerWith(page, "trellis plot", [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Viewer Type","Scatter plot"]]);
    await session.step(23, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
    await session.step(24, "And \"Global Scale\" property of trellis plot viewer should be \"false\"", () => propertyShouldBe(page, "Global Scale", el("trellis plot viewer"), "false"));
    await run.scenario("Without Global Scale the axis settings draw nothing, and Allow Zoom starts off", async () => {
      await session.step(27, "Then \"allowZoom\" inner property of trellis plot viewer should be \"false\"", () => innerPropertyShouldBe(page, "allowZoom", el("trellis plot viewer"), "false"));
      await session.step(28, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Show X Axes","Always"],["Show Y Axes","Always"],["Show Range Sliders","true"]]), [["Show X Axes","Always"],["Show Y Axes","Always"],["Show Range Sliders","true"]]);
      await session.step(32, "Then trellis plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(33, "And trellis plot viewer should not have a \"y axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "y axis"));
      await session.step(34, "And the \"x axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "x axis sliders", el("trellis plot viewer"), 0));
      await session.step(35, "And the \"y axis sliders\" reading of trellis plot viewer should be 0", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 0));
      await session.step(36, "When user sets \"Global Scale\" property of trellis plot viewer to \"true\"", () => setProperty(page, "Global Scale", el("trellis plot viewer"), "true"));
      await session.step(37, "Then trellis plot viewer should have an \"x axis\" area", () => hasArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(38, "And trellis plot viewer should have a \"y axis\" area", () => hasArea(page, el("trellis plot viewer"), "y axis"));
      await session.step(39, "And the \"y axis sliders\" reading of trellis plot viewer should be 5", () => readingIs(page, "y axis sliders", el("trellis plot viewer"), 5));
      await session.step(40, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The shared Y slider re-bounds every cell and Reset puts them back exactly", async () => {
      await session.step(43, "When user hovers over the \"cell body F | Caucasian\" area of trellis plot viewer", () => hoverArea(page, "cell body F | Caucasian", el("trellis plot viewer")));
      await session.step(44, "Then trellis plot viewer should have a \"y range slider\" area", () => hasArea(page, el("trellis plot viewer"), "y range slider"));
      await session.step(45, "And trellis plot viewer should have a \"y range slider max handle\" area", () => hasArea(page, el("trellis plot viewer"), "y range slider max handle"));
      await session.step(46, "When user remembers the \"cell signature F | Caucasian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(47, "And user remembers the \"cell signature M | Asian\" reading of trellis plot viewer", () => rememberReading(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(48, "And user drags the \"y range slider max handle\" area of trellis plot viewer by 40 pixels to the down", () => dragAreaBy(page, "y range slider max handle", el("trellis plot viewer"), 40, "down"));
      await session.step(49, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(50, "And the \"cell signature M | Asian\" reading of trellis plot viewer should not be as remembered", () => readingNotAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(51, "When user picks \"Reset Inner Range Sliders\" from the context menu of the \"view\" area of trellis plot viewer", () => pickFromAreaContextMenu(page, "Reset Inner Range Sliders", "view", el("trellis plot viewer")));
      await session.step(52, "Then the \"cell signature F | Caucasian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature F | Caucasian", el("trellis plot viewer")));
      await session.step(53, "And the \"cell signature M | Asian\" reading of trellis plot viewer should be as remembered", () => readingAsRemembered(page, "cell signature M | Asian", el("trellis plot viewer")));
      await session.step(54, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["Global Scale","false"],["Show X Axes","Auto"],["Show Y Axes","Auto"],["Show Range Sliders","false"]]), [["Global Scale","false"],["Show X Axes","Auto"],["Show Y Axes","Auto"],["Show Range Sliders","false"]]);
      await session.step(59, "Then trellis plot viewer should not have an \"x axis\" area", () => hasNoArea(page, el("trellis plot viewer"), "x axis"));
      await session.step(60, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At entry every category fits, so (+) has nothing to add and (-) pages out", async () => {
      await session.step(63, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","DIS_POP"],["Pack Categories","false"]]), [["X Column Names","DIS_POP"],["Pack Categories","false"]]);
      await session.step(66, "Then the \"x categories\" reading of trellis plot viewer should be 6", () => readingIs(page, "x categories", el("trellis plot viewer"), 6));
      await session.step(67, "And the cells of trellis plot viewer should be 6 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 6, 4));
      await session.step(68, "And x plus icon should be disabled", () => shouldBe(page, el("x plus icon"), "disabled"));
      await session.step(69, "And x minus icon should be enabled", () => shouldBe(page, el("x minus icon"), "enabled"));
      await session.step(70, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(71, "Then the cells of trellis plot viewer should be 6 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 6, 4));
      await session.step(72, "When user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(73, "Then the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(74, "And x plus icon should be enabled", () => shouldBe(page, el("x plus icon"), "enabled"));
      await session.step(75, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(76, "Then the cells of trellis plot viewer should be 6 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 6, 4));
      await session.step(77, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("At the far end (+) stops adding and (-) still pages", async () => {
      await session.step(80, "When user sets \"X Column Names\" property of trellis plot viewer to \"SEX, DIS_POP\"", () => setProperty(page, "X Column Names", el("trellis plot viewer"), "SEX, DIS_POP"));
      await session.step(81, "Then the \"x categories\" reading of trellis plot viewer should be 12", () => readingIs(page, "x categories", el("trellis plot viewer"), 12));
      await session.step(82, "And the cells of trellis plot viewer should be 5 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 5, 4));
      await session.step(83, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(84, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(85, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(86, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(87, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(88, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(89, "And user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(90, "Then the cells of trellis plot viewer should be 12 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 12, 4));
      await session.step(91, "And x plus icon should be disabled", () => shouldBe(page, el("x plus icon"), "disabled"));
      await session.step(92, "When user clicks on the \"x plus\" area of trellis plot viewer", () => clickArea(page, "x plus", el("trellis plot viewer")));
      await session.step(93, "Then the cells of trellis plot viewer should be 12 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 12, 4));
      await session.step(94, "When user clicks on the \"x minus\" area of trellis plot viewer", () => clickArea(page, "x minus", el("trellis plot viewer")));
      await session.step(95, "Then the cells of trellis plot viewer should be 11 wide and 4 tall", () => cellsWideTall(page, el("trellis plot viewer"), 11, 4));
      await session.step(96, "And x plus icon should be enabled", () => shouldBe(page, el("x plus icon"), "enabled"));
      await session.step(97, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Packing a two-column X axis drops the combinations no row holds", async () => {
      await session.step(100, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","RACE, SEVERITY"],["Y Column Names","SEX"],["Pack Categories","true"]]), [["X Column Names","RACE, SEVERITY"],["Y Column Names","SEX"],["Pack Categories","true"]]);
      await session.step(104, "Then the \"x categories\" reading of trellis plot viewer should be 20", () => readingIs(page, "x categories", el("trellis plot viewer"), 20));
      await session.step(105, "And the \"x categories packed\" reading of trellis plot viewer should be 16", () => readingIs(page, "x categories packed", el("trellis plot viewer"), 16));
      await session.step(106, "And the \"x scroll handle share\" reading of trellis plot viewer should be 0.3125", () => readingIs(page, "x scroll handle share", el("trellis plot viewer"), 0.3125));
      await session.step(107, "When user sets \"Pack Categories\" property of trellis plot viewer to \"false\"", () => setProperty(page, "Pack Categories", el("trellis plot viewer"), "false"));
      await session.step(108, "Then the \"x categories packed\" reading of trellis plot viewer should be 20", () => readingIs(page, "x categories packed", el("trellis plot viewer"), 20));
      await session.step(109, "And the \"x scroll handle share\" reading of trellis plot viewer should be 0.25", () => readingIs(page, "x scroll handle share", el("trellis plot viewer"), 0.25));
      await session.step(110, "When user sets properties of trellis plot viewer:", () => setProperties(page, el("trellis plot viewer"), [["X Column Names","SEX"],["Y Column Names","RACE"],["Pack Categories","true"]]), [["X Column Names","SEX"],["Y Column Names","RACE"],["Pack Categories","true"]]);
      await session.step(114, "Then the \"cells\" reading of trellis plot viewer should be 8", () => readingIs(page, "cells", el("trellis plot viewer"), 8));
      await session.step(115, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
