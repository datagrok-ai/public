/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/viewers/correlation-plot/correlation-plot.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [viewers.correlation-plot]
--- */
import {test} from '@playwright/test';
import '../../../bindings/spaces.js';
import '../../../bindings/tile-viewer.js';
import '../../../bindings/trellis-plot.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {correlationMatches} from '../../../bindings/correlation-plot.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {addCalculated, filterPasses, removeColumn} from '@datagrok-libraries/bdd/bindings/platform/data';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {addViewer, hasArea, hasNoArea, noErrors, painted, readingBetween, readingDoesNotRead, readingIs, readingReads, readingsDiffer, readingsEqual, repainted, setProperties, setProperty} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Correlation plot — the matrix and the numbers in it", () => {
  const session = feature(test, "features/viewers/correlation-plot/correlation-plot.feature", import.meta.url);
  test("Correlation plot — the matrix and the numbers in it", {tag: ["@journey", "@viewers", "@realizes:viewers.correlation-plot", "@known-failure"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 9, page);
    await session.step(21, "Given user is logged in", () => loggedIn(page));
    await session.step(22, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(23, "And user adds a correlation plot viewer", () => addViewer(page, "correlation plot"));
    await session.step(24, "Then 1000 rows should pass the filter", () => filterPasses(page, 1000));
    await session.step(25, "And the \"rows shown\" reading of correlation plot viewer should be 1000", () => readingIs(page, "rows shown", el("correlation plot viewer"), 1000));
    await session.step(26, "And the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
    await session.step(27, "And the \"y columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
    await session.step(28, "And the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
    await session.step(29, "And the \"correlation type\" reading of correlation plot viewer should be \"Pearson\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Pearson"));
    await session.step(30, "And the \"show pearson r\" reading of correlation plot viewer should be \"true\"", () => readingReads(page, "show pearson r", el("correlation plot viewer"), "true"));
    await session.step(31, "And correlation plot viewer should be painted", () => painted(page, el("correlation plot viewer")));
    await run.scenario("Both axes are the table's numerical columns, and only those", async () => {
      await session.step(34, "Then the \"numerical columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "numerical columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(35, "And the \"x columns\" and \"numerical columns\" readings of correlation plot viewer should be the same", () => readingsEqual(page, "x columns", "numerical columns", el("correlation plot viewer")));
      await session.step(36, "And the \"y columns\" and \"numerical columns\" readings of correlation plot viewer should be the same", () => readingsEqual(page, "y columns", "numerical columns", el("correlation plot viewer")));
      await session.step(37, "And correlation plot viewer should have a \"cell HEIGHT x AGE\" area", () => hasArea(page, el("correlation plot viewer"), "cell HEIGHT x AGE"));
      await session.step(38, "And correlation plot viewer should not have a \"cell SEX x AGE\" area", () => hasNoArea(page, el("correlation plot viewer"), "cell SEX x AGE"));
      await session.step(39, "And correlation plot viewer should not have a \"cell RACE x AGE\" area", () => hasNoArea(page, el("correlation plot viewer"), "cell RACE x AGE"));
      await session.step(40, "And the \"columns shown\" reading of correlation plot viewer should be 6", () => readingIs(page, "columns shown", el("correlation plot viewer"), 6));
      await session.step(41, "And the \"column order\" reading of correlation plot viewer should be \"__t, __name, AGE, HEIGHT, WEIGHT, STARTED\"", () => readingReads(page, "column order", el("correlation plot viewer"), "__t, __name, AGE, HEIGHT, WEIGHT, STARTED"));
      await session.step(42, "And correlation plot viewer should have a \"row header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "row header HEIGHT"));
      await session.step(43, "And correlation plot viewer should have a \"type header HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "type header HEIGHT"));
      await session.step(44, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Every off-diagonal cell holds its pair's coefficient, and DG.Stats agrees", async () => {
      await session.step(47, "Then the correlation of \"HEIGHT\" and \"AGE\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "HEIGHT", "AGE", el("correlation plot viewer"), "Pearson"));
      await session.step(48, "And the correlation of \"WEIGHT\" and \"AGE\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "WEIGHT", "AGE", el("correlation plot viewer"), "Pearson"));
      await session.step(49, "And the correlation of \"WEIGHT\" and \"HEIGHT\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "WEIGHT", "HEIGHT", el("correlation plot viewer"), "Pearson"));
      await session.step(50, "And the correlation of \"STARTED\" and \"AGE\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "STARTED", "AGE", el("correlation plot viewer"), "Pearson"));
      await session.step(51, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(52, "And the \"correlation of WEIGHT and AGE\" reading of correlation plot viewer should be between 0.0647 and 0.0649", () => readingBetween(page, "correlation of WEIGHT and AGE", el("correlation plot viewer"), 0.0647, 0.0649));
      await session.step(53, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4124 and 0.4125", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4124, 0.4125));
      await session.step(54, "And the \"correlation of STARTED and AGE\" reading of correlation plot viewer should be between -0.0090 and -0.0089", () => readingBetween(page, "correlation of STARTED and AGE", el("correlation plot viewer"), -0.009, -0.0089));
      await session.step(55, "And the \"correlation of AGE and HEIGHT\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of AGE and HEIGHT", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(56, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The diagonal draws a histogram instead of a coefficient of one", async () => {
      await session.step(59, "Then the \"cell type of AGE x AGE\" reading of correlation plot viewer should be \"histogram\"", () => readingReads(page, "cell type of AGE x AGE", el("correlation plot viewer"), "histogram"));
      await session.step(60, "And the \"cell type of HEIGHT x HEIGHT\" reading of correlation plot viewer should be \"histogram\"", () => readingReads(page, "cell type of HEIGHT x HEIGHT", el("correlation plot viewer"), "histogram"));
      await session.step(61, "And the \"cell type of HEIGHT x AGE\" reading of correlation plot viewer should be \"correlation\"", () => readingReads(page, "cell type of HEIGHT x AGE", el("correlation plot viewer"), "correlation"));
      await session.step(62, "And correlation plot viewer should have a \"cell AGE x AGE\" area", () => hasArea(page, el("correlation plot viewer"), "cell AGE x AGE"));
      await session.step(63, "And the \"text of cell AGE x AGE\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell AGE x AGE", el("correlation plot viewer"), ""));
      await session.step(64, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(65, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("What the cell drew is the coefficient rounded to two decimals", async () => {
      await session.step(68, "Then the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(69, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(70, "And the \"text of cell WEIGHT x HEIGHT\" reading of correlation plot viewer should be \"0.41\"", () => readingReads(page, "text of cell WEIGHT x HEIGHT", el("correlation plot viewer"), "0.41"));
      await session.step(71, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4124 and 0.4125", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4124, 0.4125));
      await session.step(72, "And the \"text of cell STARTED x AGE\" reading of correlation plot viewer should be \"-0.01\"", () => readingReads(page, "text of cell STARTED x AGE", el("correlation plot viewer"), "-0.01"));
      await session.step(73, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Show Pearson R empties the cells and halves their width, and moves no coefficient", async () => {
      await session.step(76, "Then the \"cell width\" reading of correlation plot viewer should be 40", () => readingIs(page, "cell width", el("correlation plot viewer"), 40));
      await session.step(77, "And the \"column width of HEIGHT\" reading of correlation plot viewer should be 40", () => readingIs(page, "column width of HEIGHT", el("correlation plot viewer"), 40));
      await session.step(78, "When user sets \"showPearsonR\" property of correlation plot viewer to \"false\"", () => setProperty(page, "showPearsonR", el("correlation plot viewer"), "false"));
      await session.step(79, "Then the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), ""));
      await session.step(80, "And the \"text of cell WEIGHT x HEIGHT\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell WEIGHT x HEIGHT", el("correlation plot viewer"), ""));
      await session.step(81, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2349 and -0.2348", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.2349, -0.2348));
      await session.step(82, "And the \"cell width\" reading of correlation plot viewer should be 20", () => readingIs(page, "cell width", el("correlation plot viewer"), 20));
      await session.step(83, "And the \"column width of HEIGHT\" reading of correlation plot viewer should be 20", () => readingIs(page, "column width of HEIGHT", el("correlation plot viewer"), 20));
      await session.step(84, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(85, "When user sets \"showPearsonR\" property of correlation plot viewer to \"true\"", () => setProperty(page, "showPearsonR", el("correlation plot viewer"), "true"));
      await session.step(86, "Then the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(87, "And the \"cell width\" reading of correlation plot viewer should be 40", () => readingIs(page, "cell width", el("correlation plot viewer"), 40));
      await session.step(88, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Spearman puts a different coefficient in the same cells", async () => {
      await session.step(91, "When user sets \"correlationType\" property of correlation plot viewer to \"Spearman\"", () => setProperty(page, "correlationType", el("correlation plot viewer"), "Spearman"));
      await session.step(92, "Then the \"correlation type\" reading of correlation plot viewer should be \"Spearman\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Spearman"));
      await session.step(93, "And the correlation of \"HEIGHT\" and \"AGE\" of correlation plot viewer should match the Spearman coefficient of the table", () => correlationMatches(page, "HEIGHT", "AGE", el("correlation plot viewer"), "Spearman"));
      await session.step(94, "And the correlation of \"WEIGHT\" and \"HEIGHT\" of correlation plot viewer should match the Spearman coefficient of the table", () => correlationMatches(page, "WEIGHT", "HEIGHT", el("correlation plot viewer"), "Spearman"));
      await session.step(95, "And the \"correlation of HEIGHT and AGE\" reading of correlation plot viewer should be between -0.2400 and -0.2399", () => readingBetween(page, "correlation of HEIGHT and AGE", el("correlation plot viewer"), -0.24, -0.2399));
      await session.step(96, "And the \"correlation of WEIGHT and AGE\" reading of correlation plot viewer should be between 0.0943 and 0.0944", () => readingBetween(page, "correlation of WEIGHT and AGE", el("correlation plot viewer"), 0.0943, 0.0944));
      await session.step(97, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4453 and 0.4454", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4453, 0.4454));
      await session.step(98, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.24\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.24"));
      await session.step(99, "And correlation plot viewer should have repainted", () => repainted(page, el("correlation plot viewer")));
      await session.step(100, "When user sets \"correlationType\" property of correlation plot viewer to \"Pearson\"", () => setProperty(page, "correlationType", el("correlation plot viewer"), "Pearson"));
      await session.step(101, "Then the \"correlation type\" reading of correlation plot viewer should be \"Pearson\"", () => readingReads(page, "correlation type", el("correlation plot viewer"), "Pearson"));
      await session.step(102, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(103, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Narrowing the axes re-tiles the matrix to their product", async () => {
      await session.step(106, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT"],["yColumnNames","AGE, HEIGHT"]]));
      await session.step(109, "Then the \"cells\" reading of correlation plot viewer should be 6", () => readingIs(page, "cells", el("correlation plot viewer"), 6));
      await session.step(110, "And the \"x columns\" reading of correlation plot viewer should be \"AGE, HEIGHT, WEIGHT\"", () => readingReads(page, "x columns", el("correlation plot viewer"), "AGE, HEIGHT, WEIGHT"));
      await session.step(111, "And the \"y columns\" reading of correlation plot viewer should be \"AGE, HEIGHT\"", () => readingReads(page, "y columns", el("correlation plot viewer"), "AGE, HEIGHT"));
      await session.step(112, "And the \"columns shown\" reading of correlation plot viewer should be 5", () => readingIs(page, "columns shown", el("correlation plot viewer"), 5));
      await session.step(113, "And correlation plot viewer should have a \"cell WEIGHT x HEIGHT\" area", () => hasArea(page, el("correlation plot viewer"), "cell WEIGHT x HEIGHT"));
      await session.step(114, "And correlation plot viewer should not have a \"cell STARTED x AGE\" area", () => hasNoArea(page, el("correlation plot viewer"), "cell STARTED x AGE"));
      await session.step(115, "And correlation plot viewer should not have a \"cell AGE x WEIGHT\" area", () => hasNoArea(page, el("correlation plot viewer"), "cell AGE x WEIGHT"));
      await session.step(116, "And the correlation of \"WEIGHT\" and \"HEIGHT\" of correlation plot viewer should match the Pearson coefficient of the table", () => correlationMatches(page, "WEIGHT", "HEIGHT", el("correlation plot viewer"), "Pearson"));
      await session.step(117, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(120, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(121, "And correlation plot viewer should have a \"cell STARTED x AGE\" area", () => hasArea(page, el("correlation plot viewer"), "cell STARTED x AGE"));
      await session.step(122, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("A column with no variance has no coefficient and its cells stay blank", async () => {
      await session.step(125, "Given user adds a calculated column \"FLAT\" with formula \"0\"", () => addCalculated(page, "FLAT", "0"));
      await session.step(126, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED, FLAT"]]));
      await session.step(129, "Then the \"cells\" reading of correlation plot viewer should be 25", () => readingIs(page, "cells", el("correlation plot viewer"), 25));
      await session.step(130, "And the \"cell type of FLAT x AGE\" reading of correlation plot viewer should be \"correlation\"", () => readingReads(page, "cell type of FLAT x AGE", el("correlation plot viewer"), "correlation"));
      await session.step(131, "And the \"text of cell FLAT x AGE\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "text of cell FLAT x AGE", el("correlation plot viewer"), ""));
      await session.step(132, "And the \"text of cell HEIGHT x AGE\" reading of correlation plot viewer should be \"-0.23\"", () => readingReads(page, "text of cell HEIGHT x AGE", el("correlation plot viewer"), "-0.23"));
      await session.step(133, "And the \"error\" reading of correlation plot viewer should be \"\"", () => readingReads(page, "error", el("correlation plot viewer"), ""));
      await session.step(134, "And no errors should have been logged", () => noErrors(page));
      await session.step(135, "When user sets properties of correlation plot viewer:", () => setProperties(page, el("correlation plot viewer"), [["xColumnNames","AGE, HEIGHT, WEIGHT, STARTED"],["yColumnNames","AGE, HEIGHT, WEIGHT, STARTED"]]));
      await session.step(138, "And user removes \"FLAT\" column", () => removeColumn(page, "FLAT"));
      await session.step(139, "Then the \"cells\" reading of correlation plot viewer should be 16", () => readingIs(page, "cells", el("correlation plot viewer"), 16));
      await session.step(140, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("Cells six times apart in coefficient are painted the same full red (GROK — color_coding.dart:322-337)", async () => {
      await session.step(153, "Then the \"correlation of WEIGHT and AGE\" reading of correlation plot viewer should be between 0.0647 and 0.0649", () => readingBetween(page, "correlation of WEIGHT and AGE", el("correlation plot viewer"), 0.0647, 0.0649));
      await session.step(154, "And the \"correlation of WEIGHT and HEIGHT\" reading of correlation plot viewer should be between 0.4124 and 0.4125", () => readingBetween(page, "correlation of WEIGHT and HEIGHT", el("correlation plot viewer"), 0.4124, 0.4125));
      await session.step(155, "And the \"color of cell AGE x WEIGHT\" and \"color of cell HEIGHT x WEIGHT\" readings of correlation plot viewer should differ", () => readingsDiffer(page, "color of cell AGE x WEIGHT", "color of cell HEIGHT x WEIGHT", el("correlation plot viewer")));
      await session.step(156, "And the \"color of cell AGE x WEIGHT\" reading of correlation plot viewer should not be \"#ff0000\"", () => readingDoesNotRead(page, "color of cell AGE x WEIGHT", el("correlation plot viewer"), "#ff0000"));
      await session.step(157, "And no errors should have been logged", () => noErrors(page));
    }, {knownFailure: true});
    run.finish();
  });
});
