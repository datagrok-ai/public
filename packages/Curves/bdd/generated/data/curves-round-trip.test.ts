/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/data/curves-round-trip.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
--- */
import {test} from '@playwright/test';
import '../../bindings/curves.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {clickOn, downloadThrough, lastDownloadCount, uploadDownloaded, watchDownloads} from '@datagrok-libraries/bdd/bindings/common/steps';
import {columnSemType, columnsExactly, everyValueContains} from '@datagrok-libraries/bdd/bindings/platform/columns';
import {rowCount} from '@datagrok-libraries/bdd/bindings/platform/data';
import {autostartsCompleted, closeAllViews, noProjectOnServer, openDataset, openProjectWithTable, packageInstalled, saveAsProject, simpleModeOff} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {noErrors} from '@datagrok-libraries/bdd/bindings/tiers/viewers/steps';
import {ds, el, feature, journey} from '@datagrok-libraries/bdd/runtime';

test.describe("Curve columns come back as curves from an exported CSV and from a project", () => {
  const session = feature(test, "features/data/curves-round-trip.feature", import.meta.url);
  test("Curve columns come back as curves from an exported CSV and from a project", {tag: ["@journey"]}, async ({browser}) => {
    const page = await session.page(browser);
    const run = journey(test, 2, page);
    await session.step(9, "Given user is logged in", () => loggedIn(page));
    await session.step(10, "And simple mode is off", () => simpleModeOff(page));
    await session.step(11, "And the package autostarts have completed", () => autostartsCompleted(page));
    await session.step(12, "And the \"Curves\" package is installed", () => packageInstalled(page, "Curves"));
    await session.step(13, "And no project named \"CurvesRoundTrip-{run}\" is on the server", () => noProjectOnServer(page, session.text("CurvesRoundTrip-{run}")));
    await run.scenario("The exported CSV opens back with its 30 rows, six columns and five curve columns", async () => {
      await session.step(16, "Given user opens curves dataset", () => openDataset(page, ds("curves")));
      await session.step(17, "And user watches downloads", () => watchDownloads(page));
      await session.step(18, "When user clicks on Export icon in toolbar", () => clickOn(page, el("Export icon in toolbar")));
      await session.step(19, "And user downloads a file through \"As CSV\" text in toolbar", () => downloadThrough(page, el("\"As CSV\" text in toolbar")));
      await session.step(20, "Then the downloaded file should contain 150 occurrences of \"Dose-Response curves\"", () => lastDownloadCount(page, 150, "Dose-Response curves"));
      await session.step(21, "When user closes all views", () => closeAllViews(page));
      await session.step(22, "And user clicks on \"Browse\" view", () => clickOn(page, el("\"Browse\" view")));
      await session.step(23, "And user uploads the downloaded file through \"Open local file\" icon inside browse toolbar", () => uploadDownloaded(page, el("\"Open local file\" icon inside browse toolbar")));
      await session.step(24, "Then the table should have 30 rows", () => rowCount(page, 30));
      await session.step(25, "And the table should have the columns \"smiles, multiple prefit, multiple styled series, styled proprtional with IC50, log min/max candlesticks, custom fit function\"", () => columnsExactly(page, "smiles, multiple prefit, multiple styled series, styled proprtional with IC50, log min/max candlesticks, custom fit function"));
      await session.step(26, "And \"multiple prefit\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple prefit", "fit"));
      await session.step(27, "And \"multiple styled series\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple styled series", "fit"));
      await session.step(28, "And \"styled proprtional with IC50\" column should have semantic type \"fit\"", () => columnSemType(page, "styled proprtional with IC50", "fit"));
      await session.step(29, "And \"log min/max candlesticks\" column should have semantic type \"fit\"", () => columnSemType(page, "log min/max candlesticks", "fit"));
      await session.step(30, "And \"custom fit function\" column should have semantic type \"fit\"", () => columnSemType(page, "custom fit function", "fit"));
      await session.step(31, "And every value of \"custom fit function\" column should contain \"series\"", () => everyValueContains(page, "custom fit function", "series"));
      await session.step(32, "And no errors should have been logged", () => noErrors(page));
    });
    await run.scenario("The curve columns are curves again in the reopened project", async () => {
      await session.step(35, "When user saves the current view as project \"CurvesRoundTrip-{run}\"", () => saveAsProject(page, session.text("CurvesRoundTrip-{run}")));
      await session.step(36, "And user closes all views", () => closeAllViews(page));
      await session.step(37, "And user opens the \"CurvesRoundTrip-{run}\" project and waits for its table", () => openProjectWithTable(page, session.text("CurvesRoundTrip-{run}")));
      await session.step(38, "Then the table should have 30 rows", () => rowCount(page, 30));
      await session.step(39, "And \"multiple prefit\" column should have semantic type \"fit\"", () => columnSemType(page, "multiple prefit", "fit"));
      await session.step(40, "And \"custom fit function\" column should have semantic type \"fit\"", () => columnSemType(page, "custom fit function", "fit"));
      await session.step(41, "And no errors should have been logged", () => noErrors(page));
    });
    run.finish();
  });
});
