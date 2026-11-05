/* eslint-disable max-len */
/* eslint-disable comma-spacing */
/* eslint-disable quotes */
/* ---
generated: features/demo/platform/dataframes.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.dg.pickers]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import '@datagrok-libraries/bdd/bindings/tiers/molecules/crux';
import {openDemoPage} from '../../../bindings/demo.js';
import {asWideAs, atLeastWidePx, atMostWidePx, leftAligned, notClipped} from '@datagrok-libraries/bdd/bindings/common/geometry';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Dataframe pickers", () => {
  const session = feature(test, "features/demo/platform/dataframes.feature", import.meta.url);
  test("Table, column and columns pickers over an open table", {tag: ["@demo", "@realizes:u2.dg.pickers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await session.step(7, "Given user is logged in", () => loggedIn(page));
    await session.step(8, "And user opens demog-1000 dataset", () => openDataset(page, ds("demog-1000")));
    await session.step(9, "And user opens the \"Dataframes\" demo page", () => openDemoPage(page, "Dataframes"));
    enter(page, "U2 Demo");
    await session.step(10, "Then editor of Column input should be left-aligned with editor of \"Open table\" input", () => leftAligned(page, el("editor of Column input"), el("editor of \"Open table\" input")));
    await session.step(11, "And editor of Columns input should be as wide as editor of Column input", () => asWideAs(page, el("editor of Columns input"), el("editor of Column input")));
    await session.step(12, "And editor of Aggregations input should be as wide as editor of Columns input", () => asWideAs(page, el("editor of Aggregations input"), el("editor of Columns input")));
    await session.step(13, "When user selects \"demog-1000\" in \"Open table\" input", () => selectIn(page, "demog-1000", el("\"Open table\" input")));
    await session.step(14, "Then value of table readout should have text \"demog-1000\"", () => shouldHaveText(page, el("value of table readout"), "demog-1000"));
    await session.step(15, "And editor of Columns input should be at least 100 px wide", () => atLeastWidePx(page, el("editor of Columns input"), 100));
    await session.step(16, "When user selects \"age\" in \"Demog column\" input", () => selectIn(page, "age", el("\"Demog column\" input")));
    await session.step(17, "Then value of \"demog column\" readout should have text \"age\"", () => shouldHaveText(page, el("value of \"demog column\" readout"), "age"));
    await session.step(18, "When user types \"hei\" into Column input", () => typeInto(page, "hei", el("Column input")));
    await session.step(19, "And user selects \"height\" in Column input", () => selectIn(page, "height", el("Column input")));
    await session.step(20, "Then value of column readout should have text \"height\"", () => shouldHaveText(page, el("value of column readout"), "height"));
    await session.step(21, "When user checks demog-1000 checkbox in Tables input", () => check(page, el("demog-1000 checkbox in Tables input")));
    await session.step(22, "Then value of tables readout should have text \"demog-1000\"", () => shouldHaveText(page, el("value of tables readout"), "demog-1000"));
    await session.step(23, "When user clicks on editor of Columns input", () => clickOn(page, el("editor of Columns input")));
    await session.step(24, "Then OK button should not be clipped", () => notClipped(page, el("OK button")));
    await session.step(25, "When user clicks on All link", () => clickOn(page, el("All link")));
    await session.step(26, "And user clicks on OK button", () => clickOn(page, el("OK button")));
    await session.step(27, "Then value of columns readout should contain text \"age, sex\"", () => shouldContainText(page, el("value of columns readout"), "age, sex"));
    await session.step(28, "When user selects \"age\" in x input in Mapping input", () => selectIn(page, "age", el("x input in Mapping input")));
    await session.step(29, "And user selects \"height\" in y input in Mapping input", () => selectIn(page, "height", el("y input in Mapping input")));
    await session.step(30, "Then value of mapping readout should contain text \"\\\"x\\\":\\\"age\\\"\"", () => shouldContainText(page, el("value of mapping readout"), "\"x\":\"age\""));
    await session.step(31, "And label of x input in Mapping input should be at most 20 px wide", () => atMostWidePx(page, el("label of x input in Mapping input"), 20));
  });
});
