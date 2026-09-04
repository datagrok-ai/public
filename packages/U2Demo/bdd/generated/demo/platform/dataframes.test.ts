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
import {openDemoPage} from '../../../bindings/demo.js';
import {loggedIn} from '@datagrok-libraries/bdd/bindings/common/session';
import {check, clickOn, selectIn, shouldContainText, shouldHaveText, typeInto} from '@datagrok-libraries/bdd/bindings/common/steps';
import {openDataset} from '@datagrok-libraries/bdd/bindings/platform/steps';
import {ds, el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("Dataframe pickers", () => {
  const session = feature(test);
  test("Table, column and columns pickers over an open table", {tag: ["@demo", "@realizes:u2.dg.pickers"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user is logged in", () => loggedIn(page));
    await test.step("And user opens demog dataset", () => openDataset(page, ds("demog")));
    await test.step("And user opens the \"Dataframes\" demo page", () => openDemoPage(page, "Dataframes"));
    enter(page, "U2 Demo");
    await test.step("When user selects \"demog\" in \"Open table\" input", () => selectIn(page, "demog", el("\"Open table\" input")));
    await test.step("Then value of table readout should have text \"demog\"", () => shouldHaveText(page, el("value of table readout"), "demog"));
    await test.step("When user selects \"age\" in \"Demog column\" input", () => selectIn(page, "age", el("\"Demog column\" input")));
    await test.step("Then value of \"demog column\" readout should have text \"age\"", () => shouldHaveText(page, el("value of \"demog column\" readout"), "age"));
    await test.step("When user types \"hei\" into Column input", () => typeInto(page, "hei", el("Column input")));
    await test.step("And user selects \"height\" in Column input", () => selectIn(page, "height", el("Column input")));
    await test.step("Then value of column readout should have text \"height\"", () => shouldHaveText(page, el("value of column readout"), "height"));
    await test.step("When user checks demog checkbox in Tables input", () => check(page, el("demog checkbox in Tables input")));
    await test.step("Then value of tables readout should have text \"demog\"", () => shouldHaveText(page, el("value of tables readout"), "demog"));
    await test.step("When user clicks on editor of Columns input", () => clickOn(page, el("editor of Columns input")));
    await test.step("And user clicks on All link", () => clickOn(page, el("All link")));
    await test.step("And user clicks on OK button", () => clickOn(page, el("OK button")));
    await test.step("Then value of columns readout should contain text \"age, sex\"", () => shouldContainText(page, el("value of columns readout"), "age, sex"));
    await test.step("When user selects \"age\" in x input in Mapping input", () => selectIn(page, "age", el("x input in Mapping input")));
    await test.step("And user selects \"height\" in y input in Mapping input", () => selectIn(page, "height", el("y input in Mapping input")));
    await test.step("Then value of mapping readout should contain text \"\\\"x\\\":\\\"age\\\"\"", () => shouldContainText(page, el("value of mapping readout"), "\"x\":\"age\""));
  });
});
