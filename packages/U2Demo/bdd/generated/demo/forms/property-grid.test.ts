/* ---
generated: features/demo/forms/property-grid.feature
generator: @datagrok-libraries/bdd — do not edit; run `grok-bdd compile` to regenerate
sub_features_covered: [u2.property-grid]
--- */
import {test} from '@playwright/test';
import '../../../bindings/elements.js';
import '@datagrok-libraries/bdd/bindings/common/kinds';
import '@datagrok-libraries/bdd/bindings/common/parameter-types';
import '@datagrok-libraries/bdd/bindings/platform/datasets';
import '@datagrok-libraries/bdd/bindings/platform/elements';
import {openDemoPage} from '../../../bindings/demo.js';
import {collapse, enterInto, expand, selectIn, shouldBe, shouldContainText, shouldHaveText, typeInto, uncheck} from '@datagrok-libraries/bdd/bindings/common/steps';
import {el, enter, feature} from '@datagrok-libraries/bdd/runtime';

test.describe("The property grid", () => {
  const session = feature(test);
  test("Editing a property replaces the value record", {tag: ["@demo", "@realizes:u2.property-grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Property grid\" demo page", () => openDemoPage(page, "Property grid"));
    enter(page, "U2 Demo");
    await test.step("Then value of \"last change\" readout should have text \"(none)\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "(none)"));
    await test.step("When user types \"Heatmap\" into title property", () => typeInto(page, "Heatmap", el("title property")));
    await test.step("Then value of \"last change\" readout should have text \"title → Heatmap\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "title → Heatmap"));
    await test.step("And value of values readout should contain text \"\\\"title\\\":\\\"Heatmap\\\"\"", () => shouldContainText(page, el("value of values readout"), "\"title\":\"Heatmap\""));
    await test.step("When user unchecks showLegend property", () => uncheck(page, el("showLegend property")));
    await test.step("Then value of \"last change\" readout should have text \"showLegend → false\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "showLegend → false"));
    await test.step("When user selects \"top\" in position property", () => selectIn(page, "top", el("position property")));
    await test.step("Then value of \"last change\" readout should have text \"position → top\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "position → top"));
    await test.step("When user enters \"0.5\" into opacity property", () => enterInto(page, "0.5", el("opacity property")));
    await test.step("Then value of \"last change\" readout should have text \"opacity → 0.5\"", () => shouldHaveText(page, el("value of \"last change\" readout"), "opacity → 0.5"));
  });
  test("Categories collapse", {tag: ["@demo", "@realizes:u2.property-grid"]}, async ({browser}) => {
    const page = await session.page(browser);
    await test.step("Given user opens the \"Property grid\" demo page", () => openDemoPage(page, "Property grid"));
    enter(page, "U2 Demo");
    await test.step("Then Layout category should be expanded", () => shouldBe(page, el("Layout category"), "expanded"));
    await test.step("And width property should be visible", () => shouldBe(page, el("width property"), "visible"));
    await test.step("When user collapses Layout category", () => collapse(page, el("Layout category")));
    await test.step("Then Layout category should be collapsed", () => shouldBe(page, el("Layout category"), "collapsed"));
    await test.step("And width property should be hidden", () => shouldBe(page, el("width property"), "hidden"));
    await test.step("When user expands Layout category", () => expand(page, el("Layout category")));
    await test.step("Then width property should be visible", () => shouldBe(page, el("width property"), "visible"));
  });
});
