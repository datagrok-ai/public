const utils = require("./utils.js");
const tests = {};
tests["filter-visual-query-mssql-test"] = async (driver, vars, opts = {}) => {
  await driver.get((new URL(`/datasets?selenium=true&q=`, BASE_URL)).href);
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div/div/input`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div/div/input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div/div[2]/input`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div/div[2]/input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@id=\'signup-login-fields\']/div[2]/button`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@id=\'signup-login-fields\']/div[2]/button`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//label[@name=\'label-Connect-to-data...\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//label[@name=\'label-Connect-to-data...\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-expander-MS-SQL\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-expander-MS-SQL\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-expander-MS-SQL---northwind\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-expander-MS-SQL---northwind\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`(//div[@name=\'tree-expander-MS-SQL---northwind\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//div[@name=\'tree-expander-MS-SQL---northwind\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---Tables---Products\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---Tables---Products\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Visual-Query...\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Visual-Query...\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.id(`db-query`)), configuration.timeout);
  await expect(driver.findElements(By.id(`db-query`))).resolves.toBePresent();
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-add-Rows\']/i[@name=\'icon-plus\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-add-Rows\']/i[@name=\'icon-plus\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-ProductName\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-ProductName\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-add-Measures\']/i[@name=\'icon-plus\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-add-Measures\']/i[@name=\'icon-plus\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Aggregation-(avg)\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Aggregation-(avg)\']`)).then(element => {
    return driver.actions({
      bridge: true
    }).move({
      origin: element
    }).perform();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Aggregation-(avg)---sum\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Aggregation-(avg)---sum\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-add-Measures\']/i[@name=\'icon-plus\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-add-Measures\']/i[@name=\'icon-plus\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-UnitsInStock\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-UnitsInStock\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-add-Filters\']/i[@name=\'icon-plus\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-add-Filters\']/i[@name=\'icon-plus\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-UnitPrice\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-UnitPrice\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.xpath(`(//input[@type=\'text\'])[3]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//input[@type=\'text\'])[3]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`(//input[@type=\'text\'])[3]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//input[@type=\'text\'])[3]`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`>15`);
    });
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-add-Filters\']/i[@name=\'icon-plus\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-add-Filters\']/i[@name=\'icon-plus\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-UnitPrice\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-UnitPrice\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`.d4-tag:nth-child(2) > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-tag:nth-child(2) > input`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`.d4-tag:nth-child(2) > input`)), configuration.timeout);
  await driver.findElement(By.css(`.d4-tag:nth-child(2) > input`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`<35`);
    });
  });
  await driver.sleep(1500);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.id(`editor`)), configuration.timeout);
  await driver.findElement(By.id(`editor`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.id(`editor`)), configuration.timeout);
  await driver.findElement(By.id(`editor`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`visual-query-filters-test`);
    });
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`button-SAVE`)), configuration.timeout);
  await driver.findElement(By.name(`button-SAVE`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-File\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-File\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-File---Close-All\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-File---Close-All\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`)).then(element => {
    return driver.actions({
      bridge: true
    }).doubleClick(element).perform();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.id(`balloon: info`)), configuration.timeout);
  await driver.findElement(By.id(`balloon: info`)).then(element => {
    return element.click();
  });
  await driver.sleep(3000);
  await driver.wait(until.elementLocated(By.id(`balloon: info`)), configuration.timeout);
  await driver.findElement(By.id(`balloon: info`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'view-handle: Products\']`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//div[@name=\'view-handle: Products\']`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'view-handle: Products\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'view-handle: Products\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.xpath(`(//div[@name=\'div-section--Columns\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//div[@name=\'div-section--Columns\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.name(`span-ProductName-(34-categories)`)), configuration.timeout);
  await expect(driver.findElements(By.name(`span-ProductName-(34-categories)`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.name(`span-sum(UnitsInStock)-(min-0,-max-123)`)), configuration.timeout);
  await expect(driver.findElements(By.name(`span-sum(UnitsInStock)-(min-0,-max-123)`))).resolves.toBePresent();
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'view-handle: Welcome\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'view-handle: Welcome\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-MS-SQL---northwind---visual-query-filters-test\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`(//div[@name=\'div-section--Actions\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//div[@name=\'div-section--Actions\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//label[contains(.,\'Delete\')]`)), configuration.timeout);
  await driver.findElement(By.xpath(`//label[contains(.,\'Delete\')]`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`button-YES`)), configuration.timeout);
  await driver.findElement(By.name(`button-YES`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
}
module.exports = tests;