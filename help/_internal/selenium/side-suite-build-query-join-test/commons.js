const utils = require("./utils.js");
const tests = {};
tests["build-query-join-test"] = async (driver, vars, opts = {}) => {
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
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-expander-PostgreSQL\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-expander-PostgreSQL\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-expander-PostgreSQL---datagrok\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-expander-PostgreSQL---datagrok\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`(//div[@name=\'tree-expander-PostgreSQL---datagrok\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//div[@name=\'tree-expander-PostgreSQL---datagrok\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'tree-PostgreSQL---datagrok---Tables---chats\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'tree-PostgreSQL---datagrok---Tables---chats\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-context-arrow-down\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Build-Query...\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Build-Query...\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.id(`query-builder`)), configuration.timeout);
  await expect(driver.findElements(By.id(`query-builder`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//input[@id=\'editor\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//input[@id=\'editor\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//input[@id=\'editor\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//input[@id=\'editor\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`(//input[@id=\'editor\'])[9]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//input[@id=\'editor\'])[9]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`(//input[@id=\'editor\'])[9]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//input[@id=\'editor\'])[9]`)).then(element => {
    return element.click();
  });
  await driver.sleep(5000);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`(//i[@name=\'icon-chevron-down\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//i[@name=\'icon-chevron-down\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Save-as-query\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Save-as-query\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//span[@name=\'span-view-name\']`)), configuration.timeout);
  await expect(driver.findElement(By.xpath(`//span[@name=\'span-view-name\']`))).resolves.toHaveText(`chats`);
  await driver.wait(until.elementLocated(By.xpath(`//i[@name=\'icon-play\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//i[@name=\'icon-play\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`(//i[@name=\'icon-chevron-down\'])[2]`)), configuration.timeout);
  await driver.findElement(By.xpath(`(//i[@name=\'icon-chevron-down\'])[2]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-Add-results-to-workspace\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-Add-results-to-workspace\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'view-handle: result\']`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//div[@name=\'view-handle: result\']`))).resolves.toBePresent();
  await expect(driver.findElements(By.xpath(`//i[@name=\'icon-exclamation-circle\']`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.xpath(`//span[@name=\'span-view-name\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//span[@name=\'span-view-name\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-section--General\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-section--General\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.xpath(`//span[contains(.,\'13\')]`)), configuration.timeout);
  await expect(driver.findElements(By.xpath(`//span[contains(.,\'13\')]`))).resolves.toBePresent();
}
module.exports = tests;