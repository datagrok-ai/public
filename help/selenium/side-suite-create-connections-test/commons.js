const utils = require("./utils.js");
const tests = {};
tests["create-connections"] = async (driver, vars, opts = {}) => {
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
  await driver.wait(until.elementLocated(By.xpath(`//label[@name=\'label-Connect-to-data...\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//label[@name=\'label-Connect-to-data...\']`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"tree-PostgreSQL\"] > div.d4-tree-view-group-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"tree-PostgreSQL\"] > div.d4-tree-view-group-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Name`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Name`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`northwind_by_Selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Server`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Server`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`localhost`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Db`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Db`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`northwind`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Login`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Login`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`postgres`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Password`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Password`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`postgres`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`button-TEST`)), configuration.timeout);
  await driver.findElement(By.name(`button-TEST`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await expect(driver.findElement(By.css(`div.d4-balloon-content`))).resolves.toHaveText(`\"northwind_by_Selenium\": connected successfully`);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await driver.findElement(By.css(`div.d4-balloon-content`)).then(element => {
    return element.click();
  });
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`button-OK`)), configuration.timeout);
  await driver.findElement(By.name(`button-OK`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`div[name=\"tree-Sparql\"] > div.d4-tree-view-group-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"tree-Sparql\"] > div.d4-tree-view-group-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Name`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Name`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Name`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Name`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`rhea_by_selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Endpoint`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Endpoint`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`https://sparql.rhea-db.org/sparql`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Requires-Server`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Requires-Server`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`button-TEST`)), configuration.timeout);
  await driver.findElement(By.name(`button-TEST`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await expect(driver.findElement(By.css(`div.d4-balloon-content`))).resolves.toHaveText(`\"rhea_by_selenium\": connected successfully`);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await driver.findElement(By.css(`div.d4-balloon-content`)).then(element => {
    return element.click();
  });
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.name(`Edit-Connection---Prefixes`)), configuration.timeout);
  await driver.findElement(By.name(`Edit-Connection---Prefixes`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`Edit-Connection---Prefixes`)), configuration.timeout);
  await driver.findElement(By.name(`Edit-Connection---Prefixes`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`PREFIX ch:<http://purl.obolibrary.org/obo/>`);
    });
  });
  await driver.sleep(2000);
  await driver.wait(until.elementLocated(By.name(`button-OK`)), configuration.timeout);
  await driver.findElement(By.name(`button-OK`)).then(element => {
    return element.click();
  });
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`div[name=\"tree-PostgresNet\"] > div.d4-tree-view-group-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"tree-PostgresNet\"] > div.d4-tree-view-group-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Add-connection...\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Name`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Name`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Name`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Name`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`posgres_net_by_selenium`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Server`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Server`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`localhost`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Db`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Db`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`northwind`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Port`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Port`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`5432`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Login`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Login`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`postgres`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`input-Edit-Connection---Password`)), configuration.timeout);
  await driver.findElement(By.name(`input-Edit-Connection---Password`)).then(element => {
    return element.clear().then(() => {
      return element.sendKeys(`postgres`);
    });
  });
  await driver.wait(until.elementLocated(By.name(`button-TEST`)), configuration.timeout);
  await driver.findElement(By.name(`button-TEST`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await expect(driver.findElement(By.css(`div.d4-balloon-content`))).resolves.toHaveText(`\"posgres_net_by_selenium\": connected successfully`);
  await driver.wait(until.elementLocated(By.css(`div.d4-balloon-content`)), configuration.timeout);
  await driver.findElement(By.css(`div.d4-balloon-content`)).then(element => {
    return element.click();
  });
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`button-OK`)), configuration.timeout);
  await driver.findElement(By.name(`button-OK`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Admin\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Admin\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Admin---Data-Connections\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Admin---Data-Connections\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.xpath(`//div[@name=\'div-section--Filters\']`)), configuration.timeout);
  await driver.findElement(By.xpath(`//div[@name=\'div-section--Filters\']`)).then(element => {
    return element.click();
  });
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.css(`label[name=\"label-Created-by-me\"]`)), configuration.timeout);
  await driver.findElement(By.css(`label[name=\"label-Created-by-me\"]`)).then(element => {
    return element.click();
  });
  await driver.sleep(1000);
  await driver.wait(until.elementLocated(By.name(`span-NorthwindBySelenium`)), configuration.timeout);
  await expect(driver.findElements(By.name(`span-NorthwindBySelenium`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.name(`span-PosgresNetBySelenium`)), configuration.timeout);
  await expect(driver.findElements(By.name(`span-PosgresNetBySelenium`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.name(`span-RheaBySelenium`)), configuration.timeout);
  await expect(driver.findElements(By.name(`span-RheaBySelenium`))).resolves.toBePresent();
  await driver.wait(until.elementLocated(By.name(`span-NorthwindBySelenium`)), configuration.timeout);
  await driver.findElement(By.name(`span-NorthwindBySelenium`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`button-YES`)), configuration.timeout);
  await driver.findElement(By.name(`button-YES`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.name(`span-PosgresNetBySelenium`)), configuration.timeout);
  await driver.findElement(By.name(`span-PosgresNetBySelenium`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`button-YES`)), configuration.timeout);
  await driver.findElement(By.name(`button-YES`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await driver.wait(until.elementLocated(By.name(`span-RheaBySelenium`)), configuration.timeout);
  await driver.findElement(By.name(`span-RheaBySelenium`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`i[name=\"icon-context-arrow-down\"]`)), configuration.timeout);
  await driver.findElement(By.css(`i[name=\"icon-context-arrow-down\"]`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)), configuration.timeout);
  await driver.findElement(By.css(`div[name=\"div-Delete\"] > div.d4-menu-item-label`)).then(element => {
    return element.click();
  });
  await driver.wait(until.elementLocated(By.name(`button-YES`)), configuration.timeout);
  await driver.findElement(By.name(`button-YES`)).then(element => {
    return element.click();
  });
  await driver.sleep(500);
  await expect(driver.findElements(By.css(`i[name=\"icon-exclamation-circle\"]`))).resolves.not.toBePresent();
  await driver.sleep(500);
  await driver.get((new URL(`/connections`, BASE_URL)).href);
  await driver.sleep(1500);
  await driver.wait(until.elementLocated(By.xpath(`//label[contains(.,\'Created by me\')]`)), configuration.timeout);
  await driver.findElement(By.xpath(`//label[contains(.,\'Created by me\')]`)).then(element => {
    return element.click();
  });
  await driver.sleep(3000);
  await expect(driver.findElements(By.name(`span-NorthwindBySelenium`))).resolves.not.toBePresent();
  await expect(driver.findElements(By.name(`span-PosgresNetBySelenium`))).resolves.not.toBePresent();
  await expect(driver.findElements(By.name(`span-RheaBySelenium`))).resolves.not.toBePresent();
}
module.exports = tests;