/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: simpleTest
//input: int num
export async function simpleTest(num) {

  alert(func());

  await initModule();

  alert(num + "^4 = " + MyModule._mySqr(num) * MyModule._mySqr(num));

  /*var factory = require('./lib.js');

  factory().then((instance) => {
      for(let i = 1; i <= 10; i++)
         console.log(i + " <--> " + instance._sqr(i));
  });*/
} // simpleTest


//name: tripleProduct
//input: int a
//input: int b
//input: int c
export async function tripleProduct(a, b, c) {
  await initNewModule();
  alert( NewModule._tripleProduct(a, b, c) );
}