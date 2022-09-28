// test.js

// Simple Node app for testing lib.js

var factory = require('./newLib');

factory().then((instance) => {
    
    // 1-st function from lib1.c
    let n = 13;
    console.log(n + "^2 = " + instance._mySqr(n));

    // 2-st function from lib1.c
    let a = 1.11;
    let b = 9.99;
    console.log(a + " + " + b + " = " + instance._myAdd(a, b));

    // 3-rd function from lib1.c
    let x = 4;
    let y = 8;
    let z = 16;
    console.log(x + " * " + y + " * " + z + " = " + instance._tripleProduct(x, y, z));

    // 1-st function from lib2.c
    x = 9;
    y = 13;
    console.log(x + " * " + y + " = " + instance._myProd(x, y));

    // 2-st function from lib2.c
    a = 1.001;
    console.log(a + "^3 = " + instance._myCube(a));

    // 3-rd function from lib1.c
    x = 2022;
    y = 1983;
    console.log(x + " - " + y + " = " + instance._myDif(x, y));
});