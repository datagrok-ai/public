// test.js

let factory = require('./lib');

factory().then((instance) => {
    for(let i = 1; i <= 10; i++)
       console.log(i + " <--> " + instance._sqr(i));
});