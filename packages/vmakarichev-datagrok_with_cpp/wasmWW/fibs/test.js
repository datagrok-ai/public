var factory = require('./fib.js');

factory().then((instance) => {
    console.log(instance._fib(10));    
});