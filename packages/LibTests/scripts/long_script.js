//name: LongScript
//language: javascript
//input: double a
//input: double b
//output: double res

const delay = (delayInms) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};
  
await delay(5000);

res = a + b;
