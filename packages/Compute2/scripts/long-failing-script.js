//name: LongFailingScript
//language: javascript
//input: double a
//input: double b
//output: double res

const delay = (delayInms) => {
  return new Promise((resolve) => setTimeout(resolve, delayInms));
};

await delay(3000);

throw new Error('Test error');

res = 1.0;
