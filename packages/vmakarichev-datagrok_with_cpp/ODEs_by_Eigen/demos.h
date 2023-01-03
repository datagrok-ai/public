// demos.h

#ifndef DEMOS_H
#define DEMOS_H

// demo of solving stiff example
int stiffExample(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount);

// demo of solving stiff example
int stiffExampleRK(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount);

// demo of solving jnj
int jnjStiff(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount);

// demo of solving jnj
int jnjRK4(float t0, float t1, float step,
	float* result, int resultRowCount, int resultColCount);

#endif // !DEMOS_H

