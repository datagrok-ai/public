# FAEinWebWorker

FAEinWebWorker is a [package](https://datagrok.ai/help/develop/develop#packages) for the [Datagrok](https://datagrok.ai) platform.

This is a demo package. It shows the difference between wasm-computations with and without webworkers.

Solution of the problem FAE using explicit method (Runge-Kutta-Cash-Karp) is provided.
There are two functions:
 
  - solveFAEexplicit(...) - all computations are performed in the main stream,
                            also, ReachFunctionView is applied;
 
  - solveFAEexplicitWebWorker(...) - all computations are performed in webworker,
                                     also, ReachFunctionView is NOT applied!

Currently, runtime-system for C++/wasm-functions call should be extended!