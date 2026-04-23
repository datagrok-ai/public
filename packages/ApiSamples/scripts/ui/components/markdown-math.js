// Render markdown with LaTeX math. Uses String.raw so backslashes reach the
// renderer intact (JS template literals process \b, \f, \r, \\ at parse time).

let document = String.raw`# Advanced

This model demonstrates a system of differential equations with inline and block math.

$$
\begin{aligned}
  \frac{dP_1}{dt} &= k_1 \cdot C_1 - k_2 \cdot P_1 \\
  \frac{dC_1}{dt} &= -k_1 \cdot C_1
\end{aligned}
$$

With $t \in [0, 10]$ and $\Delta t = 0.01$.

### Parameters

| Name | Symbol | Default |
|------|--------|---------|
| Rate constant | $k_1$ | 0.5 |
| Rate constant | $k_2$ | 0.1 |

### Constants

| Name | Symbol | Value |
|------|--------|-------|
| Initial precursor | $P_1$ | 0 |
| Initial compound | $C_1$ | 100 |

Inline Greek and subscripts: $\alpha$, $\beta$, $x_n$, $y_n$, $\sum_{i=1}^{n} x_i$.
`;

grok.shell.newView('Markdown math', [
  ui.markdown(document)
]);
