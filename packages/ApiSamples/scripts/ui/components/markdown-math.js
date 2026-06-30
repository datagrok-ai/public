// LaTeX math in ui.markdown — use String.raw so JS doesn't consume backslash
// escapes (\b, \f, \n, \r, \t, \\) before the markdown parser sees them.
// Without String.raw, '\frac' becomes just 'rac', '\beta' becomes a backspace.

let document = String.raw`# Advanced pendulum

A damped pendulum is governed by

$$
\begin{aligned}
\frac{d\theta}{dt} &= \omega \\
\frac{d\omega}{dt} &= -\frac{g}{L}\sin(\theta) - \gamma\, \omega
\end{aligned}
$$

over $t \in [0, 10]$, with Greek parameters $\alpha$, $\beta$, $\gamma$.

### Parameters

| Symbol | Meaning         | Default |
|--------|-----------------|---------|
| $L$    | Length          | 1.0     |
| $g$    | Gravity         | 9.81    |
| $\gamma$ | Damping       | 0.1     |

### Constants

Inline: $P_{1} = mgL$, $\pi \approx 3.14159$.
`;

grok.shell.newView('Markdown math example', [
  ui.markdown(document)
]);
