# Deriving the true confirmation probability $p_i$

This note derives the **true per-individual confirmation probability** $p_i$ used as
the oracle weight in the clone-weighted Kaplan–Meier estimator, and explains how it
is computed in `scripts/Simulation.R` (function `p_true_i()`).

The quantity $p_i$ is the probability that a **first positive** observed at test $i$
would be **confirmed** by a positive result at the next test $i+1$. It is exactly the
weight the cloning step needs: a censored first-positive is split into an *event*
clone with weight $p_i$ and a *censored* clone with weight $1 - p_i$.

---

## 1. The data-generating model

From `simulate_data()`:

- **True onset** time $T_0 \sim \text{Exp}(\lambda)$, with survival function
  $S(t) = P(T_0 > t) = e^{-\lambda t}$.
- **Tests** occur at times $O_1 < O_2 < \dots < O_K$. Gaps depend on the previous
  result (shorter after a positive — closer monitoring), but for this derivation we
  treat the relevant test times as given.
- **True status** at test $k$: $R_k = \mathbf{1}\{T_0 \le O_k\}$.
- **Observed result** at test $k$:
  - **Sensitivity = 1** — a truly positive person always tests positive:
    $P(RE_k = 1 \mid R_k = 1) = 1$.
  - **False positives** fire independently with rate $\varepsilon$ when truly
    negative: $P(RE_k = 1 \mid R_k = 0) = \varepsilon$.

So the only source of error is false positives among the truly negative.

---

## 2. What $p_i$ means

A **first positive at test $i$** is the event $F_i$: negative at every test before
$i$ and positive at test $i$ ($RE_1 = \dots = RE_{i-1} = 0,\ RE_i = 1$).
**Confirmation** is $RE_{i+1} = 1$. We want

$$
p_i \;=\; P\big(RE_{i+1} = 1 \,\big|\, F_i\big).
$$

A first positive can be **genuine** (the person truly developed the condition) or a
**false positive**, and the chance of confirmation is very different in the two
cases. The mix between them depends on how much *incremental prevalence* accrues
between the surrounding test times — which is why $p_i$ depends on timing.

> **Why earlier history drops out.** $F_i$ forces $T_0 > O_{i-1}$ (any earlier true
> onset would have made an earlier test positive, since sensitivity = 1). Given
> $T_0 > O_{i-1}$, the "no false positive at the earlier tests" requirement
> contributes a common factor $(1-\varepsilon)^{\,i-1}$ to every path below, so it
> cancels in the ratio. Only the structure around tests $i$ and $i+1$ matters.

---

## 3. Decomposition into paths

Define the **interval onset probabilities** (the "prevalence, evaluated at the actual
test times"):

$$
q_i = S(O_{i-1}) - S(O_i), \qquad
q_{i+1} = S(O_i) - S(O_{i+1}), \qquad
s_{i+1} = S(O_{i+1}),
$$

and note $s_i \equiv S(O_i) = q_{i+1} + s_{i+1}$.

There are exactly three disjoint ways to be a first positive at $i$ **and** be
positive at $i+1$. Each factor of $\varepsilon$ is one false positive the path
requires:

| Path | Where onset lands | Positive at $i$? | Positive at $i+1$? | Contribution |
|---|---|---|---|---|
| 1. genuine–genuine | $(O_{i-1}, O_i]$ → prob $q_i$ | yes, genuine | yes, genuine | $q_i$ |
| 2. false then genuine | $(O_i, O_{i+1}]$ → prob $q_{i+1}$ | false ($\times\varepsilon$) | yes, genuine | $\varepsilon\, q_{i+1}$ |
| 3. false then false | $> O_{i+1}$ → prob $s_{i+1}$ | false ($\times\varepsilon$) | false ($\times\varepsilon$) | $\varepsilon^2 s_{i+1}$ |

The denominator (all first positives) is genuine plus false-at-$i$:
$q_i + \varepsilon\, s_i$. Hence

$$
\boxed{\,p_i = \dfrac{q_i + \varepsilon\, q_{i+1} + \varepsilon^2 s_{i+1}}
{q_i + \varepsilon\, s_i}\,}
$$

Numerator and denominator differ only in the last term, so this also equals the
interpretable "one minus leakage" form

$$
p_i = 1 - \frac{\varepsilon(1-\varepsilon)\, s_{i+1}}{q_i + \varepsilon\,(q_{i+1}+s_{i+1})},
$$

where the leakage term is the single unconfirmed path: false positive now → still no
onset by the next test → no second false positive.

---

## 4. Memorylessness gives a gap-only closed form

The exponential is **memoryless**: $S(s+t) = S(s)\,S(t)$, equivalently
$P(T_0 > s+t \mid T_0 > s) = P(T_0 > t)$. Once a person is condition-free at a test,
the wait until onset restarts as $\text{Exp}(\lambda)$ — independent of how long they
have already waited.

Write the two surrounding gaps

$$
g_1 = O_i - O_{i-1} \quad (\text{gap into the first positive}), \qquad
g_2 = O_{i+1} - O_i \quad (\text{confirming gap}).
$$

Factoring out the common "reach the window" survival $S(O_{i-1}) = e^{-\lambda O_{i-1}}$:

$$
q_i = S(O_{i-1})\,a, \quad
q_{i+1} = S(O_{i-1})\,b\,c, \quad
s_{i+1} = S(O_{i-1})\,b\,d,
$$

with

$$
a = 1 - e^{-\lambda g_1}\ (\text{onset in gap 1}), \quad
b = e^{-\lambda g_1}\ (\text{survive gap 1}), \quad
c = 1 - e^{-\lambda g_2}, \quad
d = e^{-\lambda g_2}.
$$

The factor $S(O_{i-1})$ appears in **every** term of both numerator and denominator
and cancels, leaving a formula that depends only on the two gaps and $\varepsilon$
(not on the absolute test times, the index $i$, or the earlier history):

$$
\boxed{\,p_i = \dfrac{a + \varepsilon\, b\, c + \varepsilon^2\, b\, d}{a + \varepsilon\, b}\,}
$$

This position-independence is special to the exponential. For a non-memoryless onset
(e.g. Weibull with shape $\ne 1$), $S(O_{i-1})$ does not factor out and $p_i$ would
genuinely depend on *where* on the timeline the tests fall, requiring the absolute
test times.

---

## 5. Sanity checks

- **$\varepsilon \to 0$** (no false positives): leakage $\to 0$, so $p_i \to 1$.
  Every first positive is genuine and always confirmed. ✓
- **$\varepsilon \to 1$** (test always positive): $(1-\varepsilon)\to 0$, so
  $p_i \to 1$. Always confirmed. ✓
- **Low-prevalence region** ($q_i, q_{i+1} \to 0$, i.e. $a, c \to 0$): $p_i \to
  \varepsilon$. A spurious positive where no onset occurs is confirmed only by a
  second false positive, probability $\varepsilon$. ✓

---

## 6. How it is computed in the simulation

`p_true_i(O, interval, lambda, gap_next, eps)` returns a length-$n$ vector applying
the boxed gap form to each censored first-positive **at whatever interval $i$ it
occurs** (test 1 included). `interval` is `compute_double()$amb`, the per-individual
first-positive interval (`NA` where there is no censored first-positive). For an
ambiguous individual at interval $i$:

- $g_1 = O_i - O_{i-1}$ — taken from the **observed** test times, so it is exact and
  varies across individuals. We set $O_0 := 0$ (the time origin), so a
  **first positive at test 1** uses $g_1 = O_1$.
- $g_2 = $ `gap_next` — the **expected** confirming gap. The confirming test is
  unobserved precisely because the individual is a censored first-positive, so we
  plug in the after-positive gap mean (`gap_pos`), the design value for a test that
  follows a positive.

```r
p_vec <- p_true_i(d$O, doub$amb, lambda, gap_pos, eps)
cw_dt <- build_clone_weighted_data(d$O, doub$dt, doub$amb, p_vec)
```

`build_clone_weighted_data()` takes the confirmation probability `p` as an argument.
It accepts a **scalar**, a **length-$n$ vector** (one per individual, as returned by
`p_true_i()`), or a **length-(#ambiguous) vector** — so collaborators can drop in
their own weights without touching the cloning logic. The event clone is dated at
each individual's own first-positive time $O_i$ (interval-aware).

> **Why this matters (the test-1 case).** When the after-positive gap is large enough
> that the confirming test for an *early* positive crosses the cutoff (e.g.
> $O_1 + $ `gap_pos` $\approx$ `cutoff`), genuine early events become censored
> first-positives at **test 1**. An earlier, `k = 2`-only version of the cloning left
> these uncorrected, so the clone-weighted estimator inherited the double-positive
> bias at large gaps. Cloning all intervals removes that bias.

---

## 7. The data-driven alternative, and its caveats

`p_hat_empirical(O, RE, cutoff)` returns the simple plug-in estimate

$$
\hat p = \frac{\#\{\text{positive at tests 1 and 2},\ O_2 < \text{cutoff}\}}
{\#\{\text{positive at test 1},\ O_2 < \text{cutoff}\}},
$$

a single global number that estimates the *marginal* version of $p_i$ from the
test-1 → test-2 transition. It is provided for users who must estimate $p$ from data,
but note its limitations relative to the oracle $p_i$:

1. **Single global value** — it collapses the interval- and individual-specific
   $p_i$ to one number, assuming confirmation behaviour is homogeneous.
2. **Exchangeability (MAR-style)** — it is estimated from first-positives whose
   confirmation *is* observed, then applied to censored first-positives whose
   confirmation is *not* observed; validity needs the two groups to share the same
   confirmation probability.
3. **Horizon selection** — conditioning on $O_2 < \text{cutoff}$ drops late second
   tests, which biases $\hat p$ if gap correlates with confirmation.
4. **Small-sample instability** — the ratio is $0/0$ when no one is positive at test
   1 within the horizon (possible at small $\varepsilon$ or small $n$).

Using `p_true_i()` in the main pipeline lets the simulation isolate the clone-
weighting **mechanism** from the **estimation** of $p$; swapping in
`p_hat_empirical()` (or a modelled estimate) then shows the additional cost of having
to estimate the weight.
