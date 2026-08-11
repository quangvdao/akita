# Spec: coefficient-carrier openings and subring fold challenges

| Field | Value |
|---|---|
| Author(s) | Quang Dao |
| Created | 2026-08-10 |
| Status | proposed |
| PR | |
| Supersedes | The assumption that every extension-field opening first uses extension-opening reduction |
| Superseded-by | |
| Book-chapter | book/src/how/proving/root-fold-ring-switch.md |

The key words **MUST**, **MUST NOT**, **REQUIRED**, **SHOULD**, **SHOULD NOT**,
and **MAY** in this document are to be interpreted as described in
[BCP 14](https://www.rfc-editor.org/info/bcp14) when, and only when, they appear
in all capitals.

## Decision

Akita will support a direct opening mode for base-field committed tables opened
at extension-field points. The mode keeps one coefficient axis explicit as a
smaller cyclotomic **carrier ring** and contracts the other coefficient axes
directly over the extension field. Fold challenges live in that smaller ring
and embed sparsely into the ambient A ring. The challenge dimension and the
ambient A-ring dimension are independent schedule choices, subject to the
divisibility condition below.

For extension-field presets, generated schedules MUST use this direct mode at
absolute fold levels 0 and 1. The planner may change the A dimension or the
carrier dimension to make the mode feasible. It MUST NOT silently retain
extension-opening reduction (EOR) at either level. A catalog row with no
feasible direct candidate is unsupported until its geometry or audited
challenge family changes.

The direct mode MUST be representable at every later fold, including the
terminal fold. Nonterminal folds keep its coefficient planes in the ordinary
recursively committed witness. A transparent direct terminal reveals those
planes. The planner MAY compare that cost with the current EOR/Hachi terminal
while the coefficient-native terminal-compression problem below remains open.

The planner MUST NOT force later folds onto a small uniform-ring suffix. At
every fold it retains the supported A-dimension ladder allowed by the incoming
geometry, even when B and D use smaller rings. In particular, increasing `d_A`
does not increase a direct partial when `s` is unchanged. The planner SHOULD
normally keep the smallest admitted `s` while increasing `d_A` when the larger
ring reduces the secure A rank without increasing the exact output dimension.

This specification is stacked on [PR #383](https://github.com/LayerZero-Labs/akita/pull/383),
at commit `38ab14924e6539875abab97b706970baaa973ce3`. PR #383's dyadic B
slicing is part of the planning baseline. PR #368 is not a dependency.

## Why this change

Akita currently commits base-field coefficients but uses an extension field for
opening points and sum-check challenges. EOR converts that extension-field
opening into a form accepted by the ring relation. It is correct, but it adds a
degree-two sum-check and transcript-bound partial evaluations at every
extension-valued recursive opening.

The conversion is avoidable. A coefficient table can instead be viewed as a
polynomial over a smaller carrier ring whose coefficients already lie in the
extension field. This viewpoint has three useful consequences:

1. The original extension-field point can be used directly, so L0 and L1 need
   no EOR proof.
2. Each partial opening and its consistency quotient can use fewer base-field
   coordinates than a full ambient-ring element.
3. The current A relation, B commitment, setup matrices, and most NTT caches can
   remain over the ambient base-field ring.

The change is not unconditionally cheaper. A smaller carrier requires a
different sparse challenge family. Meeting the same entropy target can increase
the challenge's norm, which can enlarge the folded witness, the secure A rank,
and the `t` part of the next witness. Near the recursive tail, a direct
small-field opening may also require a larger A ring. This does not by itself
make the candidate worse: a larger ring can need proportionally fewer A rows.
The planner must compare exact field-coordinate counts and the complete suffix,
not ring dimension alone.

## Scope and notation

The direct-mode equations below describe one polynomial and one commitment
group. Existing claim batching and multi-group row coefficients apply outside
these equations and do not change them.

| Symbol | Meaning |
|---|---|
| `K` | Base coefficient field `F_q` |
| `E` | Challenge and evaluation field, with extension degree `k = [E:K]` |
| `d_A` | Ambient A-ring dimension |
| `n_A` | Secure output rank of the A matrix |
| `s` | Carrier-ring dimension selected by the schedule |
| `h` | Native-mode packing gain `d_A / (k s)` |
| `R` | Ambient ring `K[X]/(X^{d_A}+1)` |
| `S` | Challenge carrier `K[Y]/(Y^s+1)` |
| `C` | Extension carrier `E[Y]/(Y^s+1)` |
| `beta_t` | Element `t` of Akita's fixed canonical `K`-basis of `E` |

Every admitted direct candidate satisfies

```text
d_A = k h s,
h >= 1,
d_A, k, h, and s are powers of two.
```

The embedding from the carrier into the ambient ring is

```text
S -> R,       Y -> X^(k h).
```

It preserves coefficient support and the coefficient `l1`, `l2`, and `linf`
norms:

```text
c(Y) = sum_(j < s) c_j Y^j
  maps to
c(X^(k h)) = sum_(j < s) c_j X^(k h j).
```

The implementation MUST use this canonical embedding. It MUST NOT search over
coefficient permutations or alternative carrier embeddings.

### Final architecture boundary

There is one coefficient-native opening relation. At a nonterminal fold, its
partial digits are part of the recursively opened next witness. At a
transparent terminal, its `P B k s` coefficients are revealed. Merely replacing
that reveal by an SIS image and local sum-checks is unsound, as explained
below: the final hidden-witness evaluation still needs an authenticated
opening.

The native relation MUST work at the root, recursive folds, and terminal fold.
The descriptor MUST bind `k`, `s`, the extension-coordinate order, and the
selected terminal path.

The current tensor EOR followed by Hachi `psi` packing and a trace check remains
the supported compressed terminal. It can be removed only after
generated-catalog evidence shows that every supported terminal can afford its
transparent coefficient-native payload.

### Why the tail is the hard case

The current audited sparse challenge ladder starts at dimension 64. A native
candidate using that smallest challenge therefore carries

```text
k s = 64       for fp128,
k s = 128      for fp64,
k s = 256      for fp32.
```

This is independent of `d_A`. Thus fp32 can keep `s = 64` and use `d_A = 256`,
while fp128 can keep `s = 64` with any admitted `d_A` in `{64, 128, 256}`. The
larger A ring is often harmless when its secure rank falls proportionally. For
example, A dimension 64 at rank 4 and A dimension 256 at rank 1 both output 256
base-field coordinates; if their remaining exact successor geometry is the
same, the dimension-256 candidate is preferred.

The PR #383 tables also show why the full ladder, rather than a rule that always
raises `d_A`, is necessary. A representative fp32 terminal has 128 A-input
columns at coefficient bound `2^20 - 1`. The audited Q32 table gives:

| A dimension | Secure rank | A image width |
|---:|---:|---:|
| 128 | 8 | 1,024 |
| 256 | 5 | 1,280 |
| 512 | 2 | 1,024 |

At dimension 256, rank 4 supports only 124 columns, so the fifth row makes that
candidate worse. Dimension 512 recovers the same A image width as dimension
128 and can carry the native fp32 `s = 64` challenge. The exact tail comparison
is then between transparent native partials and the current EOR/Hachi path. A
forced dimension-128 suffix cannot make that comparison.

The terminal comparison is therefore between the transparent dimension-64
coefficient carrier and the current EOR/Hachi path. Smaller challenge carriers
and hidden coefficient-native terminals are outside this specification.

## Terminal compression boundary

Suppose the prover decomposes the native planes into short digits `ehat`, binds

```text
v_D = D ehat,
```

and then runs range and relation sum-checks without revealing `ehat`. Those
sum-checks end at a random point `r_sc` and require the claimed value
`ehat(r_sc)`. The verifier cannot derive that value from `v_D`.

This is not repaired by including `D ehat = v_D` inside the same sum-check. A
cheating prover can choose the last value after seeing `r_sc` and make the
univariate transcript close without committing to one global low-degree table.
Stage 1 has the same final-oracle problem. The D image gives computational
binding *if two short preimages are already known*; it is not by itself a
polynomial-opening protocol.

The current nonterminal protocol supplies the missing step by putting `ehat`
inside the recursively committed next witness and opening that commitment at
the Stage-2 point. The current transparent terminal supplies it by revealing
the witness, from which the verifier computes the final evaluation directly.

A sound no-raw-`e` terminal therefore MUST provide all of the following:

1. bind the coefficient planes before the block-fold challenges;
2. prove the bound witness is short;
3. authenticate every final witness evaluation requested by the range and
   relation proofs; and
4. enforce native carrier consistency and the exact `E`-valued opening on that
   same authenticated witness.

Until item 3 has a concrete proof, the planner MUST NOT price the 512-byte D
image and local sum-check messages as a complete replacement for raw `e`. The
previous 1,744-byte estimate omitted this opening argument and is invalid.

The terminal paths currently safe to price are:

- transparent native planes at the admitted carrier dimension; and
- the current EOR/Hachi path.

## Current protocol

For each live claim/block pair, the current prover produces a partial opening
`e_i` as a full element of `R`, hence `d_A` base-field coordinates. At a
nonterminal fold, it gadget-decomposes those coordinates into `e_hat`, commits
the D image, and absorbs that payload before sampling the sparse fold
challenges `c_i`. The transparent terminal is the exception: it drops the D
rows and reveals the reduced `e_i` coefficients as raw field elements.

The two A-native relations are, schematically,

```text
sum_i c_i e_i = a z                         in R
[A G z_hat]_r = sum_i c_i [G t_hat_i]_r     in R, for every A row r.
```

Here `a` is the current ring opening multiplier and `G` denotes the applicable
gadget recomposition. The first equation is the consistency row. The remaining
equations are the A rows. Their polynomial representatives need quotient rows
for divisibility by `X^{d_A}+1`.

When `k > 1`, EOR first changes the opening claim and protocol point. The
evaluation-trace row then uses the ring-subfield trace/Galois construction to
connect the ring partials to the original scalar opening. The ring-switch
verifier evaluates every sparse challenge once at `X = alpha` and reuses that
value in the consistency and A contractions.

This specification changes four facts: a direct partial is not a full `R`
element; the consistency row uses the carrier modulus; the scalar row uses
direct coefficient weights; and the carrier and A relations require different
evaluations of the same challenge. At the terminal, the direct path is
transparent.

## Direct coefficient packing

### Canonical coefficient layout

Write one ambient ring element at block `i` and position `x` as

```text
F_(i,x)(X)
  = sum_(j < s) sum_(a < k h) f_(i,x,a,j) X^(a + k h j).
```

Thus `a` is the low coefficient index and `j` is the carrier index. The physical
ambient coefficient index is exactly

```text
a + k h j.
```

The opening point's coefficient variables are split in the same order:

```text
r_pack     has log2(k h) coordinates and contracts a;
r_tail     has log2(s) coordinates and later contracts j.
```

The remaining existing axes are the position point `r_M` and block point
`r_B`. The point order and descriptor MUST bind this split. Prover and verifier
MUST derive it from `(k, d_A, s)`; it is not caller-selected layout metadata.

### Partial opening

For each live block, define

```text
e_i(Y)
  = sum_(x,a,j)
      eq(r_M, x) eq(r_pack, a) f_(i,x,a,j) Y^j
  in C = E[Y]/(Y^s+1).
```

No trace map is used. Each of the `s` coefficients of `e_i` is one ordinary
element of `E`. Fix the implementation's canonical `K`-basis
`beta_0, ..., beta_(k-1)` of `E` and write

```text
e_i(Y)
  = sum_(j < s) (sum_(t < k) beta_t e_(i,t,j)) Y^j.
```

The physical base-field layout is

```text
[claim][block][extension coordinate t][carrier coefficient j].
```

It contains exactly `k s` base-field coordinates per claim/block. Backends MAY
temporarily use packed `E` values, but transcript encoding, gadget
decomposition, commitment input, range checking, and witness sizing MUST use
the canonical base-field layout above.

`Y` is a formal carrier indeterminate, not an opening-point coordinate. The
scalar opening below contracts the coefficient table with `eq(r_tail,j)`.
Ring switching later evaluates the carrier polynomial at `Y = alpha`; those are
different operations with different purposes.

### Scalar opening equation

For one polynomial with claimed opening `v`, the direct equation is

```text
sum_i eq(r_B, i)
  sum_(j < s) eq(r_tail, j)
  sum_(t < k) beta_t e_(i,t,j)
  = v.
```

After gadget decomposition at opening basis `b_open`, this becomes

```text
sum_i eq(r_B, i)
  sum_(j < s) eq(r_tail, j)
  sum_(t < k) beta_t
  sum_l b_open^l e_hat_(i,l,t,j)
  = v.
```

This replaces the current evaluation-trace formula for a direct-carrier group.
It remains one logical field-level Stage-2 row, with the existing claim-batching
coefficient applied outside the displayed equation. It has no cyclotomic
quotient. The implementation SHOULD name its prepared weights as direct
coefficient-opening weights rather than trace weights.

At a nonterminal fold, this equation is evaluated from digits authenticated
through the next recursive witness. At a transparent terminal, it is evaluated
directly from the revealed coefficient planes. A grouped root has one
schedule-owned opening geometry per group, in canonical root group order. The
precommitted profile freezes the commitment geometry; the root schedule
separately freezes how that group is opened. The verifier MUST NOT apply one
group's coefficient layout or carrier dimension to another group.

## Fold challenges and the two relation rings

### Carrier challenge

Each fold challenge is sampled as

```text
c_i(Y) in S = K[Y]/(Y^s+1)
```

using the challenge configuration audited for dimension `s`, not the one for
`d_A`. The transcript sampler MUST bind the opening mode, `s`, the challenge
configuration, group identity, live block count, and claim count before
expansion.

The same challenge is used in two rings:

```text
carrier relation:  c_i(Y)          in S;
A relation:        c_i(X^(k h))    in R.
```

No second challenge is sampled. The ambient form is a coefficient-preserving
embedding of the carrier form.

### Folded source and carrier linearity

For each position `x`, the ambient folded source is

```text
Z_x(X) = sum_i c_i(X^(k h)) F_(i,x)(X)    in R.
```

Let the direct coefficient contraction be

```text
L(F_i)(Y)
  = sum_(x,a,j)
      eq(r_M, x) eq(r_pack, a) f_(i,x,a,j) Y^j.
```

Because multiplying by `Y` advances only the `j` index, and because wrapping
`j = s` gives the same minus sign as `X^{d_A} = -1`, this map is `S`-linear:

```text
L(c_i(X^(k h)) F_i) = c_i(Y) L(F_i)    in C.
```

Therefore honest witnesses satisfy the carrier consistency equation

```text
L(Z)(Y) = sum_i c_i(Y) e_i(Y)          in C.
```

This identity is the algebraic reason the direct protocol is complete.

### Carrier consistency quotient

Use ordinary polynomial representatives of degree below `s`. Define

```text
N_eval(Y) = sum_i c_i(Y) e_i(Y) - L(G z_hat)(Y).
```

The consistency equation is equivalent to the existence of one quotient

```text
Q_eval(Y) in E[Y],  degree(Q_eval) < s,

N_eval(Y) = (Y^s + 1) Q_eval(Y).
```

This is **one quotient over `C`**, not `k` independent relation rows. If

```text
Q_eval(Y) = sum_(j < s) (sum_(t < k) beta_t q_(t,j)) Y^j,
```

then its physical witness layout is the `k` base-field coordinate planes

```text
[extension coordinate t][carrier coefficient j].
```

The quotient contributes `k s` base-field coordinates before its ordinary
gadget decomposition. Relation layout types MUST distinguish:

- one logical row selector;
- carrier modulus dimension `s`; and
- physical coordinate width `k s`.

Treating the row as a base-field ring of dimension `k s` is incorrect: it would
use the modulus `Y^{k s}+1` and the denominator `alpha^{k s}+1`.

Since `L(G z_hat)` has degree below `s`, the quotient is just the high half of
the challenge products:

```text
Q_eval = high_s(sum_i c_i e_i).
```

The prover SHOULD compute this coordinatewise over the `k` base-field planes,
sharing the sparse challenge positions across all planes.

### A rows remain ambient

The A rows do not move to the carrier. They remain

```text
[A G z_hat]_r
  = sum_i c_i(X^(k h)) [G t_hat_i]_r
  in R, for every A row r.
```

They keep `n_A` logical rows, ambient dimension `d_A`, the existing A matrix,
and the existing `t_hat` layout. Only the challenge support changes: carrier
position `j` appears at ambient position `k h j`.

The sparse challenge-times-`t` quotient can be viewed as `k h` independent
length-`s` lanes. An implementation MAY exploit those lanes, but the result
MUST match multiplication by the embedded challenge in `R` exactly.

### B and D rows

B slicing from PR #383 is unchanged. B continues to bind `t_hat` using one
physical matrix reused across its selected dyadic slices.

D continues to bind the gadget digits of the partial openings at nonterminal
folds. A transparent terminal has no D role for `e`: it reveals the `k s`
base-field coordinates per claim/block. For nonterminal direct folds, the first
implementation requires

```text
d_D divides k s.
```

This avoids a second padding convention. The number of D-role subcolumns per
partial is `selected_partial_width / d_D`. D ranks, compression source widths,
and H compression geometry MUST be recomputed from that exact width. They MUST
NOT be obtained by scaling an old `d_A` price after rank selection.

## Ring switching

### Two evaluations of each challenge

The ring-switch challenge `alpha` remains one element of `E`. The verifier MUST
derive two values from every carrier challenge:

```text
c_carrier_alpha = c_i(alpha)
  = sum_(j < s) c_(i,j) alpha^j;

c_ambient_alpha = c_i(alpha^(k h))
  = sum_(j < s) c_(i,j) alpha^(k h j).
```

The carrier consistency row uses `c_carrier_alpha`. Every A row uses
`c_ambient_alpha`.

The current single `c_alphas` cache MUST be split or typed so that these values
cannot be interchanged. Computing one and reusing it for both relations is a
protocol error except in the degenerate case `k h = 1`.

### Evaluating the carrier quotient

In native mode, the consistency check at `Y = alpha` is

```text
sum_i c_i(alpha) e_i(alpha)
  - L(G z_hat)(alpha)
  - (alpha^s + 1) Q_eval(alpha)
  = 0 in E.
```

For a coordinate-plane representation,

```text
Q_eval(alpha)
  = sum_(t < k) beta_t sum_(j < s) q_(t,j) alpha^j.
```

This fixed basis combination does not need an additional random row-batching
challenge. Before evaluation, the `beta_t` form a `K`-basis, so a nonzero set of
coordinate polynomials gives one nonzero polynomial in `E[Y]`. Random `alpha`
then tests that single polynomial. Cancellation at a particular `alpha` is
already covered by its root bound.

The prepared relation point MUST use the carrier powers
`1, alpha, ..., alpha^(s-1)` and denominator `alpha^s+1` for this row. It MUST
continue to use ambient powers and `alpha^{d_A}+1` for A rows.

### Cyclic and negacyclic products

For any degree-below-`s` product written as

```text
c(Y)e(Y) = L(Y) + Y^s H(Y),
```

the cyclic and negacyclic reductions are

```text
cyclic     = L + H,
negacyclic = L - H,
H          = (cyclic - negacyclic) / 2.
```

These identities still apply because the base characteristic is odd. They do
not, by themselves, make a new persistent cache useful. Current sparse
challenge products already compute only the high half. Native mode SHOULD
extend that high-half kernel to `k` length-`s` coordinate planes.

The existing cyclic/negacyclic setup caches for `A z` remain useful and remain
ambient. This change MUST NOT replace them with extension-field setup matrices.
D-side cache widths change with the shorter direct partial, and setup/cache
requirements MUST be derived from the selected mode.

## Soundness requirements

This section states the security obligations introduced by the new mode. It
does not replace the existing MSIS binding proof for A, B, D, F, and H.

### Transcript order

For each native group, the transcript MUST enforce this dependency order:

1. Bind the instance, schedule, mode, dimensions, coefficient layout, group
   layout, opening point, and original commitment.
2. Bind the complete D/H payload that commits to every base-field coordinate of
   every `e_i`; at a transparent terminal, bind the canonical raw `e` segment
   instead.
3. Sample the carrier challenges `c_i` at dimension `s`.
4. Bind the challenge-dependent folded witness, A/B data, carrier quotient, and
   next-witness commitment.
5. Sample `alpha`, relation-row coefficients, and later sum-check challenges.

No coordinate of `e_i` may remain unbound when `c_i` is sampled. No coordinate
of `Q_eval`, `z_hat`, or `t_hat` may remain unbound when `alpha` is sampled.
Existing labels MAY be retained only when the serialized descriptor makes the
mode and dimensions unambiguous. Otherwise new domain-separated labels are
REQUIRED.

### Challenge entropy and unit differences

Every admitted carrier challenge configuration MUST satisfy both conditions:

1. one draw has at least the configured 128-bit Fiat-Shamir min-entropy target;
2. the difference of any two distinct challenges in the family is a unit in
   `S` under Akita's audited short-invertibility bound.

The second condition MUST be checked for the **difference** family, including
its doubled coefficient and norm bounds, not merely for one sampled challenge.
The proof and parameter checker MUST use the factorization/invertibility bound
for `Y^s+1`. Entropy validation alone is insufficient.

If `delta(Y)` is a unit in `S`, then `delta(X^(k h))` is a unit in `R`: the
carrier embedding maps the inverse of `delta` to an ambient inverse. The same
`delta` is also a unit after scalar extension to `C`.

### Forking extraction

Consider two accepting transcripts with the same pre-challenge commitments and
different challenge at one claim/block position. Let

```text
delta = c_j - c'_j.
```

After subtracting the accepted A relations,

```text
A (z - z') = delta(X^(k h)) t_j       in R^(n_A).
```

In native mode, subtracting the accepted carrier consistency relations gives

```text
L(G(z - z')) = delta(Y) e_j           in C.
```

Because `delta` is a unit in both rings, these equations determine the opened
`t_j` and `e_j` from the fork. The existing B/F binding of `t_hat`, D/H binding
of `e_hat`, A binding of the folded source, range proof for all digit planes,
and quotient checks then give the same weak-opening/MSIS reduction as the
current fold.

The implementation security note MUST spell out how the standard multi-fork
argument isolates all claim/block positions. It MUST NOT claim extraction from
challenge entropy alone.

### Ring-switch polynomial check

For an honest witness, the carrier numerator is identically zero after the
quotient is included. For a false witness, it is one nonzero polynomial over
`E` of degree at most `2s-1`. Sampling `alpha` after the quotient is bound
detects it except with probability at most

```text
(2s - 1) / |E|,
```

before accounting for the existing row batching and other sum-check errors.
The coordinate basis does not multiply this error by `k`: basis independence
shows that a nonzero coordinate vector is a nonzero coefficient in `E`, and
the verifier tests the resulting single `E` polynomial.

The final theorem statement for direct mode MUST include:

- binding of the original and partial commitments;
- the carrier challenge entropy and unit-difference assumptions;
- the carrier polynomial root bound;
- the existing A/B/D/F/H MSIS assumptions;
- the existing range and sum-check soundness errors; and
- random-oracle forking loss for the complete vector of fold challenges.

## Proof-size evidence at the PR #383 baseline

### Exact current EOR formula

The current serialized EOR payload contains challenge-field partials and a
compressed degree-two sum-check. Let

```text
k  = [E:K],
P  = total number of root polynomials,
n0 = maximum root num_vars,
W1 = field-element length entering L1.
```

All current fp32/fp64 challenge fields serialize to 16 bytes. When EOR is
enabled, the exact header-free payload is

```text
L0 bytes = 16 * (k P + 2 * (n0 - log2(k)));

L1 bytes = 16 * (k + 2 * (ceil(log2(W1)) - log2(k))).
```

For one polynomial and `k` equal to 2 or 4, these simplify to

```text
L0 bytes = 32 * n0;
L1 bytes = 32 * ceil(log2(W1)).
```

These formulas are the canonical
`extension_opening_reduction_level_bytes` calculation, which is tested against
the serialized EOR payload. Removing EOR does not remove the fold grind nonce;
the numbers below count only bytes that actually disappear with the EOR proof.

### Complete current catalog census

The table expands every fp32/fp64 generated row at PR #383's head and applies
the canonical sizing function. `Current proof` is the current planner's exact
payload estimate. `L0+L1` is the gross saving if those two EOR payloads are
removed while everything else remains fixed.

| Catalog row | Current proof | L0 EOR | L1 EOR | L0+L1 | Current proof share |
|---|---:|---:|---:|---:|---:|
| fp32 dense, nv20, P=1 | 79,840 | 640 | 672 | 1,312 | 1.64% |
| fp32 dense, nv26, P=1 | 83,172 | 832 | 768 | 1,600 | 1.92% |
| fp32 one-hot, nv14, P=1 | 66,484 | 448 | 544 | 992 | 1.49% |
| fp32 one-hot, nv16, P=1 | 67,624 | 512 | 544 | 1,056 | 1.56% |
| fp32 one-hot, nv16, P=2 | 67,688 | 576 | 544 | 1,120 | 1.65% |
| fp32 one-hot, nv20, P=1 | 74,572 | 640 | 608 | 1,248 | 1.67% |
| fp32 one-hot, nv20, P=2, two groups | 77,740 | 704 | 608 | 1,312 | 1.69% |
| fp32 one-hot, nv28, P=1 | 82,388 | 896 | 736 | 1,632 | 1.98% |
| fp32 one-hot, nv30, P=1 | 83,300 | 960 | 768 | 1,728 | 2.07% |
| fp64 dense, nv14, P=1 | 79,976 | 448 | 576 | 1,024 | 1.28% |
| fp64 dense, nv20, P=1 | 86,160 | 640 | 704 | 1,344 | 1.56% |
| fp64 dense, nv26, P=1 | 88,900 | 832 | 800 | 1,632 | 1.84% |
| fp64 one-hot, nv28, P=1 | 87,232 | 896 | 736 | 1,632 | 1.87% |
| fp64 one-hot, nv30, P=1 | 87,568 | 960 | 768 | 1,728 | 1.97% |

Thus the present catalogs spend 992 to 1,728 bytes on L0/L1 EOR, or 1.28% to
2.07% of the complete proof estimate. These are gross, schedule-local savings,
not the final planner result. Direct mode also changes the next witness, ranks,
sum-check domains, and possibly the chosen schedule.

### Carrier-coordinate savings

Before digits, one direct partial and its consistency quotient each change from

```text
d_A base-field coordinates
```

to

```text
k s = d_A / h base-field coordinates.
```

The exact reduction factor is `h`. For `B` live claim/block pairs and opening
digit depth `delta_open`, the D input changes from

```text
B * (d_A / d_D) * delta_open    D-ring elements
```

to

```text
B * (k s / d_D) * delta_open    D-ring elements.
```

The carrier quotient's base-field coordinate count changes by the same factor
before quotient-digit decomposition. Compression output payloads may have fixed
sizes, so the planner MUST propagate the shorter witness through ranks,
compression chains, relation domains, successor dimensions, and proof sizing;
it MUST NOT report `h` as an automatic proof-size factor.

### Concrete fp32, `d_A = 1024`, `k = 4`

The candidates induced by the current production challenge ladder expose the
main tradeoff. Direct-mode admission still requires the new unit-difference
certificate specified above.

| `s` | `h` | `k h` ambient stride | coordinates per partial | production sparse family at `s` | challenge `l1` mass |
|---:|---:|---:|---:|---|---:|
| 64 | 4 | 16 | 256 | 31 coefficients in `±1`, 10 in `±2` | 51 |
| 128 | 2 | 8 | 512 | 31 coefficients in `±1` | 31 |
| 256 | 1 | 4 | 1,024 | 23 coefficients in `±1` | 23 |

For the middle choice, coefficient index `a + 8j` maps carrier position `j`
to ambient position `8j`. Every partial and carrier quotient uses four
length-128 base-field coordinate planes, or 512 coordinates total, instead of
1,024. The ring-switch verifier computes `c(alpha)` for the carrier row and
`c(alpha^8)` for the A rows.

The `s=64` choice gives a fourfold smaller native partial than `s=256`, but uses
a heavier carrier challenge. At a transparent terminal that partial is raw
wire. At a nonterminal it changes the digit count, D width, and successor
witness. The planner must compare both effects rather than assuming that the
smallest `s` wins.

## Planner contract

### Schedule-owned mode

Each opening group has one frozen opening mode:

```rust
enum OpeningGroupMode {
    CoefficientCarrier { carrier_dimension: usize },
}

enum TerminalOpeningPath {
    TransparentCoefficientCarrier,
    CurrentEorHachi,
}
```

Equivalent naming is acceptable, but the mode and `carrier_dimension` MUST be
typed protocol data. They MUST appear in runtime schedules, generated rows,
canonical descriptors, catalog identity, setup-prefix identity, proof-size
reports, and transcript binding.

The terminal path is level-wide. No `CommittedCoefficientCarrier` variant is
admitted: a D image plus local sum-checks does not authenticate their final
hidden-witness evaluations. A future variant MUST specify and price that
opening argument before it can enter this enum.

A scalar or recursive fold has one opening group and therefore one entry. A
grouped root stores one entry per group in `OpeningClaimsLayout::root_group_order`.
All L0 entries are direct, but their `d_A`, `s`, and derived `h` may differ.
This is not one level-wide `s` padded to the largest group.

`h` and the ambient stride are derived values. They MUST NOT be serialized as
independent choices.

### Candidate admission

A direct candidate is admitted only when all of the following hold:

- `k`, `d_A`, and `s` are powers of two;
- `k s` divides `d_A`;
- `h = d_A/(k s)` is positive;
- an audited sparse challenge configuration exists at dimension `s`;
- that configuration meets the entropy and unit-difference requirements;
- the field/ring dispatcher supports the ambient A dimension and carrier
  kernels;
- `d_D` divides `k s` in the first implementation;
- all D/H compression sources satisfy their existing byte caps;
- A, B, and D matrix widths have secure ranks at the candidate's exact norm
  bounds; and
- the resulting next-witness and relation address geometry are representable
  without unchecked padding or allocation.

The initial `s` domain SHOULD be the small, audited subset of existing
production challenge dimensions that divides `d_A/k`. The planner MUST NOT scan
every power of two or synthesize a challenge configuration during schedule
search.

For `k = 1`, EOR is invalid and contributes no bytes. The planner MAY still use
`s < d_A` to reduce partial and quotient coordinates. The full-ring baseline is
the direct candidate `s = d_A`, `h = 1`.

A transparent-terminal candidate is admitted only if every group uses the
canonical native coefficient layout and every raw coefficient count is exactly
representable. The existing EOR/Hachi terminal remains admitted under its
current checks. The planner MUST reject any candidate that prices only a D
image and local sum-checks as a no-raw-`e` terminal.

### Level policy

For presets with `k > 1`:

- absolute levels 0 and 1 enumerate native direct candidates only;
- every nonterminal level enumerates native direct candidates only; and
- the terminal compares transparent coefficient-carrier candidates with the
  current EOR/Hachi path under the same complete-suffix objective.

The absolute level is the same level used by PR #383 B-slice eligibility. A
precommitted profile freezes its commitment shape but receives its group-local
opening mode from the grouped-root schedule. A setup prefix whose witness shape
depends on a producing fold freezes that fold's complete opening-mode vector. A
later consumer MUST validate that frozen shape; it MUST NOT reinterpret stored
partial or quotient coordinates under a different carrier.

### Independent carrier and A dimensions

For native mode, the carrier dimension `s` controls challenge entropy and the
`k s` partial and quotient widths. The A dimension `d_A` controls the ambient A
rows. Once `s` is fixed, changing `d_A` changes only the packing gain and
ambient geometry; it MUST NOT enlarge the native partial or carrier quotient.

The independent `d_A`/`s` rule also applies at the transparent terminal. The
current EOR/Hachi path is priced from its own audited geometry rather than
being modeled as a coefficient-carrier encoding.

The current `AdaptiveDimension` policy searches role-specific dimensions for
two levels and then switches to a small uniform suffix. That suffix rule MUST
NOT be used for the final carrier planner. Every nonterminal and terminal state
MUST enumerate every supported `d_A` that:

- does not exceed the incoming A-dimension ceiling;
- is compatible with the selected mode and `s`;
- admits the required A challenge embedding and kernels; and
- has an audited secure rank for the exact candidate bounds.

B and D remain independent role dimensions and MAY stay below `d_A`. The
planner MUST NOT raise them merely because a larger A ring wins.

For each retained candidate, the report MUST expose both the A image width

```text
n_A d_A
```

and the exact successor-witness field-coordinate count. A larger-`d_A`
candidate MUST remain on the frontier whenever neither count is greater than
for a smaller-`d_A` candidate with the same mode, `s`, challenge family, and
input state. If both counts and every higher-priority component of the catalog
objective tie, the planner MUST choose the larger `d_A`, hence the lower-rank
representation. This rule makes dimension 256 at rank 1 win over dimension 64
at rank 4 when the rest of the exact geometry is unchanged; it does not assume
that every larger ring wins.

### Exact pricing

For every candidate, the planner MUST recompute at least:

- removed current EOR bytes;
- transparent raw-`e` bytes;
- complete current EOR/Hachi terminal bytes;
- partial coordinate count and opening gadget depth;
- D input width, secure D rank, H source, and compression geometry;
- carrier quotient coordinates and quotient digit count;
- sparse challenge `l1`, `l2`, and `linf` bounds at `s`;
- folded `z` bounds and the secure A rank;
- `t_hat`, B input width, B slicing candidates, and F compression geometry;
- logical row count, physical row dimensions, relation address length, and
  sum-check rounds;
- setup/offloading geometry;
- successor witness length and the complete suffix; and
- terminal response geometry.

The proof-size report MUST show the selected mode, `s`, partial coordinates,
carrier quotient coordinates, and net proof/setup change at every level.
Every row MUST show `h` and removed current EOR bytes. Terminal rows MUST show
the exact raw native payload and complete current EOR/Hachi payload. The report
MUST NOT list the invalid 1,744-byte subtotal as an available candidate.

### B slicing interaction

PR #383 and this change target the same first two levels but optimize different
objects:

- coefficient carriers shorten `e`, the D/H source, and the consistency
  quotient;
- B slicing shortens the one physical B matrix, while increasing the logical B
  row stack and possibly compression/relation work.

B width depends on `t_hat`, and `t_hat` depends on the carrier challenge norm,
the selected A dimension, and secure A rank. Candidate construction MUST
therefore use this order:

1. choose opening mode, `d_A`, `s`, and the carrier challenge family;
2. derive fold bounds, A rank, `z`, and `t_hat` geometry;
3. enumerate and locally prune PR #383 B slice counts;
4. derive D/H and B/F compression plans;
5. construct the next witness and score the complete suffix.

The planner MUST NOT choose a B slice count from geometry computed before the
carrier candidate. PR #383's bounded `{1,2,4,8}` slice set and its local
profitability rule remain unchanged.

### Search-space control

This change MUST NOT create an unbounded product of mode, `s`, layout, and B
slicing choices. The planner MUST:

- derive `h` from `(k,d_A,s)` for every coefficient-carrier candidate;
- use one canonical coefficient layout and embedding;
- use a fixed audited list of `s` values;
- omit EOR entirely at L0/L1;
- enumerate the full admitted A-dimension ladder at every later level instead
  of switching to a uniform suffix;
- apply security and divisibility admission before rank lookup;
- retain only undominated `(d_A, n_A d_A, successor length)` geometries for a
  fixed mode and carrier configuration;
- apply PR #383's local B-slice pruning after A/`t` geometry is known;
- retain only the existing objective's deterministic Pareto winners at each
  bounded DP state; and
- compare the pruned result with an unpruned oracle on small generated
  fixtures.

## Implementation boundaries

### `akita-types`

- Add the schedule-owned opening mode and checked carrier geometry.
- Add canonical descriptor encoding for mode and `s`.
- Represent the carrier consistency row as one logical row with carrier
  dimension `s` and extension coordinate width `k`.
- Extend witness/address layouts for `k s` partial and quotient coordinates.
- Represent the transparent terminal payload as exactly `P B k s` canonical
  base-field coordinates.
- Generalize proof-size and successor-witness sizing by selected mode.
- Keep malformed verifier inputs on typed `AkitaError` or
  `SerializationError` paths.

### `akita-challenges`

- Reuse the signed-sparse sampler at dimension `s`.
- Add or expose a parameter certificate that covers entropy and the complete
  pairwise-difference invertibility bound.
- Bind carrier dimension and mode in the draw domain.
- Do not create a second ambient challenge draw.

### `akita-prover`

- Add dense, one-hot, and recursive kernels that compute the `s` extension
  coefficients directly from the canonical coefficient split.
- Decompose the resulting `k s` base-field coordinates into D-role `e_hat`.
- Compute `Q_eval` with shared-challenge high-half accumulation over `k`
  coordinate planes.
- Keep A quotients over `d_A` with challenges embedded at stride `k h`.
- Split carrier and ambient challenge evaluations in ring-switch preparation.
- Replace the trace-specific Stage-2 term with the direct coefficient-opening
  term for direct groups.
- Keep the direct terminal transparent. Do not construct a D-only hidden tail
  without a separately specified authenticated opening argument.
- Preserve current cyclic/negacyclic A setup caches.

### `akita-verifier`

- Reconstruct the direct scalar-opening row from `r_B`, `r_tail`, the canonical
  extension basis, and opening gadget weights.
- Evaluate carrier quotient planes with denominator `alpha^s+1`.
- For native groups, evaluate the same challenge at `alpha` and
  `alpha^(k h)` for its two roles.
- Recompute every terminal direct relation from the canonical raw coefficient
  planes. Reject a hidden coefficient-tail encoding that lacks an authenticated
  final opening.
- Reject mode/dimension/layout mismatches before allocation.
- Preserve the no-panic verifier contract.

### `akita-planner`, `akita-schedules`, and `akita-config`

- Add bounded carrier candidates and the level policy above.
- Replace the two-level A search plus uniform suffix with per-level A-ladder
  enumeration and exact dominance pruning. B and D remain role-specific.
- Recompute exact ranks, setup, compression, proof bytes, and successors.
- Regenerate every affected fp32/fp64 catalog on top of PR #383.
- Add report columns for opening mode and carrier economics.
- Report the transparent-native and current Hachi/EOR terminal totals. Keep the
  invalid D-plus-local-sum-check subtotal out of candidate generation.
- Keep fp128 EOR-free, hold `s = 64` across larger A dimensions where it wins,
  and apply the same terminal search policy as the extension-field presets.

## Acceptance criteria

### Algebra and completeness

- [ ] Checked carrier geometry accepts exactly supported `(k,d_A,s)` triples and
      derives `h` and stride without independent metadata.
- [ ] Coefficient index `a + k h j` round-trips for every supported geometry.
- [ ] Dense, one-hot, and recursive direct partials match a flat MLE reference.
- [ ] `L(c(X^(k h))F) = c(Y)L(F)` holds against a naive reference for random
      small fixtures and every supported field tier.
- [ ] Direct scalar-opening weights reproduce the claimed opening, including
      partial final blocks, multiple polynomials, and multiple groups.
- [ ] `Q_eval = high_s(sum_i c_i e_i)` satisfies the full ordinary-polynomial
      divisibility identity in `E[Y]`.
- [ ] The `k` base-field quotient planes evaluate to the same `E` value as a
      packed-extension reference.
- [ ] Honest direct proofs verify at L0 and L1 for fp32 degree 4 and fp64 degree
      2.
- [ ] Transparent-terminal parsing contains exactly `P B k s` canonical
      base-field coordinates and reproduces the direct extension-field opening
      for every claim and group.

### Soundness and transcript

- [ ] Every admitted challenge family has a reviewable 128-bit entropy and
      pairwise-difference unit certificate for `S`.
- [ ] The certificate checks the difference family's exact coefficient/norm
      envelope.
- [ ] Nonterminal partial D/H payloads and transparent raw `e` are
      transcript-bound before their carrier challenge draws.
- [ ] Carrier quotient and next-witness data are bound before `alpha`.
- [ ] Mode, `s`, challenge configuration, coefficient layout, and group identity
      change the descriptor/transcript bytes.
- [ ] The verifier computes distinct `c(alpha)` and `c(alpha^(k h))` values and
      tests fail when either is substituted for the other.
- [ ] A nonzero coordinate-plane numerator is detected by the packed
      `E[Y]` ring-switch oracle.
- [ ] A verifier test rejects any terminal descriptor that claims a hidden
      D-bound coefficient tail without an authenticated opening proof.
- [ ] The direct-mode theorem adds no `1/|K|` coordinate-projection term.
- [ ] Multi-fork extraction and total soundness-error accounting are documented
      alongside the implementation.
- [ ] Malformed mode/dimension/coordinate counts return typed errors without
      panic or unbounded allocation.

### Planner and sizing

- [ ] Generated fp32/fp64 schedules contain no EOR mode at absolute L0 or L1.
- [ ] Every nonterminal uses the native direct relation; the terminal compares
      transparent native candidates with the current EOR/Hachi baseline.
- [ ] The planner has no D-only hidden coefficient-tail candidate.
- [ ] Later folds enumerate the full admitted A-dimension ladder rather than a
      uniform suffix, while retaining independent B and D dimensions.
- [ ] Reports show `n_A d_A` and exact successor length for every retained A
      geometry and select the larger `d_A` on an otherwise exact tie.
- [ ] `d_D` not dividing the selected native or hidden-digit width rejects
      before matrix/rank construction.
- [ ] Exact D/H and A/B/F ranks are recomputed from carrier geometry and norms.
- [ ] PR #383 B slicing is enumerated only after carrier-derived A/`t` geometry.
- [ ] Bounded DP output matches an unpruned oracle on small search fixtures.
- [ ] Proof reports reproduce the PR #383 baseline EOR census in this spec and
      show the new gross and net changes.
- [ ] At least one fp32 and one fp64 production row demonstrate the expected
      L0/L1 EOR removal in actual serialized proof breakdowns.
- [ ] No generated catalog silently drops a previously supported row; any row
      made unsupported by mandatory L0/L1 direct mode is listed with its failed
      admission condition for review.

### Performance and caches

- [ ] Direct partial and carrier quotient allocations contain exactly `k s`
      base-field coordinates per semantic item before digits.
- [ ] Transparent-terminal allocations and serialization contain exactly
      `P B k s` base-field coordinates.
- [ ] Carrier high-half construction does not materialize full extension-field
      convolution tables.
- [ ] Existing A cyclic/negacyclic setup caches remain shared and correct.
- [ ] D/H cache requirements use the selected carrier width and do not retain
      old `d_A`-wide buffers.
- [ ] Profile output records prover time, verifier time, peak memory, setup
      field elements, proof bytes, and per-level witness sizes for transparent
      native versus current EOR/Hachi terminal candidates.
- [ ] A packed-`E` verifier Horner loop is adopted only if it beats the canonical
      coordinate-plane loop without changing bytes or arithmetic results.

### Repository validation

- [ ] Generated schedule tables are clean after regeneration.
- [ ] Focused algebra, prover, verifier, planner, and catalog tests pass.
- [ ] All required feature-graph Clippy jobs pass.
- [ ] `./scripts/check-doc-guardrails.sh` passes.

## Non-goals

- No implementation code lands in the initial spec-only commit.
- No dependency on PR #368.
- No arbitrary carrier dimensions or coefficient-layout search.
- No second independent challenge for the A rows.
- No pure extension-field commitment or setup matrix.
- No claim that the smallest carrier is always optimal.
- No claim that flattening a native `k s` carrier and adding four public rows
  yields a sound `s`-coordinate prechallenge carrier.
- No claim that a D image plus local sum-checks forms a complete opening proof;
  short-preimage binding does not authenticate a final multilinear evaluation.
- No use of the rejected ring-valued interpolation described below: its
  opening operators are only `K`-linear and do not preserve degree over `S`.
- No requirement that a production schedule retain EOR when a fully priced
  transparent native terminal candidate wins.
- No change to PR #383's B-slice count set, dyadic partition, or 8 KiB
  compression-source limit.
- No backward-compatible decoding of schedules or proofs that predate this
  mode. Akita remains in development; affected catalogs and descriptors are
  regenerated rather than aliased.

## Documentation follow-up

The implementation PR must fold stable protocol prose into:

- `book/src/how/proving/root-fold-ring-switch.md` for the carrier relation and
  two challenge evaluations;
- `book/src/how/proving/extension-opening-reduction.md` for L0/L1 removal and
  the unresolved terminal migration boundary;
- `book/src/how/configuration.md` for planner candidates and reports;
- `book/src/foundations/rings-and-fields.md` for the carrier embedding and unit
  condition; and
- `book/src/how/security.md` for the forking and polynomial-root arguments.

Once those chapters own the durable explanation and the implementation ships,
this spec moves through `implemented` to the normal archive workflow in
[`specs/PRUNING.md`](PRUNING.md).

## References

- [B commitment slicing](commitment-slicing.md), PR #383 baseline and B-planner
  interaction.
- [Extension-field opening batching](extension-field-opening-batching.md),
  tensor EOR and the transformed-commitment soundness boundary.
- [Ring-dimension and challenge cutover](ring-dim-challenge-cutover.md), current
  production sparse families and role dimensions.
- [EOR streamed prover](eor-streamed-prover.md), current EOR prover path and
  performance context.
- [`crates/akita-types/src/layout/proof_size.rs`](../crates/akita-types/src/layout/proof_size.rs),
  canonical current EOR byte formula.
- [`crates/akita-prover/src/protocol/ring_relation/relation_quotient.rs`](../crates/akita-prover/src/protocol/ring_relation/relation_quotient.rs),
  current high-half, consistency, and A quotient construction.
- [`crates/akita-prover/src/protocol/ring_switch/relation_weights.rs`](../crates/akita-prover/src/protocol/ring_switch/relation_weights.rs),
  current structured relation weights and challenge reuse.
- [`crates/akita-verifier/src/protocol/ring_switch.rs`](../crates/akita-verifier/src/protocol/ring_switch.rs),
  current `c_alphas` preparation.
- [`crates/akita-verifier/src/protocol/evaluation_trace.rs`](../crates/akita-verifier/src/protocol/evaluation_trace.rs),
  current trace-based scalar-opening contraction.
