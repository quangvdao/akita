# Terminal compression boundary for coefficient carriers

This note records the smallest terminal design currently justified by Akita's
protocol and binding arguments. The normative protocol and planner requirements
live in
[`specs/subring-coefficient-carrier.md`](specs/subring-coefficient-carrier.md).

## Established direct terminal

For fp32 with extension degree four and carrier dimension 64, one partial
opening has four ordinary coefficient planes:

```text
e_(i,0), e_(i,1), e_(i,2), e_(i,3)
    in S = K[Y] / (Y^64 + 1).
```

Together they represent

```text
e_i(Y) = sum_(t < 4) beta_t e_(i,t)(Y)
    in E[Y] / (Y^64 + 1).
```

This representation supports both required checks:

1. The carrier-consistency relation uses the same subring challenge in every
   coordinate plane.
2. The scalar opening uses the canonical `K`-coordinates of `E` and the
   ordinary MLE equality weights.

At a nonterminal fold, the next recursively committed witness authenticates
the gadget digits of all four planes. At a direct terminal, the prover reveals
all four planes, and the verifier computes both checks directly. No trace map
or extension-field projection is used.

## Why a D image is not a terminal opening

Suppose the prover hides the digit table `ehat`, sends

```text
v_D = D ehat,
```

and uses sum-check for its range and relation claims. The final sum-check
equation still needs `ehat(r_sc)` at a verifier challenge `r_sc`. The verifier
cannot derive that value from `v_D`.

Including `D ehat = v_D` in the same sum-check does not fix the problem. That
check also ends at the unauthenticated value `ehat(r_sc)`. MSIS binding says
that two known short preimages of the same D image yield a forbidden short
kernel vector. A single claimed MLE evaluation is not such a preimage.

A sound hidden terminal therefore needs a separate opening argument for the
committed digit table. Akita does not currently have such a coefficient-native
terminal argument.

## Why one projected plane is not enough

An `S`-linear projection such as

```text
p_i = e_(i,0) + eta e_(i,1) + eta^2 e_(i,2) + eta^3 e_(i,3)
```

does preserve carrier consistency:

```text
sum_i c_i p_i
  = L_0(z) + eta L_1(z) + eta^2 L_2(z) + eta^3 L_3(z)
  in S.
```

It does not make the extension-field MLE opening computable from `p_i`.
Negacyclic multiplication by `eta` acts on the carrier coefficient index. The
extension-field opening weights act on the MLE indices and mix the four field
coordinates. In general,

```text
Open(eta e) != eta Open(e).
```

Therefore the projected plane checks only a projection of the fold relation.
It is not a terminal opening protocol, and the planner must not price it as
one.

## Current terminal choices

For one claim, six live blocks, extension degree four, carrier dimension 64,
and four-byte base-field elements, a transparent direct terminal sends

```text
1 * 6 * 4 * 64 * 4 = 6,144 bytes
```

of raw `e`. The representative current EOR/Hachi path sends 3,072 raw `e`
bytes and a 544-byte EOR proof, or 3,616 opening-side bytes.

| Established terminal path | Opening-side bytes |
|---|---:|
| Transparent coefficient carrier at `s = 64` | 6,144 |
| Current representative EOR/Hachi path | 3,616 |

The planner may compare these complete paths. It must not compare either path
with an incomplete D-only or one-plane construction.

## Implementation boundary

The coefficient-carrier work can proceed at the root and at nonterminal folds:
it removes EOR, keeps every coordinate plane explicit, and authenticates those
planes through the ordinary recursive witness. At the terminal, the supported
coefficient-native path remains transparent. The current EOR/Hachi path remains
the compressed baseline.

Any future hidden or projected terminal must specify the prior commitment, the
authenticated final-opening equation, the extractor, the complete soundness
error, proof bytes, prover time, verifier time, setup cost, and planner
admission rules before it becomes a protocol candidate.
