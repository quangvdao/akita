# Spec: Ideal Setup-Offload Pareto Planner

| Field         | Value                                      |
|---------------|--------------------------------------------|
| Author(s)     |                                            |
| Created       | 2026-07-19                                 |
| Status        | experimental design tool                   |
| PR            |                                            |
| Supersedes    |                                            |
| Superseded-by |                                            |
| Book-chapter  | book/src/roadmap/verifier-offloading.md    |

## Purpose

The ideal setup-offload planner explores the protocol we want, not the subset
accepted by today's runtime planner. Its first target is a 128-bit field
polynomial of length `2^30`, followed through the first two recursive setup
offloads. The implementation is
`crates/akita-sis-estimator/examples/ideal_setup_offload_planner.rs`.

The central design correction is that the source decomposition used to form
the inner committed witness `s` is independent of the checked decomposition
used by the opening proof. The source basis is therefore unrestricted by the
range-check basis. The planner searches it independently, including values far
above the currently shipped range-check ceiling.

This tool is deliberately separate from runtime schedule resolution. It
defines the target parameter model and exposes the frontier that should guide
the refactor. It does not claim that every returned schedule can be replayed by
the current proof format.

The planner also supports the binary one-hot root used by the CI profile. With
`--onehot-chunk-size 256`, the root source is the certified `2^num_vars`-entry
binary vector itself: its source digit count and infinity norm are both one,
and the chunk size is preserved as source identity. The one-hot root uses an
exact-support Chernoff tail specialized to the certified one-in-`K` structure,
rather than the dense `||c||_2^2 ||s||_inf^2` proxy described below.
The recursive witness and any setup prefix produced after that root are
ordinary dense sources. In
particular, the one-hot switch must not apply the full-field source-basis search
to the binary root or pretend that its `1` coefficients require 128-bit limbs.
It still searches the root `(B,P)` split, checked bases, role dimensions,
module ranks, and B/D slicing independently.

## Search Variables

At each fold the planner searches:

- A-, B-, and D-role ring dimensions `(d_A, d_B, d_D)`, subject only to
  `d_D | d_B | d_A`;
- the source decomposition basis `ell_s` independently for every committed
  input group;
- the checked outer and opening bases `(ell_t, ell_e)`, normally capped at 4
  for the first two folds but configurable without a protocol ceiling;
- every power-of-two positions-per-block value `P`, with live block count
  `B = ceil(R / P)` and no requirement that `B` be a power of two;
- the minimum SIS-secure module ranks `(n_A, n_B, n_D)` up to a configurable
  ceiling.

Increasing `d_A` remains useful until `n_A = 1`: beyond rank one it no longer
reduces the A row count. B/D dimensions remain independent storage and
security choices, but changing them only regroups the same flat witness. It
must not make `t_hat` or `e_hat` artificially smaller.

The first level has one full-field input. A recursive level has two groups: the
previous folded witness and the previous setup prefix. Each group gets its own
source basis and `(B,P)` split. The fold-level checked bases and role dimensions
are shared because they define one opening relation and one shared D matrix.

## Geometry and Cost Model

For one single-claim group with source length `N`, source digit count
`delta_s`, and `R = ceil(N / d_A)` source rings:

```text
A_width  = P * delta_s
A_fields = d_A * n_A * A_width

B_logical_width = (d_A / d_B) * B * n_A * delta_t
B_fields(full)  = d_B * n_B * B_logical_width

D_logical_width(group) = (d_A / d_D) * B * delta_e
```

The bridge factors are mandatory: `d_B * B_logical_width` and
`d_D * D_logical_width` recover the same flat coefficient counts represented
under the A ring. For a multi-group fold, A and B matrices are independent per
group and D is shared.

Commitment slicing partitions each logical B/D vector into consecutive,
value-aligned slices. B slice endpoints are multiples of
`(d_A/d_B)*delta_t`; D partitions restart at group boundaries and their
endpoints are multiples of `(d_A/d_D)*delta_e`. If `m_B` and `m_D` are the
largest physical slice widths, the reusable role matrices occupy:

```text
B_fields        = d_B * n_B * m_B
D_fields        = d_D * n_D * m_D
level_envelope  = max(max_group A_fields,
                      max_group B_fields,
                      D_fields)
prefix_fields   = next_power_of_two(level_envelope)
```

The same physical matrix is reused across all slices. SIS rank is therefore
estimated at physical width `m`, not total logical width. The current
`a-capped` slicing policy chooses the widest aligned slice whose matrix fits
under the associated A footprint; if one value run already exceeds A, it uses
that irreducible width.

The outgoing witness is `e_hat || t_hat || z_hat`. The ideal flat sizes are:

```text
e_fields = d_A * B * delta_e
t_fields = d_A * B * n_A * delta_t
z_fields = d_A * P * delta_s * delta_z
```

Thus B/D ring dimension does not alter the flat `e||t||z` size. Slicing does:
the complete stacked B image has `S_B*n_B*d_B` field coefficients and the D
image has `S_D*n_D*d_D`. Those images are decomposed for the first F/H
compression map, and the resulting digits join the next witness. The planner
currently prices the dominant first-layer suffix as

```text
B_cmp_fields = S_B * n_B * d_B * delta_cmp_B
D_cmp_fields = S_D * n_D * d_D * delta_cmp_D
```

using the role's checked basis for `delta_cmp`. This is the crucial cost that
prevents arbitrary over-slicing. Later F/H-chain digit layers and their matrix
views are not yet searched, so this is a first-layer model rather than a full
compression-chain proof-byte estimate.

The basis-sensitive witness price is:

```text
witness_bits = ell_e * (e_fields + z_fields) + ell_t * t_fields
```

This is preferable to recursive field-coordinate count: field-coordinate
count alone rewards inflating the checked basis. Exact proof bytes should
replace this proxy once the target witness encoding is implemented.

Matrix work charges every matrix application: A once, and B/D over their full
logical widths even when one physical slice matrix is reused many times. Role
storage instead uses the physical slice footprint and takes the maximum, since
the three matrices occupy one shared setup envelope. Across offload levels the
planner tracks the maximum role envelope and cumulative witness bits and work.
F/H setup views must be added before the reported envelope can be treated as
the complete protocol setup envelope.

## Canonical Folded-Response and Certified A Pricing

The planner follows the runtime fold-cap heuristic. For a flat sparse challenge:

```text
t_star^2 = 2 * B * rho_2 * ||s||_inf^2 * ln_term
```

Here `rho_2` is the production challenge's per-block squared L2 bound and the
integer `ln_term` follows the current conservative union-bound convention over
the number of folded coefficients. The pre-snap honest cap is
`min(beta_inf,t_star)`; for the searched single-claim flat profiles `t_star` is
the active side. Starting from the digit depth that covers this cap, the
planner walks `delta_z` downward while the next smaller positive clean-digit
envelope still retains at least `t_star/2`. Honest grinding uses the resulting
cap.

Security must cover every balanced digit string accepted by the verifier, not
only the honest grind cap. If the certified digit alphabet has base `b_prime`,
the response recomposes in base `b`, and there are `delta_z` planes, its exact
certified *difference* envelope is

```text
delta_cert = (b_prime - 1) * (b^delta_z - 1) / (b - 1)
```

This is the diameter seen by two accepting responses, not twice the farther
negative-side reach. In the searched matched-base schedules, `b_prime = b`, so
the tight value is exactly `b^delta_z - 1`. For the 128-bit base-field path the
ring-subfield embedding norm is one, so A-role weak binding is priced at:

```text
beta_A = 4 * ||c||_1 * delta_cert
```

B and D collision bounds remain `2^ell - 1`. `--a-collision honest-cap` is
retained only as an experimental control. Shipped candidates must use
`certified-difference`. The runtime's current `rounded_up_role_a_inf_norm`
still doubles the negative one-sided reach; it must be cut over to this same
primitive before schedule generation and verification can share one security
contract.

### Exact one-hot average-case tail

For the binary one-hot root, let `q = ceil(D/K)` be the maximum number of hot
coefficients in one ring. Fixing any output coordinate and any adversarially
chosen hot locations, the contribution of one block is the sum of the `q`
challenge coefficients at those locations. A production challenge chooses
exactly `k1` magnitude-one and `k2` magnitude-two positions uniformly without
replacement and gives every selected position an independent random sign.
Consequently the numbers `(a,b)` hit among the `q` locations have the exact
multivariate-hypergeometric probability

```text
p(a,b) = C(q,a) C(q-a,b)
         C(D-q,k1-a) C(D-q-(k1-a),k2-b)
         / (C(D,k1) C(D-k1,k2)).
```

The exact one-block moment-generating function is therefore

```text
M_q(lambda) = sum_(a,b) p(a,b) cosh(lambda)^a cosh(2 lambda)^b.
```

Independent block challenges give `M_q(lambda)^B`. For `N=P*D` folded
coefficients and the shipped grind acceptance target `1/8`, the planner uses
the best two-sided Chernoff-plus-union threshold

```text
t_onehot = ceil(inf_(lambda>0)
                  (B log M_q(lambda) + log(16 N / 7)) / lambda).
```

This is uniform over the witness's hot indices; it averages only over the
transcript-derived challenge support and signs. The optimizer solves the
unique monotone stationarity equation for `lambda`. It then caps by the exact
worst-case ring-product bound
`B * min(||c||_1, ||c||_inf*q)` and applies the same clean-digit snap-down as
the dense path, retaining at least half of the unsnapped tail threshold. A-role
security continues to price the exact certified digit *difference* envelope
after this snap, not the honest tail cap.

## Security Estimation

Every distinct `(role, ring dimension, width, coefficient bound)` cell is sent
directly to the Rust `akita-sis-estimator` with the field's exact Q32, Q64, or
Q128 modulus, infinity norm, ADPS16 quantum cost, and the estimator's default
reduction optimizer. Discovery caches the minimum secure rank. A-role pricing
also multiplies by the ring-subfield embedding norm: one for the Q128
base-field path and two for the Q64/Q32 trace-subfield paths.

The local-minimum profile is a discovery oracle, not final table
certification. Before a selected design becomes a shipped schedule, regenerate
the distinct boundary cells with the existing exhaustive width-table workflow:
certify the accepted rank and its rejected predecessor. Exhaustively searching
beta and zeta for every dominated PCS candidate is unnecessary and is the
source of the previously observed multi-minute estimates.

## Pareto State and Dominance

A partial schedule is retained when no other schedule is weakly better in all
of:

- maximum setup-envelope fields over all completed folds;
- cumulative decomposed witness bits;
- cumulative A+B+D matrix-work fields;
- outgoing folded-witness field length;
- outgoing padded setup-prefix field length.

The outgoing source bounds are part of the dynamic-programming state, so a
locally smaller proof cannot discard a schedule that leaves a cheaper future
commitment. Equal metric points are deduplicated. Output rows also report the
maximum A/B/D rank and whether every A matrix is rank one.

This is the primary frontier. Useful views are projections of the same rows:

- minimum storage envelope;
- minimum cumulative witness bits under an envelope budget;
- minimum matrix work under an envelope budget;
- the all-A-rank-one subfrontier;
- checked-basis caps 2, 3, and 4 for the early folds;
- honest-response versus certified-difference-envelope norm pricing.

## Running the Experiment

The default command searches the full two-fold `2^30` experiment:

```bash
cargo run -p akita-sis-estimator --release \
  --example ideal_setup_offload_planner
```

The important controls are:

```text
--num-vars 30
--sources LEN:BITS,...       # optional continuation-state override
--setup-offload true         # false omits the setup prefix from successor sources
--a-collision certified-difference
--levels 2
--source-basis-min 2 --source-basis-max 32
--checked-basis-min 2 --checked-basis-max 4
--fixed-checked-basis 2 --fixed-checked-basis-levels 2
--a-dims 64,128,256,512
--bd-dims 16,32,64,128
--slicing a-capped
--max-slices-per-matrix 0
--max-rank 20
--print-limit 0
```

`--slicing none` gives the corrected unsliced geometry. A nonzero slice cap is
an optional projection aid, not a substitute for the compression-suffix cost.

CSV output contains the global envelope in fields and bytes, cumulative witness
bits and bytes, matrix work, final state sizes, maximum ranks, the all-A-rank-one
flag, and a complete compact description of every level and group. It also
records which `(level, role)` attains the global envelope, the exact B/A or D/A
ratio when one of those roles dominates, and every group-local A/B matrix plus
the shared D matrix in both field elements and bytes. A print limit of zero
emits the entire frontier; a positive value emits that many rows in
storage-first order.

## Slicing Experiment on `2^30`

The corrected one-level sweep retains 485 nondominated schedules. The useful
view is minimum witness cost under progressively larger role-envelope budgets:

| envelope fields | witness bits | B slices | D slices | dominant role(s) |
|---:|---:|---:|---:|:---|
| 524,288 | 4,324,327,424 | 1,024 | 1,024 | A=B=D |
| 851,968 | 2,173,108,224 | 316 | 316 | A=B=D |
| 3,407,872 | 612,335,616 | 20 | 40 | A=B=D |
| 13,631,488 | 351,338,496 | 5 | 3 | A=B=D |
| 67,108,864 | 327,303,232 | 1 | 1 | A=B |

This is the expected shape. Moderate slicing cuts the role envelope by roughly
20x (67.1M to 3.41M fields) for less than 2x witness bits. The next 4x storage
reduction to 0.85M fields costs another 3.5x in witness bits. The extreme
storage endpoint with one value run per slice reaches 2,048 fields but carries
tens of millions of slices and is not a plausible schedule.

The 3.41M-field knee uses `(d_A,d_B,d_D)=(512,64,64)`, source basis 10,
checked bases `(4,2)`, `(B,P)=(4096,512)`, and rank one for all roles. Its B
matrix is reused over 20 slices and D over 40. The 0.85M point uses the same
role dimensions and source/checked basis 10/4/4, with `(B,P)=(16384,128)` and
316 slices for each role.

A two-level projection capped at 1,024 slices per role retains 7,127 schedules.
Representative minimum-witness points are:

| global envelope fields | cumulative witness bits | L0 B/D slices | L1 B slices (witness/prefix) | L1 D slices |
|---:|---:|:---:|:---:|---:|
| 524,288 | 4,611,686,400 | 1,024 / 1,024 | 129 / 8 | 131 |
| 1,048,576 | 2,268,614,976 | 316 / 316 | 9 / 2 | 13 |
| 3,407,872 | 658,799,200 | 20 / 20 | 2 / 1 | 2 |
| 13,631,488 | 409,051,136 | 5 / 3 | 1 / 2 | 2 |

At the 1,048,576-field point, L0 has A=B=D=851,968 fields and L1 has
A=B=1,048,576 fields with D=1,034,752. At 3,407,872 fields, L0 has all roles
equal at 3,407,872 while L1 is A-dominant at 2,981,888 fields. Thus the user's
desired A-dominant/equalized envelope does appear in the moderate-slicing
region, and the g1 setup prefix can remain at or below the same global band.

These are design estimates under honest-response and first-compression-layer
pricing. They are not certified schedule constants. In particular, complete
F/H chain matrices may enlarge the setup envelope and must be included in the
next experiment.

## Non-Offloaded Comparison

`--setup-offload false` performs a fresh schedule search with the same source
and checked bases, splits, role dimensions, slicing, SIS ranks, envelope, and
matrix-work objectives. It differs only in not appending the current setup
prefix as a full-field source for the next fold.

The unbounded two-level non-offloaded search retains 4,368 nondominated
multidimensional schedules. Their envelopes range from 2,048 to 67,108,864
fields. Projecting that set onto envelope versus cumulative witness bits gives
56 strict lower-hull breakpoints. Selected knees are:

| envelope fields | cumulative bits | matrix work | final witness fields |
|---:|---:|---:|---:|
| 524,288 | 4,605,698,048 | 2,850,816 | 139,923,456 |
| 851,968 | 2,313,519,104 | 4,915,200 | 35,102,720 |
| 1,703,936 | 1,158,851,328 | 9,826,176 | 12,120,320 |
| 3,407,872 | 641,442,352 | 17,910,272 | 6,473,232 |
| 6,815,744 | 440,975,360 | 25,001,984 | 5,607,424 |
| 11,010,048 | 402,974,320 | 28,063,296 | 5,029,072 |
| 13,631,488 | 372,195,328 | 45,056,000 | 5,214,208 |
| 31,457,280 | 350,505,008 | 80,340,928 | 5,107,024 |
| 67,108,864 | 344,393,520 | 145,934,016 | 5,055,056 |

After deriving that frontier, an offloaded comparison search with a 2,048
slice ceiling is sufficient to cover these practical knees. Minimum
cumulative witness bits under equal global-envelope budgets compare as follows:

| envelope budget | non-offloaded bits | offloaded bits | overhead | final witness fields, non/off |
|---:|---:|---:|---:|---:|
| 524,288 | 4,605,698,048 | 4,613,341,184 | 0.17% | 139,923,456 / 143,745,024 |
| 851,968 | 2,313,519,104 | 2,322,251,776 | 0.38% | 35,102,720 / 48,240,640 |
| 1,703,936 | 1,158,851,328 | 1,171,688,256 | 1.11% | 12,120,320 / 16,049,088 |
| 3,407,872 | 641,442,352 | 658,799,200 | 2.71% | 6,473,232 / 11,646,496 |
| 6,815,744 | 440,975,360 | 467,206,144 | 5.95% | 5,607,424 / 12,165,120 |
| 11,010,048 | 402,974,320 | 439,676,928 | 9.11% | 5,029,072 / 20,750,336 |
| 13,631,488 | 372,195,328 | 409,051,136 | 9.90% | 5,214,208 / 23,279,616 |
| 31,457,280 | 350,505,008 | 402,964,576 | 14.97% | 5,107,024 / 21,020,000 |
| 67,108,864 | 344,393,520 | 402,964,576 | 17.01% | 5,055,056 / 21,020,000 |

At the attractive 13.6M point both modes independently choose the same L0
schedule. Non-offloaded L1 uses `(d_A,d_B,d_D)=(128,64,64)`, checked bases
`(4,4)`, and an A-dominant 2,097,152-field envelope. Offloaded L1 uses the same
dimensions but checked bases `(4,2)` and an A=B-dominant 2,883,584-field
envelope. Consequently setup offloading increases L1 witness bits from
20,856,832 to 57,712,640 and L1 output fields from 5,214,208 to 23,279,616.
Since L0 dominates both global envelopes, this recursive overhead buys setup
offloading without increasing the schedule's maximum A/B/D storage footprint.

## Tight Certified-Difference Resweep

The certified resweep uses the exact difference envelope above and the Rust SIS
estimator for every affected boundary cell. It retains 4,196 non-offloaded and
7,463 offloaded Pareto schedules. Selected equal-envelope projections are:

| envelope budget | non-offloaded bits | offloaded bits | offload overhead |
|---:|---:|---:|---:|
| 3,407,872 | 627,458,048 | 644,251,648 | 2.68% |
| 6,815,744 | 421,721,088 | 446,330,016 | 5.84% |
| 11,010,048 | 407,377,920 | 441,999,360 | 8.50% |
| 13,631,488 | 354,299,904 | 388,351,136 | 9.61% |
| 31,457,280 | 319,856,640 | 370,667,520 | 15.89% |

The correction changes discrete ranks only at some boundaries: relative to the
naive certified-negative-reach scan, the 6.8M projection improves by 0.132%
without offloading and 0.136% with offloading, while the offloaded 13.6M point
improves by 0.118%. Other displayed knees retain the same ranks.

Dedicated first-fold checked-base `(2,2)` scans preserve the two practical
schedules:

| L0 envelope | L0 bits | best non-offloaded two-fold bits | best offloaded two-fold bits |
|---:|---:|---:|---:|
| 6,815,744 | 404,914,176 | 427,046,912 | 451,657,952 |
| 13,631,488 | 337,707,008 | 358,266,880 | 392,321,248 |

At 13.6M, L0 remains `(d_A,d_B,d_D)=(256,64,64)`, source basis 10,
checked bases `(2,2)`, split `(B,P)=(2048,2048)`, `n_A=2`, and
`delta_z=10`. Its exact A collision bound raises the reported SIS security
margin but does not cross a rank boundary. The best unconstrained continuation
uses checked bases `(4,4)` without offloading and `(4,3)` with offloading.

## Q64 and Q32 Resweep

The small-field scan keeps the root polynomial length at `2^30`, selects the
exact Q64/Q32 SIS modulus, and includes trace-subfield embedding norm two in
every A-role collision. It searches source bases through the full field width,
checked bases 2--4, `d_A` through 2,048, and `d_B,d_D` through 512. The latter
bounds reach rank one, so larger dimensions have no rank-reduction advantage.

The main source-decomposition result is:

| field | proof-efficient source log basis | rank-one source log basis | 2/2-throughout source log basis |
|:---:|---:|---:|---:|
| Q128 | 10--11 | 10 | 10 at the 13.6M knee |
| Q64 | 11 | 11 | 16 |
| Q32 | 11 | 11 | 8 |

For the unconstrained Q64/Q32 scans, log basis 11 dominates the practical
near-proof-minimum region. Q64's best 2/2 schedule instead uses four base-16-bit
source digits, while Q32's uses four base-8-bit source digits. These choices
are direct search results, not scaled Q128 parameters.

Representative non-offloaded schedules are:

| field/profile | `(d_A,d_B,d_D)` | `n_A` | `(B,P)` | source digits | checked bases | two-fold bits | envelope bytes |
|:---|:---:|---:|:---:|---:|:---:|---:|---:|
| Q64 proof minimum | `(64,32,32)` | 19 | `(1024,16384)` | 6 at log 11 | `(4,3)` | 228,602,048 | 956,301,312 |
| Q64 rank-one knee | `(1024,128,128)` | 1 | `(1024,1024)` | 6 at log 11 | `(4,4)` | 272,592,896 | 50,331,648 |
| Q64 2/2 throughout | `(128,128,32)` | 12 | `(1024,8192)` | 4 at log 16 | `(2,2)` | 233,777,152 | 402,653,184 |
| Q32 proof minimum | `(128,128,64)` | 14 | `(1024,8192)` | 3 at log 11 | `(4,3)` | 138,272,448 | 176,160,768 |
| Q32 rank-one knee | `(2048,256,256)` | 1 | `(1024,512)` | 3 at log 11 | `(4,4)` | 207,814,144 | 20,971,520 |
| Q32 2/2 throughout | `(128,128,32)` | 12 | `(1024,8192)` | 4 at log 8 | `(2,2)` | 139,899,904 | 201,326,592 |

The rank-one geometry exhibits exact inverse scaling with field width:

```text
(d_A,d_B,d_D):  (512,64,64) -> (1024,128,128) -> (2048,256,256)
field bits:       128         ->   64          ->   32
d_A * field bits: 65,536      -> 65,536        -> 65,536
```

The proof-minimum schedules do not follow this rule: they trade much larger
module rank for fewer decomposed coordinates. Relative to Q128, minimum
non-offloaded witness bits fall to 71.5% at Q64 and 43.2% at Q32, rather than
the naive 50% and 25%.

At the rank-one knees, setup offloading raises two-fold bits from 272,592,896
to 288,337,920 for Q64 (5.78%) and from 207,814,144 to 214,903,808 for Q32
(3.41%). At the more proof-aggressive approximately 100 MB storage band, the
overheads are about 10.1% and 15.9%, respectively.

## Dense Q128 Offload-Crossover Sweep

The dense Q128 crossover sweep varies the source length from `2^22` through
`2^30`. It uses the ideal mixed-dimension search above, source log bases 2--32,
checked bases 2--4, A dimensions 64--2,048, B/D dimensions 16--512, rank at
most 20, A-capped slicing with at most 2,048 slices, and the exact certified
digit-difference bound. Tensor challenges and distributed proving are disabled.

An offload round is called *profitable* only when both

```text
outgoing_witness_bits * 3 <= entering_recursive_witness_bits
next_matrix_envelope < current_matrix_envelope.
```

The witness inequality uses the exact bit length carried by the dynamic state,
not `field_len * max_checked_basis`. This distinction matters for schedules
whose outgoing groups use different checked bases.

Let `W1` be the balanced witness produced from the dense root, `W2` the witness
after opening the first setup prefix, and `W3` the witness after opening the
second setup prefix. Let `Ei` be the largest flattened A/B/D matrix at level
`i`. All following sizes are MiB, using 16 bytes per Q128 field element for
matrix storage.

The complete local first-offload projection onto `(E1,W2)` is:

| `nv` | Pareto pairs `E1 / W2` |
|---:|:---|
| 22 | `1 / 2.742`; `2.25 / 1.333` |
| 23 | `2.25 / 1.849`; `4.5 / 1.455`; `4.875 / 1.386` |
| 24 | `2.25 / 3.608`; `2.5 / 2.880`; `4.5 / 1.707` |
| 25 | `2 / 17.460`; `2.25 / 6.391`; `2.5 / 6.390`; `4.5 / 2.218`; `5 / 2.151` |
| 26 | `2.25 / 17.299`; `2.5 / 11.919`; `4.5 / 3.958`; `5 / 2.655`; `20 / 2.570` |
| 27 | `2 / 66.955`; `2.25 / 33.788`; `2.5 / 33.689`; `3 / 33.498`; `3.25 / 33.403`; `4 / 9.566`; `4.5 / 9.373`; `5 / 3.670`; `10 / 3.158`; `44 / 3.151` |
| 28 | `2.25 / 66.589`; `2.5 / 66.492`; `3 / 66.107`; `3.25 / 66.007`; `4 / 17.676`; `4.5 / 17.479`; `5 / 17.477`; `6 / 17.254`; `8 / 6.192`; `12 / 5.595`; `16 / 3.992`; `24 / 3.336` |
| 29 | `4 / 33.864`; `5 / 33.669`; `6 / 33.432`; `8 / 10.355`; `12 / 9.631`; `16 / 5.335`; `24 / 4.013` |
| 30 | `4 / 66.413`; `6 / 65.827`; `8 / 18.402`; `12 / 17.732`; `16 / 7.452`; `20 / 7.146`; `24 / 5.352`; `40 / 5.166`; `48 / 4.637` |

No qualifying first offload exists at `nv=20` or `nv=21`; `nv=22` is the
exact crossover. The `nv=22` and `nv=23` envelope reductions are modest, so
the first operationally compelling region starts around `nv=24`--`nv=26`.

For every first-offload point, the second round was freshly searched with the
first round's exact outgoing geometry and setup prefix. The maximum contraction
under `E2 < E1` is:

| `nv` | maximum `W2/W3` | result |
|---:|---:|:---|
| 22 | 1.50 | fails |
| 23 | 1.51 | fails |
| 24 | 2.54 | fails |
| 25 | 4.05 | passes |
| 26 | 5.04 | passes |
| 27 | 3.96 | passes |
| 28 | 7.80 | passes |
| 29 | 10.17 | passes |
| 30 | 6.98 | passes |

Thus `nv=25` is the exact local second-offload crossover. Maximum contraction
alone is misleading, however, because it can retain a very large intermediate
witness. Projecting all twice-profitable schedules onto initial envelope `E0`
versus total early witness `W1+W2+W3` gives:

| `nv` | `E0 / E1 / E2` | `W1 / W2 / W3` | total witness |
|---:|:---|:---|---:|
| 25 | `3 / 2.5 / 2.25` | `33.246 / 6.390 / 1.579` | 41.214 |
| 26 | `3 / 2.5 / 2` | `66.000 / 11.919 / 2.364` | 80.283 |
| 27 | `3 / 2.5 / 2` | `129.850 / 33.689 / 9.411` | 172.950 |
| 27 | `3.75 / 3.25 / 3` | `129.654 / 33.403 / 9.169` | 172.227 |
| 27 | `6 / 4.5 / 4` | `65.365 / 9.373 / 2.366` | 77.104 |
| 28 | `3.75 / 3.25 / 3` | `258.721 / 66.007 / 17.411` | 342.139 |
| 28 | `6 / 5 / 4.5` | `129.698 / 17.477 / 2.242` | 149.417 |
| 28 | `13 / 8 / 4.875` | `66.186 / 6.192 / 1.510` | 73.888 |
| 29 | `6.5 / 5 / 4` | `258.349 / 33.669 / 3.312` | 295.329 |
| 29 | `15 / 12 / 5` | `130.611 / 9.631 / 2.325` | 142.567 |
| 29 | `52 / 16 / 12` | `40.145 / 5.335 / 1.697` | 47.177 |
| 30 | `7 / 6 / 4.5` | `516.037 / 65.827 / 9.435` | 591.299 |
| 30 | `15 / 12 / 5` | `258.879 / 17.732 / 2.640` | 279.251 |
| 30 | `30 / 16 / 12` | `132.822 / 7.452 / 1.847` | 142.121 |
| 30 | `104 / 20 / 12` | `48.270 / 7.146 / 2.333` | 57.749 |

The second offload is therefore mathematically available at `nv=25`, but it
does not become a broadly attractive global tradeoff until approximately
`nv=29`. At `nv=30`, the minimum-total twice-offloaded schedule uses 57.749 MiB
of early witness instead of 44.187 MiB for the minimum-total once-offloaded
schedule, while reducing the last native envelope from 192 MiB to 12 MiB.

Offload surcharge is measured against a fresh one-source fold optimized under
the same next-envelope cap. At the minimum-total twice-offloaded points for
`nv=25` through `nv=30`, the second setup prefix adds respectively 0.574,
0.661, 0.872, 0.741, 1.003, and 1.503 MiB to the outgoing witness. The much
larger global premium at small `nv` comes from choosing a larger `W2` upstream,
not from the prefix's direct contribution at the second opening.

One concrete low-total `nv=30` schedule has envelopes `104 / 20 / 12` MiB and
A/B/D matrices `104/52/52`, `20/16/20`, and `12/12/8` MiB. Its role dimensions
are `(512,64,64)`, `(512,64,64)`, and `(128,64,64)`. The root uses source log
basis 10 and split `(B,P)=(2048,1024)`; the later balanced/setup source bases
are `(4,13)` and `(4,8)`.

## One-Hot Q128 Offload-Crossover Sweep

The corresponding one-hot sweep uses a certified one-in-256 binary root of
length `2^nv` and the exact-support MGF bound above. It scans every integer
`nv` from 24 through 32 with the same dimensions, basis ranges, slicing limit,
rank ceiling, exact bit accounting, and certified-difference collision pricing
as the dense crossover sweep. The one-hot specialization applies only at the
root; every recursive witness and setup prefix is an ordinary dense source.

One-hot exposes three distinct notions of profitability:

| requirement | first offload | both first and second offloads |
|:---|:---:|:---:|
| each witness contracts by at least 3x | `nv=24` | `nv=25` |
| each witness contracts by at least 3x and envelopes never grow | `nv=26` | `nv=29` |
| each witness contracts by at least 3x and every envelope strictly shrinks | `nv=28` | `nv=33` |

The first two entries in the witness-only column are lower bounds within the
requested scan, since it starts at `nv=24`. The strict first-offload crossover
is exact within the scanned neighborhood: `nv=24`--`nv=27` all fail
`E1 < E0`, while `nv=28` passes.

The complete strict first-offload projection onto `(E1,W2)` is:

| `nv` | Pareto pairs `E1 / W2` |
|---:|:---|
| 24--27 | none |
| 28 | `4.5 / 1.432`; `4.875 / 1.351` |
| 29 | `4.5 / 1.683`; `12 / 1.668` |
| 30 | `4.5 / 2.187`; `5 / 2.097`; `12 / 1.949` |
| 31 | `5 / 2.601`; `10 / 2.547`; `12 / 2.451`; `26.25 / 2.385`; `28.5 / 2.322` |
| 32 | `10 / 3.049`; `24 / 2.771`; `28.5 / 2.760` |

Projecting the same strict schedules onto initial envelope `E0` versus total
one-off witness `W1+W2` gives the more useful global view:

| `nv` | `E0 / E1` | `W1 / W2` | total witness |
|---:|:---|:---|---:|
| 28 | `8 / 4.875` | `4.409 / 1.416` | 5.825 |
| 29 | `8 / 4.5` | `8.444 / 1.770` | 10.213 |
| 29 | `16 / 12` | `5.008 / 1.668` | 6.676 |
| 30 | `8 / 4.5` | `16.562 / 2.187` | 18.750 |
| 30 | `16 / 12` | `8.800 / 2.036` | 10.836 |
| 31 | `8 / 5` | `32.625 / 3.198` | 35.823 |
| 31 | `16 / 12` | `17.031 / 2.451` | 19.482 |
| 31 | `32 / 28.5` | `9.541 / 2.386` | 11.927 |
| 32 | `32 / 28.5` | `18.016 / 2.760` | 20.775 |
| 32 | `64 / 44` | `11.036 / 3.099` | 14.135 |

Against a fresh one-source continuation optimized under the same `E1` cap,
the setup prefix adds 0.732, 1.003, 1.002, 1.470, and 2.031 MiB at the
minimum-total strict points for `nv=28` through `nv=32`. These are 16.6%,
20.0%, 11.4%, 15.4%, and 18.4% of the corresponding entering witnesses.
Storage-aggressive points have smaller relative surcharges; for example the
`nv=31`, `E0/E1=8/5` point pays only 0.767 MiB, or 2.3% of `W1`.

The full three-level dynamic program is necessary for the second-offload
study. Continuing only the local `(E1,W2)` projection misses schedules whose
first-round pair is dominated but whose exact outgoing digit geometry is
cheaper at level two. The authoritative twice-offloaded frontier first has any
3x/3x witness schedule at `nv=25`; its early shape is:

| `nv` | representative `E0 / E1 / E2` | `W1 / W2 / W3` |
|---:|:---|:---|
| 25 | `0.75 / 1 / 2` | `8.350 / 2.541 / 0.771` |
| 26 | `1 / 1 / 2` | `8.312 / 2.525 / 0.771` |
| 27 | `2 / 1 / 2` | `8.250 / 2.701 / 0.794` |
| 28 | `2 / 2 / 4` | `16.375 / 2.675 / 0.784` |

These are witness-profitable but not setup-envelope-profitable. Requiring both
envelopes to be non-growing delays the crossover to `nv=29`. The Pareto
frontier in initial envelope versus total witness is then:

| `nv` | `E0 / E1 / E2` | `W1 / W2 / W3` | total witness |
|---:|:---|:---|---:|
| 29 | `2 / 2 / 2` | `33.125 / 6.206 / 1.389` | 40.720 |
| 30 | `2 / 2 / 2` | `66.107 / 11.792 / 2.095` | 79.994 |
| 30 | `4 / 4 / 4` | `32.500 / 5.027 / 1.508` | 39.035 |
| 31 | `2 / 2 / 2` | `130.156 / 34.138 / 6.380` | 170.674 |
| 31 | `4 / 4 / 4` | `65.484 / 6.437 / 1.328` | 73.249 |
| 32 | `2 / 2 / 2` | `260.188 / 66.763 / 17.552` | 344.502 |
| 32 | `4 / 4 / 4` | `129.312 / 17.224 / 2.019` | 148.555 |
| 32 | `8 / 8 / 4.875` | `64.750 / 5.710 / 1.463` | 71.923 |
| 32 | `16 / 4 / 4` | `33.062 / 3.846 / 1.176` | 38.084 |
| 32 | `32 / 4 / 4` | `18.016 / 3.908 / 1.184` | 23.108 |

No schedule through `nv=32` makes *both* envelope inequalities strict. At
`nv=32`, however, the last row makes the first envelope shrink eightfold and
holds the second flat while contracting witnesses by 4.61x and 3.30x. Thus the
next integer is the natural strict-crossover candidate, while `nv=32` is the
practical equality boundary.

Relative to the minimum-total once-offloaded points, the minimum-total
non-growing twice-offloaded schedules cost 510%, 260%, 514%, and 63% more
early witness at `nv=29` through `nv=32`, while reducing the last native
envelope from respectively 12, 12, 28.5, and 44 MiB to 2, 4, 4, and 4 MiB.
Only the `nv=32` point is an immediately compelling global tradeoff; the
earlier points buy a small verifier envelope by deliberately retaining a much
larger intermediate witness.

The concrete `nv=32`, `32/4/4` schedule uses exact A/B/D sizes
`32/32/32`, `4/4/4`, and `4/4/4` MiB. Its role dimensions are
`(256,64,64)`, `(256,64,64)`, and `(256,32,32)`. The root split is
`(B,P)=(2048,8192)` with checked bases 4/4; the next two folds use balanced
source log basis 4 and checked bases 4/4 then 2/2.

### Strict second-offload crossover at `nv=33`

The full-state `nv=33` continuation retains 362 schedules after both exact
3x contraction gates. Eighteen also satisfy `E1 < E0` and `E2 < E1`, proving
that `nv=33` is the exact strict second-offload crossover immediately after
the failing `nv=32` scan. Projecting these schedules onto initial envelope
`E0` versus total early witness gives four knees:

| `E0 / E1 / E2` | `W1 / W2 / W3` | contractions | total witness |
|:---|:---|:---|---:|
| `16 / 8 / 4.5` | `65.125 / 6.130 / 1.534` | `10.62x / 4.00x` | 72.789 |
| `32 / 5.5 / 4.5` | `34.031 / 4.462 / 1.423` | `7.63x / 3.13x` | 39.917 |
| `64 / 8 / 4.5` | `20.008 / 4.891 / 1.455` | `4.09x / 3.36x` | 26.354 |
| `128 / 11 / 10` | `16.004 / 5.258 / 1.738` | `3.04x / 3.02x` | 23.000 |

The direct second-prefix surcharges under matched `E2` caps are 0.764, 0.810,
0.801, and 1.050 MiB, respectively, or 12.5%, 18.2%, 16.4%, and 20.0% of the
entering `W2`. The minimum-total once-offloaded `nv=33` schedule uses 20.069 MiB
of early witness with envelopes `128/88`. The minimum-total twice-offloaded
schedule raises early witness to 23.000 MiB, a 14.6% premium, while reducing
the last native envelope from 88 MiB to 10 MiB. This is the first strict
two-offload point and is already a compelling global tradeoff.

## One-Hot-16 Q128 Offload-Crossover Sweep

The one-in-16 sweep changes only the certified root chunk size from 256 to 16
and scans every integer `nv` from 24 through 35. All field, basis, dimension,
rank, slicing, security, exact-bit, and 3x contraction settings remain
identical to the one-in-256 sweep. The denser hot support increases the root
response bound and leaves more continuation geometries alive, but the same
discrete setup-envelope shelves remain visible.

The crossover summary is:

| requirement | first offload | both first and second offloads |
|:---|:---:|:---:|
| each witness contracts by at least 3x | `nv=24` | `nv=25` |
| each witness contracts by at least 3x and envelopes never grow | `nv=26` | `nv=29` |
| each witness contracts by at least 3x and every envelope strictly shrinks | `nv=28` | `nv=32` |

The first two witness-only entries are lower bounds within the requested scan,
which starts at `nv=24`. The strict first crossover is the same as one-in-256,
but the strict second crossover occurs one power earlier: `nv=32` instead of
`nv=33`.

The strict first-offload envelope-versus-total-witness frontier has the
following endpoints. When the minimum-envelope and minimum-total points differ,
both are shown:

| `nv` | `E0 / E1` | `W1 / W2` | total witness |
|---:|:---|:---|---:|
| 24--27 | none | | |
| 28 | `8 / 4.875` | `4.516 / 1.351` | 5.866 |
| 29 | `8 / 4.5` | `8.531 / 1.683` | 10.214 |
| 29 | `16 / 12` | `5.008 / 1.668` | 6.676 |
| 30 | `8 / 4.5` | `16.700 / 2.375` | 19.075 |
| 30 | `16 / 12` | `9.016 / 1.949` | 10.965 |
| 31 | `8 / 5` | `32.837 / 3.567` | 36.404 |
| 31 | `32 / 28.5` | `10.008 / 2.322` | 12.330 |
| 32 | `32 / 24` | `18.331 / 2.889` | 21.220 |
| 32 | `64 / 48` | `12.004 / 3.009` | 15.013 |
| 33 | `32 / 24` | `34.412 / 3.477` | 37.889 |
| 33 | `128 / 88` | `16.004 / 4.065` | 20.069 |
| 34 | `12 / 8` | `513.709 / 34.387` | 548.096 |
| 34 | `256 / 96` | `21.034 / 5.477` | 26.511 |
| 35 | `12 / 8` | `1027.042 / 66.730` | 1093.772 |
| 35 | `512 / 192` | `30.096 / 7.697` | 37.793 |

The full-state three-level search retains no strict twice-offloaded schedule
through `nv=31`. At `nv=32` exactly one survives:

```text
E0/E1/E2 = 16 / 8 / 5 MiB
W1/W2/W3 = 33.325 / 4.497 / 1.428 MiB
contraction = 7.41x / 3.15x
```

This establishes the exact strict crossover. Its A/B/D matrices occupy
`16/16/15.96`, `8/8/8`, and `5/5/4.5` MiB, with role dimensions
`(256,64,64)`, `(512,64,64)`, and `(256,64,64)`. The root uses checked bases
4/3 and split `(B,P)=(4096,4096)`; both recursive folds use balanced source
log basis 4 and checked bases 4/4.

Above the crossover, the strict initial-envelope-versus-total-witness Pareto
frontier is:

| `nv` | `E0 / E1 / E2` | `W1 / W2 / W3` | total witness |
|---:|:---|:---|---:|
| 32 | `16 / 8 / 5` | `33.325 / 4.497 / 1.428` | 39.249 |
| 33 | `12 / 8 / 4.5` | `257.042 / 18.427 / 2.305` | 277.774 |
| 33 | `16 / 8 / 4.5` | `65.438 / 8.247 / 1.660` | 75.344 |
| 33 | `32 / 5.5 / 4.5` | `34.412 / 4.871 / 1.447` | 40.731 |
| 33 | `64 / 8 / 4.5` | `20.572 / 5.189 / 1.471` | 27.232 |
| 33 | `128 / 11 / 10` | `16.004 / 5.258 / 1.738` | 23.000 |
| 34 | `12 / 8 / 5` | `513.709 / 34.521 / 3.323` | 551.553 |
| 34 | `24 / 16 / 12` | `257.084 / 11.851 / 2.122` | 271.057 |
| 34 | `32 / 16 / 12` | `66.594 / 6.466 / 1.784` | 74.844 |
| 34 | `64 / 16 / 12` | `36.644 / 5.169 / 1.680` | 43.492 |
| 34 | `128 / 16 / 12` | `25.069 / 6.073 / 1.763` | 32.905 |
| 34 | `256 / 32 / 28.5` | `21.034 / 6.822 / 2.170` | 30.026 |
| 35 | `32 / 16 / 12` | `259.250 / 11.913 / 2.139` | 273.302 |
| 35 | `48 / 32 / 28.5` | `257.667 / 9.603 / 2.299` | 269.568 |
| 35 | `64 / 16 / 12` | `69.047 / 7.401 / 1.847` | 78.295 |
| 35 | `256 / 32 / 28.5` | `34.066 / 7.232 / 2.188` | 43.486 |
| 35 | `512 / 44 / 24` | `30.034 / 8.821 / 2.935` | 41.790 |

The `nv=32` crossover is mathematically strict but not yet the minimum-total
global choice: it uses 39.249 MiB of early witness versus 15.013 MiB for the
minimum-total once-offloaded schedule, a 161% premium, while reducing the last
native envelope from 48 MiB to 5 MiB. At `nv=33`, twice offloading costs only
14.6% more early witness and reduces 88 MiB to 10 MiB. The corresponding
premiums at `nv=34` and `nv=35` are 13.3% and 10.6%, reducing the final native
envelopes from 96 to 28.5 MiB and from 192 to 24 MiB. Thus strict feasibility
begins at `nv=32`, while broad global attractiveness begins at `nv=33`.

For the minimum-total strict points at `nv=32` through `nv=35`, the direct
second-prefix surcharges under matched `E2` caps are 0.814, 1.050, 1.471, and
2.142 MiB, respectively. These are 18.1%, 20.0%, 21.6%, and 24.3% of the
entering `W2`; most of the remaining global tradeoff comes from choosing the
upstream geometry needed to keep the next setup envelope small.

### Forcing checked `2/2` in the first two folds

This constrained scan fixes both checked bases, `lt=lo=2`, at `L0` and `L1`.
The arbitrary source decomposition basis remains fully searched from 2 through
32. In a three-level continuation, `L2` is deliberately left free over checked
bases 2 through 4: it is the fold after the two folds covered by the policy.

The crossover changes are:

| requirement | unrestricted | first two folds checked `2/2` |
|:---|:---:|:---:|
| first offload: 3x witness contraction | `nv=24` | `nv=24` |
| first offload: 3x and non-growing envelope | `nv=26` | `nv=27` |
| first offload: 3x and strict envelope shrink | `nv=28` | `nv=30` |
| two offloads: both witnesses contract 3x | `nv=25` | `nv=25` |
| two offloads: 3x and both envelopes non-growing | `nv=29` | `nv=29` |
| two offloads: 3x and both envelopes strictly shrink | `nv=32` | `nv=32` |

Thus the constraint materially delays the first strict offload, by two powers
of two, but does not delay the strict two-offload crossover provided the third
fold can choose its own checked basis. The minimum-total strict points are:

| `nv` | once: `E0/E1` | once: `W1/W2` | once total | twice: `E0/E1/E2` | twice: `W1/W2/W3` | twice total |
|---:|:---|:---|---:|:---|:---|---:|
| 30 | `16/9` | `9.281/2.364` | 11.645 | none | | |
| 31 | `32/20` | `10.516/2.613` | 13.129 | none | | |
| 32 | `64/20` | `13.008/3.357` | 16.365 | `16/5.5/5` | `33.375/5.760/1.876` | 41.011 |
| 33 | `128/90` | `17.006/4.525` | 21.531 | `128/16/12` | `17.006/5.668/1.851` | 24.524 |
| 34 | `128/48` | `29.012/5.065` | 34.077 | `128/16/12` | `29.012/6.044/1.886` | 36.942 |
| 35 | `768/224` | `36.002/10.853` | 46.855 | `384/52/20` | `42.006/9.861/3.202` | 55.069 |

All sizes are MiB. Relative to the unrestricted minimum-total strict frontier,
the once-offloaded totals rise by 6.2%, 6.5%, 9.0%, 7.3%, 28.5%, and 24.0%
at `nv=30` through `nv=35`. The twice-offloaded totals rise by 4.5%, 6.6%,
23.0%, and 31.8% at `nv=32` through `nv=35`.

The `nv=32` strict two-offload boundary is the important geometry. Its exact
A/B/D sizes are `16/8/8`, `5.5/5.5/5.5`, and `5/5/5` MiB. `L0` and `L1`
use role dimensions `256/32/32` with checked `2/2`; their main splits are
respectively `(B,P)=(4096,4096)` and `(534,1024)`. `L2` selects dimensions
`256/64/64` and checked `4/4`. At `nv=33` through `nv=35`, the minimum-total
continuations instead choose `L2` checked `4/2`. This later freedom is what
preserves the `nv=32` crossover.

If checked `2/2` is also forced at `L2`, no strict twice-offloaded schedule
survives at `nv=32`; the crossover moves to `nv=33`. Above that boundary the
minimum totals barely change (24.525, 36.943, and 55.072 MiB), showing that
the main benefit of a freer third fold is recovering the marginal `nv=32`
geometry rather than broadly shrinking the later frontier.

For the minimum-total strict twice-offloaded points at `nv=32` through `nv=35`,
the direct second-prefix surcharges under matched `E2` caps are 0.831, 1.002,
1.002, and 2.112 MiB. These are 14.4%, 17.7%, 16.6%, and 21.4% of entering
`W2`. The larger total-frontier premiums therefore mostly come from the
upstream splits and ring dimensions needed to sustain checked `2/2`, not from
the full-field prefix by itself.

## Refactor Target

The runtime cutover should proceed in this order:

1. Split source decomposition from checked opening decomposition in durable
   parameter types and transcript identity.
2. Make source digit depth group-local and allow arbitrary source bases.
3. Represent B/D role-ring bridges explicitly and test that changing `d_B` or
   `d_D` preserves the corresponding flat witness size.
4. Add value-aligned B/D partitions, physical widths, and slice counts to the
   authenticated schedule; derive ranks from physical widths.
5. Carry and range-check the complete ordered F/H compression suffix, and
   include every F/H matrix in setup-envelope and work accounting.
6. Make `(d_A,d_B,d_D)` and A/B/D/F/H ranks planner outputs rather than preset
   inputs; stop raising `d_A` after rank one unless another frontier metric
   improves.
7. Replace the single local winner with Pareto state propagation over envelope,
   encoded witness size, matrix work, and outgoing state.
8. Add exact target-format proof-byte accounting while retaining the
   basis-sensitive bit metric as an auditable lower-level component.
9. Exhaustively certify only the SIS boundary cells used by selected shipped
   schedules, then emit those cells and schedules through the canonical table
   generator.
10. Replace runtime A-role pricing with the exact certified response-difference
   envelope used by this frontier, then regenerate the shipped SIS cells before
   enabling verification.

## Relationship to the Current Planner Spec

`specs/setup-offloading-planner.md` describes the constrained first rollout:
uniform D64, a fixed recursion window, local suffix minimization, and current
generated catalogs. This document does not amend those implementation rules.
It defines the desired successor design once those rollout constraints are
removed.
