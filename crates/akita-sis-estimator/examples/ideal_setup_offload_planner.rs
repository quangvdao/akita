//! Ideal-state Pareto planner for recursive setup-prefix offloading.
//!
//! This intentionally does not inherit the runtime planner's current D64,
//! uniform-role-dimension, shared-basis, or generated-table restrictions.  It
//! is an offline design tool: every listed role dimension, checked basis,
//! source basis, and block split is searched, while SIS ranks are obtained
//! directly from the Rust estimator and cached.

use akita_sis_estimator::{
    akita_q128, akita_q32, akita_q64, estimate, Adps16Mode, Bound, CostValue, EstimateConfig,
    ReductionCostModel, SisNorm, SisParameters,
};
use akita_types::SisMatrixRole;
use num_bigint::BigUint;
use num_traits::ToPrimitive;
use std::collections::HashMap;
use std::env;

const TARGET_BITS: f64 = 128.0;

#[derive(Clone, Debug)]
struct Config {
    num_vars: u32,
    initial_sources: Option<Vec<Source>>,
    setup_offload: bool,
    setup_offload_levels: Option<u32>,
    min_offload_contraction: u32,
    a_collision: ACollisionMode,
    offload_levels: usize,
    field_bits: u32,
    onehot_chunk_size: u32,
    witness_chunks: u32,
    witness_chunk_levels: Option<u32>,
    tensor_levels: u32,
    tensor_onehot_bound: TensorOnehotBound,
    source_basis_min: u32,
    source_basis_max: u32,
    checked_basis_min: u32,
    checked_basis_max: u32,
    fixed_checked_basis: Option<u32>,
    fixed_checked_basis_levels: u32,
    a_dims: Vec<u32>,
    bd_dims: Vec<u32>,
    slicing: SlicingMode,
    max_slices_per_matrix: u128,
    max_rank: u32,
    print_limit: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            num_vars: 30,
            initial_sources: None,
            setup_offload: true,
            setup_offload_levels: None,
            min_offload_contraction: 1,
            a_collision: ACollisionMode::CertifiedDifference,
            offload_levels: 2,
            field_bits: 128,
            onehot_chunk_size: 0,
            witness_chunks: 1,
            witness_chunk_levels: None,
            tensor_levels: 0,
            tensor_onehot_bound: TensorOnehotBound::SparseProxy,
            source_basis_min: 2,
            source_basis_max: 32,
            checked_basis_min: 2,
            checked_basis_max: 4,
            fixed_checked_basis: None,
            fixed_checked_basis_levels: 0,
            a_dims: vec![64, 128, 256, 512],
            bd_dims: vec![16, 32, 64, 128],
            slicing: SlicingMode::ACapped,
            max_slices_per_matrix: 0,
            max_rank: 20,
            print_limit: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SlicingMode {
    None,
    ACapped,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ACollisionMode {
    HonestCap,
    CertifiedDifference,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TensorOnehotBound {
    Generic,
    SparseProxy,
}

impl TensorOnehotBound {
    const fn name(self) -> &'static str {
        match self {
            Self::Generic => "generic",
            Self::SparseProxy => "onehot-sparse-proxy",
        }
    }
}

impl ACollisionMode {
    const fn name(self) -> &'static str {
        match self {
            Self::HonestCap => "honest-cap",
            Self::CertifiedDifference => "certified-difference",
        }
    }
}

impl SlicingMode {
    const fn name(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::ACapped => "a-capped",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Source {
    field_len: u128,
    value_bits: u32,
    bit_len: u128,
    onehot_chunk_size: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct RoleDims {
    a: u32,
    b: u32,
    d: u32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct RankKey {
    role: SisMatrixRole,
    d: u32,
    width: u64,
    bound: u128,
}

#[derive(Clone, Copy, Debug)]
struct RankEstimate {
    rank: u32,
    rop_log2: f64,
}

#[derive(Clone, Copy, Debug)]
struct SlicedMatrixPlan {
    physical_width: u128,
    values_per_slice: u128,
    rank: u32,
    rop_log2: f64,
    fields: u128,
}

struct RankOracle {
    config: EstimateConfig,
    modulus: BigUint,
    max_rank: u32,
    cache: HashMap<RankKey, Option<RankEstimate>>,
    estimator_calls: u64,
}

impl RankOracle {
    fn new(max_rank: u32, field_bits: u32) -> Self {
        let modulus = match field_bits {
            32 => akita_q32(),
            64 => akita_q64(),
            128 => akita_q128(),
            _ => panic!("unsupported Akita SIS field width {field_bits}"),
        };
        Self {
            config: EstimateConfig {
                red_cost_model: ReductionCostModel::Adps16 {
                    mode: Adps16Mode::Quantum,
                },
                ..EstimateConfig::default()
            },
            modulus,
            max_rank,
            cache: HashMap::new(),
            estimator_calls: 0,
        }
    }

    fn min_rank(&mut self, key: RankKey) -> Option<RankEstimate> {
        if let Some(cached) = self.cache.get(&key) {
            return *cached;
        }
        let mut result = None;
        for rank in 1..=self.max_rank {
            let Some(n) = rank.checked_mul(key.d) else {
                break;
            };
            let Some(m) = key.width.checked_mul(u64::from(key.d)) else {
                break;
            };
            let Ok(params) = SisParameters::try_new(
                n,
                self.modulus.clone(),
                Some(m),
                Bound::Integer(BigUint::from(key.bound)),
                SisNorm::Infinity,
            ) else {
                break;
            };
            self.estimator_calls += 1;
            let Ok(cost) = estimate(&params, &self.config) else {
                break;
            };
            if secure(cost.rop) {
                result = Some(RankEstimate {
                    rank,
                    rop_log2: cost.rop.log2().unwrap_or(f64::INFINITY),
                });
                break;
            }
        }
        self.cache.insert(key, result);
        result
    }
}

#[allow(clippy::too_many_arguments)]
fn fit_sliced_matrix(
    role: SisMatrixRole,
    d: u32,
    logical_width: u128,
    value_run: u128,
    bound: u128,
    target_fields: u128,
    slicing: SlicingMode,
    oracle: &mut RankOracle,
) -> Option<SlicedMatrixPlan> {
    if logical_width == 0 || value_run == 0 || !logical_width.is_multiple_of(value_run) {
        return None;
    }

    let total_values = logical_width / value_run;
    let estimate = |physical_width: u128, oracle: &mut RankOracle| -> Option<SlicedMatrixPlan> {
        let width = u64::try_from(physical_width).ok()?;
        let rank = oracle.min_rank(RankKey {
            role,
            d,
            width,
            bound,
        })?;
        Some(SlicedMatrixPlan {
            physical_width,
            values_per_slice: physical_width / value_run,
            rank: rank.rank,
            rop_log2: rank.rop_log2,
            fields: physical_width
                .saturating_mul(u128::from(d))
                .saturating_mul(u128::from(rank.rank)),
        })
    };

    if slicing == SlicingMode::None {
        return estimate(logical_width, oracle);
    }

    let mut rank_guess = 1u32;
    let mut best = None;
    for _ in 0..=oracle.max_rank {
        let denominator = u128::from(d)
            .saturating_mul(u128::from(rank_guess))
            .saturating_mul(value_run);
        let values = target_fields
            .checked_div(denominator)
            .unwrap_or(0)
            .min(total_values);
        if values == 0 {
            break;
        }
        let candidate = estimate(values.saturating_mul(value_run), oracle)?;
        if candidate.fields <= target_fields
            && best.is_none_or(|current: SlicedMatrixPlan| {
                candidate.physical_width > current.physical_width
            })
        {
            best = Some(candidate);
        }
        if candidate.rank == rank_guess {
            break;
        }
        rank_guess = candidate.rank;
    }

    best.or_else(|| estimate(value_run, oracle))
}

#[derive(Clone, Debug)]
struct GroupPlan {
    source: Source,
    log_source: u32,
    log_outer: u32,
    log_open: u32,
    positions: u128,
    blocks: u128,
    digits_source: u128,
    digits_outer: u128,
    digits_open: u128,
    digits_fold: u128,
    unsnapped_fold_cap: u128,
    honest_fold_cap: u128,
    tensor_low_len: u128,
    n_a: u32,
    n_b: u32,
    a_rop_log2: f64,
    b_rop_log2: f64,
    a_fields: u128,
    b_fields: u128,
    b_logical_width: u128,
    b_physical_width: u128,
    b_slices: u128,
    b_compression_fields: u128,
    d_logical_width: u128,
    next_fields: u128,
    witness_bits: u128,
    matrix_work_fields: u128,
}

#[derive(Clone, Debug)]
struct PartialLevel {
    groups: Vec<GroupPlan>,
    max_a_fields: u128,
    max_b_fields: u128,
    d_logical_width: u128,
    max_d_logical_width: u128,
    next_fields: u128,
    witness_bits: u128,
    matrix_work_fields: u128,
}

#[derive(Clone, Debug)]
struct LevelPlan {
    dims: RoleDims,
    log_outer: u32,
    log_open: u32,
    groups: Vec<GroupPlan>,
    n_d: u32,
    d_rop_log2: f64,
    a_fields: u128,
    b_fields: u128,
    d_fields: u128,
    d_logical_width: u128,
    d_physical_width: u128,
    d_slices: u128,
    compression_suffix_fields: u128,
    envelope_fields: u128,
    prefix_fields: u128,
    next_witness: Source,
    witness_bits: u128,
    matrix_work_fields: u128,
}

#[derive(Clone, Debug)]
struct State {
    sources: Vec<Source>,
    levels: Vec<LevelPlan>,
    global_envelope_fields: u128,
    cumulative_witness_bits: u128,
    cumulative_matrix_work_fields: u128,
}

#[derive(Clone, Copy, Debug)]
struct ChallengeProfile {
    ring_d: u32,
    count_pm1: u32,
    count_pm2: u32,
    l1: u128,
    l2_sq: u128,
}

fn main() {
    let config = parse_args();
    validate_config(&config);
    let mut oracle = RankOracle::new(config.max_rank, config.field_bits);
    let root_sources = config.initial_sources.clone().unwrap_or_else(|| {
        let root_len = 1u128
            .checked_shl(config.num_vars)
            .expect("num-vars must fit u128");
        let value_bits = if config.onehot_chunk_size == 0 {
            config.field_bits
        } else {
            1
        };
        vec![Source {
            field_len: root_len,
            value_bits,
            bit_len: root_len.saturating_mul(u128::from(value_bits)),
            onehot_chunk_size: config.onehot_chunk_size,
        }]
    });
    let mut states = vec![State {
        sources: root_sources,
        levels: Vec::new(),
        global_envelope_fields: 0,
        cumulative_witness_bits: 0,
        cumulative_matrix_work_fields: 0,
    }];

    for level in 0..config.offload_levels {
        let mut next_states = Vec::new();
        for state in &states {
            for plan in enumerate_level(&state.sources, level, &config, &mut oracle) {
                if state.sources.len() > 1 {
                    let entering_recursive_bits = state.sources[0].bit_len;
                    if plan
                        .witness_bits
                        .saturating_mul(u128::from(config.min_offload_contraction))
                        > entering_recursive_bits
                    {
                        continue;
                    }
                }
                let mut next_sources = vec![plan.next_witness];
                if config.setup_offload
                    && level
                        < config
                            .setup_offload_levels
                            .unwrap_or(config.offload_levels as u32)
                            as usize
                {
                    next_sources.push(Source {
                        field_len: plan.prefix_fields,
                        value_bits: config.field_bits,
                        bit_len: plan
                            .prefix_fields
                            .saturating_mul(u128::from(config.field_bits)),
                        onehot_chunk_size: 0,
                    });
                }
                let mut levels = state.levels.clone();
                levels.push(plan.clone());
                next_states.push(State {
                    sources: next_sources,
                    levels,
                    global_envelope_fields: state.global_envelope_fields.max(plan.envelope_fields),
                    cumulative_witness_bits: state
                        .cumulative_witness_bits
                        .saturating_add(plan.witness_bits),
                    cumulative_matrix_work_fields: state
                        .cumulative_matrix_work_fields
                        .saturating_add(plan.matrix_work_fields),
                });
            }
        }
        states = prune_states(next_states);
        eprintln!(
            "level={} frontier={} rank_keys={} estimator_calls={}",
            level,
            states.len(),
            oracle.cache.len(),
            oracle.estimator_calls
        );
    }

    states.sort_by_key(state_sort_key);
    let field_bytes = u128::from(config.field_bits).div_ceil(8);
    println!(
        "field_bits,onehot_chunk_size,witness_chunks,witness_chunk_levels,tensor_levels,tensor_onehot_bound,subfield_embedding_norm,setup_offload,setup_offload_levels,min_offload_contraction,a_collision,slicing,slice_cap,global_envelope_fields,global_envelope_bytes,envelope_dominant,envelope_dominant_over_a,matrix_fields,matrix_bytes,slice_geometry,cumulative_witness_bits,cumulative_witness_bytes,cumulative_matrix_work_fields,final_witness_fields,final_prefix_fields,max_n_a,max_n_b,max_n_d,all_a_rank_one,levels"
    );
    let print_limit = if config.print_limit == 0 {
        states.len()
    } else {
        config.print_limit
    };
    for state in states.iter().take(print_limit) {
        let (max_n_a, max_n_b, max_n_d) = state_max_ranks(state);
        let (envelope_dominant, envelope_dominant_over_a) = state_envelope_dominance(state);
        let matrix_fields = format_matrix_sizes(&state.levels, 1);
        let matrix_bytes = format_matrix_sizes(&state.levels, field_bytes);
        let slice_geometry = format_slice_geometry(&state.levels);
        println!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},\"{}\",\"{}\",\"{}\",\"{}\",\"{}\",{},{},{},{},{},{},{},{},{},\"{}\"",
            config.field_bits,
            config.onehot_chunk_size,
            config.witness_chunks,
            config
                .witness_chunk_levels
                .unwrap_or(config.offload_levels as u32),
            config.tensor_levels,
            config.tensor_onehot_bound.name(),
            subfield_embedding_norm(config.field_bits),
            config.setup_offload,
            if config.setup_offload {
                config
                    .setup_offload_levels
                    .unwrap_or(config.offload_levels as u32)
            } else {
                0
            },
            config.min_offload_contraction,
            config.a_collision.name(),
            config.slicing.name(),
            config.max_slices_per_matrix,
            state.global_envelope_fields,
            state.global_envelope_fields.saturating_mul(field_bytes),
            envelope_dominant,
            envelope_dominant_over_a,
            matrix_fields,
            matrix_bytes,
            slice_geometry,
            state.cumulative_witness_bits,
            state.cumulative_witness_bits.div_ceil(8),
            state.cumulative_matrix_work_fields,
            state.sources.first().map_or(0, |source| source.field_len),
            state.sources.get(1).map_or(0, |source| source.field_len),
            max_n_a,
            max_n_b,
            max_n_d,
            max_n_a == 1,
            format_levels(&state.levels),
        );
    }
}

fn enumerate_level(
    sources: &[Source],
    level: usize,
    config: &Config,
    oracle: &mut RankOracle,
) -> Vec<LevelPlan> {
    let mut plans = Vec::new();
    let witness_chunks = if level
        < config
            .witness_chunk_levels
            .unwrap_or(config.offload_levels as u32) as usize
    {
        config.witness_chunks
    } else {
        1
    };
    for &d_a in &config.a_dims {
        let challenge = challenge_profile(d_a);
        for &d_b in &config.bd_dims {
            if d_b > d_a || !d_a.is_multiple_of(d_b) {
                continue;
            }
            for &d_d in &config.bd_dims {
                if d_d > d_b || !d_b.is_multiple_of(d_d) {
                    continue;
                }
                let dims = RoleDims {
                    a: d_a,
                    b: d_b,
                    d: d_d,
                };
                let checked_basis_min = if level < config.fixed_checked_basis_levels as usize {
                    config
                        .fixed_checked_basis
                        .unwrap_or(config.checked_basis_min)
                } else {
                    config.checked_basis_min
                };
                let checked_basis_max = if level < config.fixed_checked_basis_levels as usize {
                    config
                        .fixed_checked_basis
                        .unwrap_or(config.checked_basis_max)
                } else {
                    config.checked_basis_max
                };
                for log_open in checked_basis_min..=checked_basis_max {
                    for log_outer in checked_basis_min..=checked_basis_max {
                        let mut group_sets = Vec::with_capacity(sources.len());
                        let mut feasible = true;
                        for &source in sources {
                            let candidates = enumerate_group(
                                source,
                                dims,
                                log_outer,
                                log_open,
                                challenge,
                                level < config.tensor_levels as usize,
                                witness_chunks,
                                config,
                                oracle,
                            );
                            if candidates.is_empty() {
                                feasible = false;
                                break;
                            }
                            group_sets.push(candidates);
                        }
                        if !feasible {
                            continue;
                        }
                        let mut partials = vec![PartialLevel {
                            groups: Vec::new(),
                            max_a_fields: 0,
                            max_b_fields: 0,
                            d_logical_width: 0,
                            max_d_logical_width: 0,
                            next_fields: 0,
                            witness_bits: 0,
                            matrix_work_fields: 0,
                        }];
                        for group_set in group_sets {
                            let mut expanded = Vec::new();
                            for partial in &partials {
                                for group in &group_set {
                                    let mut groups = partial.groups.clone();
                                    groups.push(group.clone());
                                    expanded.push(PartialLevel {
                                        groups,
                                        max_a_fields: partial.max_a_fields.max(group.a_fields),
                                        max_b_fields: partial.max_b_fields.max(group.b_fields),
                                        d_logical_width: partial
                                            .d_logical_width
                                            .saturating_add(group.d_logical_width),
                                        max_d_logical_width: partial
                                            .max_d_logical_width
                                            .max(group.d_logical_width),
                                        next_fields: partial
                                            .next_fields
                                            .saturating_add(group.next_fields),
                                        witness_bits: partial
                                            .witness_bits
                                            .saturating_add(group.witness_bits),
                                        matrix_work_fields: partial
                                            .matrix_work_fields
                                            .saturating_add(group.matrix_work_fields),
                                    });
                                }
                            }
                            partials = prune_partials(expanded);
                        }
                        for partial in partials {
                            let d_value_run = u128::from(d_a / d_d).saturating_mul(u128::from(
                                full_field_digits(config.field_bits, log_open),
                            ));
                            let Some(d_plan) = fit_sliced_matrix(
                                SisMatrixRole::D,
                                d_d,
                                partial.max_d_logical_width,
                                d_value_run,
                                gadget_bound(log_open),
                                partial.max_a_fields,
                                config.slicing,
                                oracle,
                            ) else {
                                continue;
                            };
                            let d_slices: u128 = partial
                                .groups
                                .iter()
                                .map(|group| {
                                    let values = group.d_logical_width / d_value_run;
                                    values.div_ceil(d_plan.values_per_slice)
                                })
                                .sum();
                            if config.max_slices_per_matrix != 0
                                && d_slices > config.max_slices_per_matrix
                            {
                                continue;
                            }
                            // First H-layer witness: decompose every coefficient of
                            // the complete stacked D image.  Later H layers are a
                            // deliberately separate refinement of this experiment.
                            let d_compression_fields = d_slices
                                .saturating_mul(u128::from(d_plan.rank))
                                .saturating_mul(u128::from(d_d))
                                .saturating_mul(u128::from(full_field_digits(
                                    config.field_bits,
                                    log_open,
                                )));
                            let d_compression_bits =
                                d_compression_fields.saturating_mul(u128::from(log_open));
                            let b_compression_fields = partial
                                .groups
                                .iter()
                                .map(|group| group.b_compression_fields)
                                .sum::<u128>();
                            let d_fields = d_plan.fields;
                            let envelope_fields =
                                partial.max_a_fields.max(partial.max_b_fields).max(d_fields);
                            let prefix_fields = next_power_of_two(envelope_fields);
                            plans.push(LevelPlan {
                                dims,
                                log_outer,
                                log_open,
                                groups: partial.groups,
                                n_d: d_plan.rank,
                                d_rop_log2: d_plan.rop_log2,
                                a_fields: partial.max_a_fields,
                                b_fields: partial.max_b_fields,
                                d_fields,
                                d_logical_width: partial.d_logical_width,
                                d_physical_width: d_plan.physical_width,
                                d_slices,
                                compression_suffix_fields: b_compression_fields
                                    .saturating_add(d_compression_fields),
                                envelope_fields,
                                prefix_fields,
                                next_witness: Source {
                                    field_len: partial
                                        .next_fields
                                        .saturating_add(d_compression_fields),
                                    value_bits: log_outer.max(log_open),
                                    bit_len: partial
                                        .witness_bits
                                        .saturating_add(d_compression_bits),
                                    onehot_chunk_size: 0,
                                },
                                witness_bits: partial
                                    .witness_bits
                                    .saturating_add(d_compression_bits),
                                // Storage uses the one reusable physical slice,
                                // but the prover applies it to every logical D
                                // coefficient across all groups.
                                matrix_work_fields: partial.matrix_work_fields.saturating_add(
                                    partial
                                        .d_logical_width
                                        .saturating_mul(u128::from(d_plan.rank))
                                        .saturating_mul(u128::from(d_d)),
                                ),
                            });
                        }
                    }
                }
            }
        }
    }
    prune_levels(plans)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_group(
    source: Source,
    dims: RoleDims,
    log_outer: u32,
    log_open: u32,
    challenge: ChallengeProfile,
    tensor: bool,
    witness_chunks: u32,
    config: &Config,
    oracle: &mut RankOracle,
) -> Vec<GroupPlan> {
    let ring_elems = source.field_len.div_ceil(u128::from(dims.a));
    let max_position_exp = bit_length(next_power_of_two(ring_elems)).saturating_sub(1);
    let digits_outer = u128::from(full_field_digits(config.field_bits, log_outer));
    let digits_open = u128::from(full_field_digits(config.field_bits, log_open));
    let mut candidates = Vec::new();
    let source_bases: Vec<u32> = if source.onehot_chunk_size == 0 {
        (config.source_basis_min..=config.source_basis_max).collect()
    } else {
        vec![1]
    };
    for log_source in source_bases {
        let (digits_source, source_linf) = if source.onehot_chunk_size == 0 {
            let digits_source = u128::from(num_digits_for_bound(
                source.value_bits,
                config.field_bits,
                log_source,
            ));
            let digit_bits = source.value_bits.min(log_source);
            let Some(source_linf) = 1u128.checked_shl(digit_bits.saturating_sub(1)) else {
                continue;
            };
            (digits_source, source_linf)
        } else {
            // The root commits the original certified binary one-hot vector
            // directly. It has one digit plane and coefficient norm one;
            // there is no checked source decomposition basis to optimize.
            (1, 1)
        };
        for position_exp in 0..=max_position_exp {
            let Some(positions) = 1u128.checked_shl(position_exp) else {
                continue;
            };
            let blocks = ring_elems.div_ceil(positions);
            let a_width_u128 = positions.saturating_mul(digits_source);
            let num_fold_coeffs = a_width_u128.saturating_mul(u128::from(dims.a));
            let chunks = u128::from(witness_chunks);
            let max_chunk_blocks = blocks.div_ceil(chunks);
            let response_coeffs = num_fold_coeffs.saturating_mul(chunks);
            let ln_arg = 16u128.saturating_mul(response_coeffs).div_ceil(7);
            let ln_term = ceil_natural_log(ln_arg);
            let tensor_low_len = if tensor {
                optimal_tensor_low_len(blocks)
            } else {
                0
            };
            let effective_challenge_l1 = if tensor {
                challenge.l1.saturating_mul(challenge.l1)
            } else {
                challenge.l1
            };
            let tstar = if tensor {
                let witness_linf_sq = source_linf.saturating_mul(source_linf);
                let sparse_ratio = if source.onehot_chunk_size > 0
                    && config.tensor_onehot_bound == TensorOnehotBound::SparseProxy
                {
                    Some((
                        u128::from(dims.a).div_ceil(u128::from(source.onehot_chunk_size)),
                        u128::from(dims.a),
                    ))
                } else {
                    None
                };
                let tail = tensor_fold_cap(
                    blocks,
                    chunks,
                    tensor_low_len,
                    response_coeffs,
                    challenge,
                    witness_linf_sq,
                    sparse_ratio,
                );
                let witness_l1 = if source.onehot_chunk_size > 0 {
                    u128::from(dims.a)
                        .div_ceil(u128::from(source.onehot_chunk_size))
                        .saturating_mul(source_linf)
                } else {
                    u128::from(dims.a).saturating_mul(source_linf)
                };
                let beta_per_block = challenge
                    .l1
                    .saturating_mul(challenge_infinity_norm(challenge))
                    .saturating_mul(witness_l1)
                    .min(effective_challenge_l1.saturating_mul(source_linf));
                tail.min(max_chunk_blocks.saturating_mul(beta_per_block))
            } else if source.onehot_chunk_size == 0 {
                let tstar_sq = BigUint::from(2u8)
                    * BigUint::from(max_chunk_blocks)
                    * BigUint::from(challenge.l2_sq)
                    * BigUint::from(source_linf)
                    * BigUint::from(source_linf)
                    * BigUint::from(ln_term);
                let tstar_big = sqrt_ceil_big(&tstar_sq);
                let Some(tstar) = tstar_big.to_u128() else {
                    continue;
                };
                tstar
            } else {
                let nonzeros_per_ring =
                    u128::from(dims.a).div_ceil(u128::from(source.onehot_chunk_size));
                let beta_per_block = challenge
                    .l1
                    .min(challenge_infinity_norm(challenge).saturating_mul(nonzeros_per_ring));
                let beta = max_chunk_blocks.saturating_mul(beta_per_block);
                exact_onehot_chernoff_cap(
                    max_chunk_blocks,
                    response_coeffs,
                    nonzeros_per_ring,
                    challenge,
                )
                .min(beta)
            };
            let (digits_fold, honest_z_cap, certified_z_difference) =
                snapped_fold_digit_plan(tstar, config.field_bits, log_open);
            let embedding_norm = u128::from(subfield_embedding_norm(config.field_bits));
            let a_bound_big = match config.a_collision {
                // Experimental control: the old honest-response heuristic also
                // pays the naive response-difference factor of two.
                ACollisionMode::HonestCap => {
                    BigUint::from(honest_z_cap)
                        * BigUint::from(
                            8u128
                                .saturating_mul(effective_challenge_l1)
                                .saturating_mul(embedding_norm),
                        )
                }
                // The certified digit interval is already a difference
                // envelope, so it must not be doubled again.
                ACollisionMode::CertifiedDifference => {
                    BigUint::from(certified_z_difference)
                        * BigUint::from(
                            4u128
                                .saturating_mul(effective_challenge_l1)
                                .saturating_mul(embedding_norm),
                        )
                }
            };
            let Some(a_bound) = a_bound_big.to_u128() else {
                continue;
            };
            let Ok(a_width) = u64::try_from(a_width_u128) else {
                continue;
            };
            let Some(a_rank) = oracle.min_rank(RankKey {
                role: SisMatrixRole::A,
                d: dims.a,
                width: a_width,
                bound: a_bound,
            }) else {
                continue;
            };
            let b_value_run = u128::from(dims.a / dims.b).saturating_mul(digits_outer);
            let b_logical_width = blocks
                .saturating_mul(u128::from(a_rank.rank))
                .saturating_mul(b_value_run);
            let Some(b_plan) = fit_sliced_matrix(
                SisMatrixRole::B,
                dims.b,
                b_logical_width,
                b_value_run,
                gadget_bound(log_outer),
                a_width_u128
                    .saturating_mul(u128::from(a_rank.rank))
                    .saturating_mul(u128::from(dims.a)),
                config.slicing,
                oracle,
            ) else {
                continue;
            };
            let digits_fold = u128::from(digits_fold);
            let e_fields = blocks
                .saturating_mul(digits_open)
                .saturating_mul(u128::from(dims.a));
            let t_fields = blocks
                .saturating_mul(u128::from(a_rank.rank))
                .saturating_mul(digits_outer)
                .saturating_mul(u128::from(dims.a));
            let z_fields = positions
                .saturating_mul(digits_source)
                .saturating_mul(digits_fold)
                .saturating_mul(u128::from(dims.a))
                .saturating_mul(chunks);
            let mut next_fields = e_fields.saturating_add(t_fields).saturating_add(z_fields);
            let mut witness_bits = e_fields
                .saturating_mul(u128::from(log_open))
                .saturating_add(t_fields.saturating_mul(u128::from(log_outer)))
                .saturating_add(z_fields.saturating_mul(u128::from(log_open)));
            let a_fields = a_width_u128
                .saturating_mul(u128::from(a_rank.rank))
                .saturating_mul(u128::from(dims.a));
            let b_slices =
                (blocks.saturating_mul(u128::from(a_rank.rank))).div_ceil(b_plan.values_per_slice);
            if config.max_slices_per_matrix != 0 && b_slices > config.max_slices_per_matrix {
                continue;
            }
            // First F-layer witness: decompose every coefficient of the complete
            // stacked B image.  This is the cost that prevents free over-slicing.
            let b_compression_fields = b_slices
                .saturating_mul(u128::from(b_plan.rank))
                .saturating_mul(u128::from(dims.b))
                .saturating_mul(u128::from(full_field_digits(config.field_bits, log_outer)));
            next_fields = next_fields.saturating_add(b_compression_fields);
            witness_bits = witness_bits
                .saturating_add(b_compression_fields.saturating_mul(u128::from(log_outer)));
            let d_logical_width = blocks
                .saturating_mul(u128::from(dims.a / dims.d))
                .saturating_mul(digits_open);
            candidates.push(GroupPlan {
                source,
                log_source,
                log_outer,
                log_open,
                positions,
                blocks,
                digits_source,
                digits_outer,
                digits_open,
                digits_fold,
                unsnapped_fold_cap: tstar,
                honest_fold_cap: honest_z_cap,
                tensor_low_len,
                n_a: a_rank.rank,
                n_b: b_plan.rank,
                a_rop_log2: a_rank.rop_log2,
                b_rop_log2: b_plan.rop_log2,
                a_fields,
                b_fields: b_plan.fields,
                b_logical_width,
                b_physical_width: b_plan.physical_width,
                b_slices,
                b_compression_fields,
                d_logical_width,
                next_fields,
                witness_bits,
                // B storage is the reusable physical slice. Work charges all
                // applications over the complete logical vector.
                matrix_work_fields: a_fields.saturating_add(
                    b_logical_width
                        .saturating_mul(u128::from(b_plan.rank))
                        .saturating_mul(u128::from(dims.b)),
                ),
            });
        }
    }
    prune_groups(candidates)
}

fn prune_groups(mut candidates: Vec<GroupPlan>) -> Vec<GroupPlan> {
    candidates.sort_by_key(group_sort_key);
    pareto_insert(candidates, group_dominates)
}

fn prune_partials(mut candidates: Vec<PartialLevel>) -> Vec<PartialLevel> {
    candidates.sort_by_key(partial_sort_key);
    pareto_insert(candidates, partial_dominates)
}

fn prune_levels(mut candidates: Vec<LevelPlan>) -> Vec<LevelPlan> {
    candidates.sort_by_key(level_sort_key);
    pareto_insert(candidates, level_dominates)
}

fn prune_states(mut candidates: Vec<State>) -> Vec<State> {
    candidates.sort_by_key(state_sort_key);
    pareto_insert(candidates, state_dominates)
}

fn pareto_insert<T>(items: Vec<T>, dominates: fn(&T, &T) -> bool) -> Vec<T> {
    let mut frontier = Vec::new();
    for item in items {
        if frontier.iter().any(|current| dominates(current, &item)) {
            continue;
        }
        frontier.retain(|current| !dominates(&item, current));
        frontier.push(item);
    }
    frontier
}

fn group_dominates(a: &GroupPlan, b: &GroupPlan) -> bool {
    le_all(
        &[
            a.a_fields,
            a.b_fields,
            a.d_logical_width,
            a.next_fields,
            a.witness_bits,
            a.matrix_work_fields,
        ],
        &[
            b.a_fields,
            b.b_fields,
            b.d_logical_width,
            b.next_fields,
            b.witness_bits,
            b.matrix_work_fields,
        ],
    )
}

fn partial_dominates(a: &PartialLevel, b: &PartialLevel) -> bool {
    le_all(
        &[
            a.max_a_fields,
            a.max_b_fields,
            a.d_logical_width,
            a.max_d_logical_width,
            a.next_fields,
            a.witness_bits,
            a.matrix_work_fields,
        ],
        &[
            b.max_a_fields,
            b.max_b_fields,
            b.d_logical_width,
            b.max_d_logical_width,
            b.next_fields,
            b.witness_bits,
            b.matrix_work_fields,
        ],
    )
}

fn level_dominates(a: &LevelPlan, b: &LevelPlan) -> bool {
    le_all(
        &[
            a.envelope_fields,
            a.prefix_fields,
            a.next_witness.field_len,
            a.next_witness.bit_len,
            u128::from(a.next_witness.value_bits),
            a.witness_bits,
            a.matrix_work_fields,
        ],
        &[
            b.envelope_fields,
            b.prefix_fields,
            b.next_witness.field_len,
            b.next_witness.bit_len,
            u128::from(b.next_witness.value_bits),
            b.witness_bits,
            b.matrix_work_fields,
        ],
    )
}

fn state_dominates(a: &State, b: &State) -> bool {
    if a.sources.len() != b.sources.len() {
        return false;
    }
    let source_le = a.sources.iter().zip(&b.sources).all(|(lhs, rhs)| {
        lhs.field_len <= rhs.field_len
            && lhs.bit_len <= rhs.bit_len
            && lhs.value_bits <= rhs.value_bits
    });
    source_le
        && le_all(
            &[
                a.global_envelope_fields,
                a.cumulative_witness_bits,
                a.cumulative_matrix_work_fields,
            ],
            &[
                b.global_envelope_fields,
                b.cumulative_witness_bits,
                b.cumulative_matrix_work_fields,
            ],
        )
}

fn le_all(a: &[u128], b: &[u128]) -> bool {
    a.iter().zip(b).all(|(lhs, rhs)| lhs <= rhs)
}

fn group_sort_key(plan: &GroupPlan) -> (u128, u128, u128, u128) {
    (
        plan.a_fields.max(plan.b_fields),
        plan.witness_bits,
        plan.next_fields,
        plan.matrix_work_fields,
    )
}

fn partial_sort_key(plan: &PartialLevel) -> (u128, u128, u128, u128) {
    (
        plan.max_a_fields.max(plan.max_b_fields),
        plan.witness_bits,
        plan.next_fields,
        plan.matrix_work_fields,
    )
}

fn level_sort_key(plan: &LevelPlan) -> (u128, u128, u128, u128) {
    (
        plan.envelope_fields,
        plan.witness_bits,
        plan.next_witness.field_len,
        plan.matrix_work_fields,
    )
}

fn state_sort_key(state: &State) -> (u128, u128, u128, u128) {
    (
        state.global_envelope_fields,
        state.cumulative_witness_bits,
        state.cumulative_matrix_work_fields,
        state.sources.iter().map(|source| source.field_len).sum(),
    )
}

fn state_max_ranks(state: &State) -> (u32, u32, u32) {
    let max_n_a = state
        .levels
        .iter()
        .flat_map(|level| &level.groups)
        .map(|group| group.n_a)
        .max()
        .unwrap_or(0);
    let max_n_b = state
        .levels
        .iter()
        .flat_map(|level| &level.groups)
        .map(|group| group.n_b)
        .max()
        .unwrap_or(0);
    let max_n_d = state
        .levels
        .iter()
        .map(|level| level.n_d)
        .max()
        .unwrap_or(0);
    (max_n_a, max_n_b, max_n_d)
}

fn level_envelope_dominance(level: &LevelPlan) -> (String, String) {
    let mut roles = Vec::new();
    let mut ratios = Vec::new();
    for (role, fields) in [
        ("A", level.a_fields),
        ("B", level.b_fields),
        ("D", level.d_fields),
    ] {
        if fields != level.envelope_fields {
            continue;
        }
        roles.push(role);
        if role != "A" {
            let decimal = if level.a_fields == 0 {
                "undefined".to_string()
            } else {
                format!("{:.6}", fields as f64 / level.a_fields as f64)
            };
            ratios.push(format!("{role}/A={fields}/{}={decimal}", level.a_fields));
        }
    }
    (roles.join("+"), ratios.join(";"))
}

fn state_envelope_dominance(state: &State) -> (String, String) {
    let mut dominators = Vec::new();
    let mut ratios = Vec::new();
    for (level_index, level) in state.levels.iter().enumerate() {
        if level.envelope_fields != state.global_envelope_fields {
            continue;
        }
        let (roles, level_ratios) = level_envelope_dominance(level);
        dominators.push(format!("L{level_index}:{roles}"));
        if !level_ratios.is_empty() {
            ratios.push(format!("L{level_index}:{level_ratios}"));
        }
    }
    (dominators.join("|"), ratios.join("|"))
}

fn format_matrix_sizes(levels: &[LevelPlan], scale: u128) -> String {
    levels
        .iter()
        .enumerate()
        .map(|(level_index, level)| {
            let a = level
                .groups
                .iter()
                .enumerate()
                .map(|(group_index, group)| {
                    format!("g{group_index}={}", group.a_fields.saturating_mul(scale))
                })
                .collect::<Vec<_>>()
                .join(";");
            let b = level
                .groups
                .iter()
                .enumerate()
                .map(|(group_index, group)| {
                    format!("g{group_index}={}", group.b_fields.saturating_mul(scale))
                })
                .collect::<Vec<_>>()
                .join(";");
            format!(
                "L{level_index}:A[{a}],B[{b}],D={}",
                level.d_fields.saturating_mul(scale)
            )
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn format_slice_geometry(levels: &[LevelPlan]) -> String {
    levels
        .iter()
        .enumerate()
        .map(|(level_index, level)| {
            let b = level
                .groups
                .iter()
                .enumerate()
                .map(|(group_index, group)| {
                    format!(
                        "g{group_index}={}/{}x{}",
                        group.b_physical_width, group.b_logical_width, group.b_slices
                    )
                })
                .collect::<Vec<_>>()
                .join(";");
            format!(
                "L{level_index}:B[{b}],D={}/{}x{}",
                level.d_physical_width, level.d_logical_width, level.d_slices
            )
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn format_levels(levels: &[LevelPlan]) -> String {
    levels
        .iter()
        .enumerate()
        .map(|(level_index, level)| {
            let (dominant, dominant_over_a) = level_envelope_dominance(level);
            let dominance = if dominant_over_a.is_empty() {
                format!("dom={dominant}")
            } else {
                format!("dom={dominant}({dominant_over_a})")
            };
            let groups = level
                .groups
                .iter()
                .enumerate()
                .map(|(group_index, group)| {
                    format!(
                        "g{group_index}:src={}b{}k{} ls/lt/lo={}/{}/{} B={} P={} tensorL={} da={} db={} na={} nb={} ds/dt/de/dz={}/{}/{}/{} zcap={}->{} A/B={}/{} Bslice={}/{}x{} Bcmp={} secA/B={:.2}/{:.2}",
                        group.source.field_len,
                        group.source.value_bits,
                        group.source.onehot_chunk_size,
                        group.log_source,
                        group.log_outer,
                        group.log_open,
                        group.blocks,
                        group.positions,
                        group.tensor_low_len,
                        level.dims.a,
                        level.dims.b,
                        group.n_a,
                        group.n_b,
                        group.digits_source,
                        group.digits_outer,
                        group.digits_open,
                        group.digits_fold,
                        group.unsnapped_fold_cap,
                        group.honest_fold_cap,
                        group.a_fields,
                        group.b_fields,
                        group.b_physical_width,
                        group.b_logical_width,
                        group.b_slices,
                        group.b_compression_fields,
                        group.a_rop_log2,
                        group.b_rop_log2,
                    )
                })
                .collect::<Vec<_>>()
                .join(";");
            format!(
                "L{level_index}[dims={}/{}/{} lt/lo={}/{} nd={} A/B/D={}/{}/{} Dslice={}/{}x{} cmp={} env={} {} prefix={} next={} bits={} work={} secD={:.2};{}]",
                level.dims.a,
                level.dims.b,
                level.dims.d,
                level.log_outer,
                level.log_open,
                level.n_d,
                level.a_fields,
                level.b_fields,
                level.d_fields,
                level.d_physical_width,
                level.d_logical_width,
                level.d_slices,
                level.compression_suffix_fields,
                level.envelope_fields,
                dominance,
                level.prefix_fields,
                level.next_witness.field_len,
                level.witness_bits,
                level.matrix_work_fields,
                level.d_rop_log2,
                groups,
            )
        })
        .collect::<Vec<_>>()
        .join(" | ")
}

fn challenge_profile(d: u32) -> ChallengeProfile {
    let (count_pm1, count_pm2) = match d {
        64 => (31, 10),
        128 => (31, 0),
        256 => (23, 0),
        512 => (19, 0),
        1024 => (16, 0),
        2048 => (14, 0),
        _ => (minimum_pm1_weight_for_128_bits(d), 0),
    };
    ChallengeProfile {
        ring_d: d,
        count_pm1,
        count_pm2,
        l1: u128::from(count_pm1) + 2 * u128::from(count_pm2),
        l2_sq: u128::from(count_pm1) + 4 * u128::from(count_pm2),
    }
}

fn challenge_infinity_norm(challenge: ChallengeProfile) -> u128 {
    if challenge.count_pm2 == 0 {
        1
    } else {
        2
    }
}

fn optimal_tensor_low_len(blocks: u128) -> u128 {
    let capacity = blocks.max(1).next_power_of_two();
    let mut low = 1u128;
    let mut best = (u128::MAX, 1u128);
    loop {
        let work = low.saturating_add(blocks.div_ceil(low));
        if (work, low) < best {
            best = (work, low);
        }
        if low == capacity {
            break;
        }
        low = low.saturating_mul(2);
    }
    best.1
}

fn max_chunk_tensor_high_len(blocks: u128, chunks: u128, low_len: u128) -> u128 {
    let base = blocks / chunks;
    let extra = blocks % chunks;
    let mut start = 0u128;
    let mut max_high = 0u128;
    for chunk in 0..chunks {
        let len = base + u128::from(chunk < extra);
        if len > 0 {
            let end = start + len;
            let first_high = start / low_len;
            let last_high = (end - 1) / low_len;
            max_high = max_high.max(last_high - first_high + 1);
            start = end;
        }
    }
    max_high
}

/// Two-stage tensor-chaos cap. `sparse_ratio = Some(q, D)` is the proposed
/// one-hot refinement; the generic path matches the certified dense tensor
/// envelope and deliberately does not claim the sparsity improvement.
fn tensor_fold_cap(
    blocks: u128,
    chunks: u128,
    low_len: u128,
    response_coeffs: u128,
    challenge: ChallengeProfile,
    witness_linf_sq: u128,
    sparse_ratio: Option<(u128, u128)>,
) -> u128 {
    let high_len = max_chunk_tensor_high_len(blocks, chunks, low_len);
    let factor_l2_sq = challenge.l2_sq;
    let factor_support = u128::from(challenge.count_pm1 + challenge.count_pm2);
    if let Some((sparse_num, sparse_den)) = sparse_ratio {
        let outer_ln = (32.0 * response_coeffs as f64 / 7.0).ln();
        let high_ln =
            (32.0 * response_coeffs as f64 * high_len as f64 * factor_support as f64 / 7.0).ln();
        let low_ln =
            (32.0 * response_coeffs as f64 * low_len as f64 * factor_support as f64 / 7.0).ln();
        let variance = 4.0
            * high_len as f64
            * low_len as f64
            * factor_l2_sq as f64
            * factor_l2_sq as f64
            * witness_linf_sq as f64
            * outer_ln
            * high_ln.min(low_ln)
            * sparse_num as f64
            / sparse_den as f64;
        return (variance.sqrt() * (1.0 + 1e-12) + 1e-9).ceil() as u128;
    }

    let outer_ln = ceil_natural_log(32u128.saturating_mul(response_coeffs).div_ceil(7));
    let high_ln = ceil_natural_log(
        32u128
            .saturating_mul(response_coeffs)
            .saturating_mul(high_len)
            .saturating_mul(factor_support)
            .div_ceil(7),
    );
    let low_ln = ceil_natural_log(
        32u128
            .saturating_mul(response_coeffs)
            .saturating_mul(low_len)
            .saturating_mul(factor_support)
            .div_ceil(7),
    );
    let variance = BigUint::from(4u8)
        * BigUint::from(high_len)
        * BigUint::from(low_len)
        * BigUint::from(factor_l2_sq)
        * BigUint::from(factor_l2_sq)
        * BigUint::from(witness_linf_sq)
        * BigUint::from(outer_ln)
        * BigUint::from(high_ln.min(low_ln));
    sqrt_ceil_big(&variance).to_u128().unwrap_or(u128::MAX)
}

fn log_binomial(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    let k = k.min(n - k);
    (1..=k)
        .map(|i| (f64::from(n - k + i) / f64::from(i)).ln())
        .sum()
}

fn log_cosh(value: f64) -> f64 {
    let abs = value.abs();
    abs + (-2.0 * abs).exp().ln_1p() - std::f64::consts::LN_2
}

/// Exact one-block log-MGF and its derivative for a fixed output coordinate
/// of a certified one-hot ring. The challenge's fixed-weight support makes the
/// numbers of ±1 and ±2 coefficients hit by the witness a multivariate
/// hypergeometric draw; random signs then contribute powers of cosh.
fn onehot_log_mgf_and_derivative(
    lambda: f64,
    nonzeros_per_ring: u32,
    challenge: ChallengeProfile,
) -> (f64, f64) {
    let d = challenge.ring_d;
    let q = nonzeros_per_ring.min(d);
    let denominator = log_binomial(d, challenge.count_pm1)
        + log_binomial(d - challenge.count_pm1, challenge.count_pm2);
    let log_cosh_one = log_cosh(lambda);
    let log_cosh_two = log_cosh(2.0 * lambda);
    let mut terms = Vec::new();
    for hit_pm1 in 0..=q.min(challenge.count_pm1) {
        for hit_pm2 in 0..=(q - hit_pm1).min(challenge.count_pm2) {
            let outside_pm1 = challenge.count_pm1 - hit_pm1;
            let outside_pm2 = challenge.count_pm2 - hit_pm2;
            if outside_pm1 > d - q {
                continue;
            }
            let outside_after_pm1 = d - q - outside_pm1;
            if outside_pm2 > outside_after_pm1 {
                continue;
            }
            let log_probability = log_binomial(q, hit_pm1)
                + log_binomial(q - hit_pm1, hit_pm2)
                + log_binomial(d - q, outside_pm1)
                + log_binomial(outside_after_pm1, outside_pm2)
                - denominator;
            let log_term = log_probability
                + f64::from(hit_pm1) * log_cosh_one
                + f64::from(hit_pm2) * log_cosh_two;
            let derivative = f64::from(hit_pm1) * lambda.tanh()
                + 2.0 * f64::from(hit_pm2) * (2.0 * lambda).tanh();
            terms.push((log_term, derivative));
        }
    }
    let max_log = terms
        .iter()
        .map(|(log_term, _)| *log_term)
        .fold(f64::NEG_INFINITY, f64::max);
    let scaled_sum: f64 = terms
        .iter()
        .map(|(log_term, _)| (*log_term - max_log).exp())
        .sum();
    let derivative = terms
        .iter()
        .map(|(log_term, derivative)| (*log_term - max_log).exp() * derivative)
        .sum::<f64>()
        / scaled_sum;
    (max_log + scaled_sum.ln(), derivative)
}

/// Smallest conservative integer cap obtained from the exact signed-sparse
/// one-hot Chernoff exponent, followed by the caller's clean-digit snap.
///
/// For one block, `M(lambda)` is exact. Independent block challenges give
/// `M(lambda)^blocks`; a two-sided union bound over all folded coefficients and
/// the shipped 1/8 grind target requires
/// `2*N*exp(B log M(lambda) - lambda*t) <= 7/8`.
fn exact_onehot_chernoff_cap(
    blocks: u128,
    num_fold_coeffs: u128,
    nonzeros_per_ring: u128,
    challenge: ChallengeProfile,
) -> u128 {
    let blocks_f = blocks as f64;
    let union_log = (16.0 * num_fold_coeffs as f64 / 7.0).ln();
    let q = u32::try_from(nonzeros_per_ring).unwrap_or(u32::MAX);

    // The derivative of g(lambda)=(B log M(lambda)+L)/lambda changes sign
    // exactly when lambda*B*M'/M - B*log(M) - L does. This expression is
    // monotone because its derivative is lambda*B times the tilted variance.
    let stationarity = |lambda: f64| {
        let (log_mgf, derivative) = onehot_log_mgf_and_derivative(lambda, q, challenge);
        lambda * blocks_f * derivative - blocks_f * log_mgf - union_log
    };
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    while stationarity(hi) < 0.0 && hi < 64.0 {
        hi *= 2.0;
    }
    for _ in 0..96 {
        let mid = (lo + hi) * 0.5;
        if stationarity(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let lambda = hi.max(f64::EPSILON);
    let (log_mgf, _) = onehot_log_mgf_and_derivative(lambda, q, challenge);
    let threshold = (blocks_f * log_mgf + union_log) / lambda;
    // Float evaluation is used only by this offline discovery tool. Bias the
    // final ceiling upward so roundoff cannot select a smaller integer cap.
    (threshold * (1.0 + 1e-12) + 1e-9).ceil() as u128
}

fn minimum_pm1_weight_for_128_bits(d: u32) -> u32 {
    for weight in 1..=d {
        let mut support_bits = f64::from(weight);
        for i in 1..=weight {
            support_bits += (f64::from(d - weight + i) / f64::from(i)).log2();
        }
        if support_bits >= 128.0 {
            return weight;
        }
    }
    d
}

fn secure(cost: CostValue) -> bool {
    match cost {
        CostValue::Finite(cost) => cost.log2 >= TARGET_BITS,
        CostValue::ProvenAboveTarget(lower) => lower.log2 >= TARGET_BITS,
        CostValue::Infinity => true,
    }
}

fn gadget_bound(log_basis: u32) -> u128 {
    1u128
        .checked_shl(log_basis)
        .unwrap_or(u128::MAX)
        .saturating_sub(1)
}

fn full_field_digits(field_bits: u32, log_basis: u32) -> u32 {
    field_bits.div_ceil(log_basis).max(1)
}

const fn subfield_embedding_norm(field_bits: u32) -> u32 {
    match field_bits {
        128 => 1,
        32 | 64 => 2,
        _ => panic!("unsupported Akita field width"),
    }
}

fn num_digits_for_bound(log_bound: u32, field_bits: u32, log_basis: u32) -> u32 {
    if log_bound >= field_bits {
        return full_field_digits(field_bits, log_basis);
    }
    let mut digits = log_bound.div_ceil(log_basis).max(1);
    let required_positive = if log_bound == 0 {
        0
    } else {
        (1u128 << (log_bound - 1)) - 1
    };
    if balanced_digit_max(log_basis, digits) < required_positive {
        digits += 1;
    }
    digits
}

fn balanced_digit_max(log_basis: u32, digits: u32) -> u128 {
    let Some(base) = 1u128.checked_shl(log_basis) else {
        return u128::MAX;
    };
    let max_digit = base / 2 - 1;
    let mut power = 1u128;
    for _ in 0..digits {
        power = power.saturating_mul(base);
    }
    max_digit.saturating_mul(power.saturating_sub(1) / (base - 1))
}

fn certified_digit_difference_envelope(
    log_digit_basis: u32,
    log_recomposition_basis: u32,
    digits: u32,
) -> u128 {
    let Some(digit_base) = 1u128.checked_shl(log_digit_basis) else {
        return u128::MAX;
    };
    let Some(recomposition_base) = 1u128.checked_shl(log_recomposition_basis) else {
        return u128::MAX;
    };
    let mut power = 1u128;
    for _ in 0..digits.max(1) {
        power = power.saturating_mul(recomposition_base);
    }
    // Two balanced base-b' digits differ by at most b'-1. Recomposition
    // in base b therefore has exact certified diameter
    // (b'-1)(b^k-1)/(b-1).
    digit_base
        .saturating_sub(1)
        .saturating_mul(power.saturating_sub(1) / recomposition_base.saturating_sub(1))
}

fn snapped_fold_digit_plan(tstar: u128, field_bits: u32, log_basis: u32) -> (u32, u128, u128) {
    let mut digits =
        num_digits_for_bound(bit_length(tstar).saturating_add(1), field_bits, log_basis);
    let floor = (tstar / 2).max(1);
    let mut honest_cap = tstar;
    while digits > 1 {
        let positive_lower = balanced_digit_max(log_basis, digits - 1);
        if positive_lower < floor {
            break;
        }
        digits -= 1;
        honest_cap = honest_cap.min(positive_lower);
    }
    (
        digits,
        honest_cap,
        certified_digit_difference_envelope(log_basis, log_basis, digits),
    )
}

fn ceil_natural_log(x: u128) -> u128 {
    if x <= 1 {
        return 0;
    }
    let ceil_log2 = u128::from(128u32.saturating_sub((x - 1).leading_zeros()));
    ceil_log2.saturating_mul(71).div_ceil(100)
}

fn sqrt_ceil_big(value: &BigUint) -> BigUint {
    let root = value.sqrt();
    if &root * &root == *value {
        root
    } else {
        root + BigUint::from(1u8)
    }
}

fn next_power_of_two(value: u128) -> u128 {
    value
        .max(1)
        .checked_next_power_of_two()
        .unwrap_or(u128::MAX)
}

fn bit_length(value: u128) -> u32 {
    128 - value.leading_zeros()
}

fn parse_args() -> Config {
    let mut config = Config::default();
    let mut args = env::args().skip(1);
    while let Some(flag) = args.next() {
        if flag == "--help" || flag == "-h" {
            usage();
        }
        let value = args
            .next()
            .unwrap_or_else(|| panic!("missing value for {flag}"));
        match flag.as_str() {
            "--num-vars" => config.num_vars = parse(&value, &flag),
            "--sources" => config.initial_sources = Some(parse_sources(&value, &flag)),
            "--setup-offload" => config.setup_offload = parse(&value, &flag),
            "--setup-offload-levels" => config.setup_offload_levels = Some(parse(&value, &flag)),
            "--min-offload-contraction" => config.min_offload_contraction = parse(&value, &flag),
            "--a-collision" => {
                config.a_collision = match value.as_str() {
                    "honest" | "honest-cap" => ACollisionMode::HonestCap,
                    "certified" | "certified-digits" | "certified-difference" => {
                        ACollisionMode::CertifiedDifference
                    }
                    _ => panic!("invalid value {value} for {flag}"),
                }
            }
            "--levels" => config.offload_levels = parse(&value, &flag),
            "--field-bits" => config.field_bits = parse(&value, &flag),
            "--onehot-chunk-size" => config.onehot_chunk_size = parse(&value, &flag),
            "--witness-chunks" => config.witness_chunks = parse(&value, &flag),
            "--witness-chunk-levels" => config.witness_chunk_levels = Some(parse(&value, &flag)),
            "--tensor-levels" => config.tensor_levels = parse(&value, &flag),
            "--tensor-onehot-bound" => {
                config.tensor_onehot_bound = match value.as_str() {
                    "generic" => TensorOnehotBound::Generic,
                    "sparse" | "onehot-sparse-proxy" => TensorOnehotBound::SparseProxy,
                    _ => panic!("invalid value {value} for {flag}"),
                }
            }
            "--source-basis-min" => config.source_basis_min = parse(&value, &flag),
            "--source-basis-max" => config.source_basis_max = parse(&value, &flag),
            "--checked-basis-min" => config.checked_basis_min = parse(&value, &flag),
            "--checked-basis-max" => config.checked_basis_max = parse(&value, &flag),
            "--fixed-checked-basis" => config.fixed_checked_basis = Some(parse(&value, &flag)),
            "--fixed-checked-basis-levels" => {
                config.fixed_checked_basis_levels = parse(&value, &flag)
            }
            "--a-dims" => config.a_dims = parse_list(&value, &flag),
            "--bd-dims" => config.bd_dims = parse_list(&value, &flag),
            "--slicing" => {
                config.slicing = match value.as_str() {
                    "none" => SlicingMode::None,
                    "a-capped" | "a_capped" => SlicingMode::ACapped,
                    _ => panic!("invalid value {value} for {flag}"),
                }
            }
            "--max-slices-per-matrix" => config.max_slices_per_matrix = parse(&value, &flag),
            "--max-rank" => config.max_rank = parse(&value, &flag),
            "--print-limit" => config.print_limit = parse(&value, &flag),
            _ => panic!("unknown argument {flag}"),
        }
    }
    config
}

fn parse<T: std::str::FromStr>(value: &str, flag: &str) -> T {
    value
        .parse()
        .unwrap_or_else(|_| panic!("invalid value {value} for {flag}"))
}

fn parse_list(value: &str, flag: &str) -> Vec<u32> {
    value.split(',').map(|item| parse(item, flag)).collect()
}

fn parse_sources(value: &str, flag: &str) -> Vec<Source> {
    value
        .split(',')
        .map(|item| {
            let (field_len, value_bits) = item
                .split_once(':')
                .unwrap_or_else(|| panic!("invalid value {item} for {flag}; expected LEN:BITS"));
            let field_len = parse(field_len, flag);
            let value_bits = parse(value_bits, flag);
            Source {
                field_len,
                value_bits,
                bit_len: field_len.saturating_mul(u128::from(value_bits)),
                onehot_chunk_size: 0,
            }
        })
        .collect()
}

fn validate_config(config: &Config) {
    assert!(config.num_vars < 128);
    assert!(config.initial_sources.as_ref().is_none_or(|sources| {
        !sources.is_empty()
            && sources
                .iter()
                .all(|source| source.field_len > 0 && source.value_bits > 0)
    }));
    assert!(config.offload_levels > 0);
    assert!(config.min_offload_contraction > 0);
    assert!(config
        .setup_offload_levels
        .is_none_or(|levels| levels <= config.offload_levels as u32));
    assert!(matches!(config.field_bits, 32 | 64 | 128));
    assert!(config.onehot_chunk_size == 0 || config.onehot_chunk_size.is_power_of_two());
    assert!(config.witness_chunks.is_power_of_two());
    assert!(config
        .witness_chunk_levels
        .is_none_or(|levels| levels <= config.offload_levels as u32));
    assert!(config.tensor_levels <= config.offload_levels as u32);
    assert!(config.source_basis_min >= 2);
    assert!(config.source_basis_min <= config.source_basis_max);
    assert!(config.source_basis_max < 128);
    assert!(config.checked_basis_min >= 2);
    assert!(config.checked_basis_min <= config.checked_basis_max);
    assert!(config.checked_basis_max < 128);
    if let Some(basis) = config.fixed_checked_basis {
        assert!(basis >= 2);
        assert!(basis < 128);
    }
    assert!(config.fixed_checked_basis.is_some() || config.fixed_checked_basis_levels == 0);
    assert!(config.max_rank > 0);
    assert!(!config.a_dims.is_empty());
    assert!(!config.bd_dims.is_empty());
    assert!(config.a_dims.iter().all(|d| d.is_power_of_two()));
    assert!(config.bd_dims.iter().all(|d| d.is_power_of_two()));
}

fn usage() -> ! {
    println!(
        "usage: ideal_setup_offload_planner [--num-vars 30] [--sources LEN:BITS,...] [--setup-offload true] [--setup-offload-levels LEVELS] [--min-offload-contraction 1] [--a-collision certified-difference] [--levels 2] \
         [--field-bits 128] [--onehot-chunk-size 0] [--witness-chunks 1] [--witness-chunk-levels LEVELS] [--tensor-levels 0] \
         [--tensor-onehot-bound generic|sparse] [--source-basis-min 2] [--source-basis-max 32] \
         [--checked-basis-min 2] [--checked-basis-max 4] \
         [--fixed-checked-basis BASIS] [--fixed-checked-basis-levels LEVELS] \
         [--a-dims 64,128,256,512] [--bd-dims 16,32,64,128] \
         [--slicing none|a-capped] [--max-slices-per-matrix 0 (unbounded)] \
         [--max-rank 20] [--print-limit 0 (all)]"
    );
    std::process::exit(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mixed_role_envelope_is_max_of_flat_footprints() {
        let level = LevelPlan {
            dims: RoleDims {
                a: 256,
                b: 64,
                d: 32,
            },
            log_outer: 3,
            log_open: 3,
            groups: Vec::new(),
            n_d: 1,
            d_rop_log2: 130.0,
            a_fields: 100,
            b_fields: 90,
            d_fields: 80,
            d_logical_width: 10,
            d_physical_width: 10,
            d_slices: 1,
            compression_suffix_fields: 0,
            envelope_fields: 100,
            prefix_fields: 128,
            next_witness: Source {
                field_len: 10,
                value_bits: 3,
                bit_len: 30,
                onehot_chunk_size: 0,
            },
            witness_bits: 30,
            matrix_work_fields: 270,
        };
        assert_eq!(level.envelope_fields, 100);
        assert_eq!(level.prefix_fields, 128);
        assert_eq!(
            level_envelope_dominance(&level),
            ("A".into(), String::new())
        );

        let mut b_dominant = level.clone();
        b_dominant.b_fields = 150;
        b_dominant.envelope_fields = 150;
        assert_eq!(
            level_envelope_dominance(&b_dominant),
            ("B".into(), "B/A=150/100=1.500000".into())
        );

        let mut b_d_tie = b_dominant;
        b_d_tie.d_fields = 150;
        assert_eq!(
            level_envelope_dominance(&b_d_tie),
            (
                "B+D".into(),
                "B/A=150/100=1.500000;D/A=150/100=1.500000".into()
            )
        );
    }

    #[test]
    fn dominance_keeps_storage_proof_tradeoff() {
        let source = Source {
            field_len: 8,
            value_bits: 3,
            bit_len: 24,
            onehot_chunk_size: 0,
        };
        let base = State {
            sources: vec![source],
            levels: Vec::new(),
            global_envelope_fields: 100,
            cumulative_witness_bits: 100,
            cumulative_matrix_work_fields: 100,
        };
        let mut proof_better = base.clone();
        proof_better.global_envelope_fields = 110;
        proof_better.cumulative_witness_bits = 90;
        assert!(!state_dominates(&base, &proof_better));
        assert!(!state_dominates(&proof_better, &base));
    }

    #[test]
    fn recursive_witness_bound_collapses_source_digits() {
        assert_eq!(num_digits_for_bound(4, 128, 10), 1);
        assert_eq!(num_digits_for_bound(128, 128, 10), 13);
    }

    #[test]
    fn small_field_a_pricing_includes_trace_subfield_embedding_norm() {
        assert_eq!(subfield_embedding_norm(128), 1);
        assert_eq!(subfield_embedding_norm(64), 2);
        assert_eq!(subfield_embedding_norm(32), 2);
    }

    #[test]
    fn certified_difference_envelope_is_the_exact_recomposed_diameter() {
        assert_eq!(certified_digit_difference_envelope(2, 2, 10), (1 << 20) - 1);
        assert_eq!(certified_digit_difference_envelope(3, 3, 5), (1 << 15) - 1);
        assert_eq!(
            certified_digit_difference_envelope(2, 3, 5),
            3 * ((1 << 15) - 1) / 7
        );
    }

    #[test]
    fn snapped_fold_plan_retains_half_tstar_and_prices_certified_difference() {
        for log_basis in 2..=4 {
            for tstar in [17u128, 257, 4097, 65_537] {
                let (digits, honest_cap, certified_difference) =
                    snapped_fold_digit_plan(tstar, 128, log_basis);
                assert!(honest_cap >= tstar / 2);
                assert!(certified_difference >= honest_cap);
                if digits > 1 {
                    assert!(balanced_digit_max(log_basis, digits - 1) < tstar / 2);
                }
            }
        }
    }

    #[test]
    fn onehot_exact_mgf_matches_single_coordinate_distribution() {
        let challenge = challenge_profile(64);
        for lambda in [0.05, 0.5, 1.25] {
            let (actual, _) = onehot_log_mgf_and_derivative(lambda, 1, challenge);
            let expected = ((23.0 / 64.0)
                + (31.0 / 64.0) * lambda.cosh()
                + (10.0 / 64.0) * (2.0 * lambda).cosh())
            .ln();
            assert!((actual - expected).abs() < 1e-12);
        }
    }

    #[test]
    fn exact_onehot_tail_improves_generic_l2_bound_before_snap() {
        let blocks = 2_048u128;
        let num_fold_coeffs = 32_768u128 * 64;
        let challenge = challenge_profile(64);
        let exact = exact_onehot_chernoff_cap(blocks, num_fold_coeffs, 1, challenge);
        let ln_arg = 16u128.saturating_mul(num_fold_coeffs).div_ceil(7);
        let generic_sq = BigUint::from(2u8)
            * BigUint::from(blocks)
            * BigUint::from(challenge.l2_sq)
            * BigUint::from(ceil_natural_log(ln_arg));
        let generic = sqrt_ceil_big(&generic_sq).to_u128().unwrap();
        assert!(exact < generic);

        let (digits, honest_cap, _) = snapped_fold_digit_plan(exact, 128, 2);
        assert!(honest_cap >= exact / 2);
        if digits > 1 {
            assert!(balanced_digit_max(2, digits - 1) < exact / 2);
        }
    }

    #[test]
    fn tensor_shape_and_onehot_sparse_proxy_match_reference_case() {
        assert_eq!(optimal_tensor_low_len(1 << 15), 1 << 7);
        assert_eq!(optimal_tensor_low_len(1 << 16), 1 << 8);
        assert_eq!(max_chunk_tensor_high_len(1 << 16, 4, 1 << 8), 64);

        let challenge = challenge_profile(512);
        let response_coeffs = (1u128 << 19) * 512;
        let generic = tensor_fold_cap(1 << 15, 1, 1 << 8, response_coeffs, challenge, 1, None);
        let sparse = tensor_fold_cap(
            1 << 15,
            1,
            1 << 8,
            response_coeffs,
            challenge,
            1,
            Some((2, 512)),
        );
        assert_eq!(sparse, 10_541);
        assert!(generic > sparse);
    }

    #[test]
    fn role_ring_bridge_preserves_flat_witness_size() {
        let d_a = 512u128;
        let blocks = 604u128;
        let n_a = 1u128;
        let digits_outer = 32u128;
        let digits_open = 64u128;
        for d_b in [16u128, 32, 64, 128] {
            let logical_b = blocks * n_a * digits_outer * (d_a / d_b);
            assert_eq!(logical_b * d_b, blocks * n_a * digits_outer * d_a);
        }
        for d_d in [16u128, 32, 64, 128] {
            let logical_d = blocks * digits_open * (d_a / d_d);
            assert_eq!(logical_d * d_d, blocks * digits_open * d_a);
        }
    }

    #[test]
    fn higher_checked_basis_does_not_dominate_a_cheaper_future_bound() {
        let level = LevelPlan {
            dims: RoleDims {
                a: 256,
                b: 64,
                d: 32,
            },
            log_outer: 2,
            log_open: 2,
            groups: Vec::new(),
            n_d: 1,
            d_rop_log2: 130.0,
            a_fields: 100,
            b_fields: 90,
            d_fields: 80,
            d_logical_width: 10,
            d_physical_width: 10,
            d_slices: 1,
            compression_suffix_fields: 0,
            envelope_fields: 100,
            prefix_fields: 128,
            next_witness: Source {
                field_len: 10,
                value_bits: 2,
                bit_len: 30,
                onehot_chunk_size: 0,
            },
            witness_bits: 30,
            matrix_work_fields: 270,
        };
        let mut higher_basis = level.clone();
        higher_basis.next_witness.value_bits = 4;
        higher_basis.next_witness.field_len = 9;
        assert!(!level_dominates(&higher_basis, &level));
        assert!(!level_dominates(&level, &higher_basis));
    }
}
