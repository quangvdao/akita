//! Canonical batched prove/verify schedule resolution.

use crate::CommitmentConfig;
use akita_field::{AkitaError, ExtField, FieldCore, MulBaseUnreduced};
use akita_types::{
    folded_root_supports_opening_shape, root_direct_schedule, root_tensor_projection_enabled,
    schedule_root_fold_step, FpExtEncoding, OpeningClaimsLayout, Schedule,
};

/// Resolve the runtime schedule used by batched prove and verify.
///
/// Grouped roots (`num_groups() > 1`) are rewritten to root-direct with the
/// summed hypercube witness length. Scalar openings may fall back from a
/// fold-first catalog schedule to root-direct when the opening shape does not
/// support folded roots and tensor projection is disabled.
pub fn effective_batched_schedule<Cfg, const D: usize>(
    opening_batch: &OpeningClaimsLayout,
    opening_point: &[Cfg::ExtField],
) -> Result<Schedule, AkitaError>
where
    Cfg: CommitmentConfig,
    Cfg::Field: FieldCore,
    Cfg::ExtField: FpExtEncoding<Cfg::Field> + ExtField<Cfg::Field> + MulBaseUnreduced<Cfg::Field>,
{
    let num_vars = opening_batch.max_num_vars();
    let root_direct_witness_len = opening_batch.root_direct_witness_len()?;
    let mut schedule = Cfg::get_params_for_prove(opening_batch)?;
    if opening_batch.num_groups() > 1 {
        let commit_params = Cfg::grouped_root_commit_params(&schedule)?;
        schedule = root_direct_schedule(root_direct_witness_len, commit_params)?;
    }
    if let Some(root_step) = schedule_root_fold_step(&schedule) {
        let alpha_bits = root_step.params.ring_dimension.trailing_zeros() as usize;
        if !folded_root_supports_opening_shape::<Cfg::Field, Cfg::ExtField, D>(
            std::slice::from_ref(&opening_point),
            &root_step.params,
            alpha_bits,
        ) && !root_tensor_projection_enabled::<Cfg::Field, Cfg::ExtField, D>(num_vars)
        {
            let commit_params = Cfg::grouped_root_commit_params(&schedule)?;
            schedule = root_direct_schedule(root_direct_witness_len, commit_params)?;
        }
    }

    Ok(schedule)
}
