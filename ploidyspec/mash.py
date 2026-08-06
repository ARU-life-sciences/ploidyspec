"""
Mash-style (Ondov et al. 2016) k-mer-size correction for canonical-k-mer Jaccard
distance, plus an explicit noise-floor check for when a given k has run out of
resolving power.

Model: under a simple per-site substitution model with divergence p, a random
k-mer survives identical between two sequences with probability q = (1-p)^k,
giving J ~= q / (2 - q). Inverting:

    distance(J, k) = -1/k * ln(2J / (1+J))     -> 0 as J -> 1, -> inf as J -> 0

That correction is only meaningful while the observed shared-kmer count is well
above what two *unrelated* sequences would share purely by chance, given the
finite canonical k-mer space (4**k / 2). Below that floor, the "distance" is
noise, not a measurement -- see resolution_z / summarize_pair_across_k.
"""

import math

RESOLUTION_Z_THRESHOLD = 5.0


def jaccard_to_distance(jaccard, k):
    """Mash distance from canonical-kmer Jaccard. None if undefined (jaccard <= 0)."""
    if jaccard is None or jaccard <= 0:
        return None
    if jaccard >= 1:
        return 0.0
    return -1.0 / k * math.log(2.0 * jaccard / (1.0 + jaccard))


def containment_to_p(containment, k):
    """Secondary divergence estimate from the (asymmetric) containment index.
    More robust than the Jaccard-based estimate when the two sequences have very
    different k-mer counts, since it doesn't assume comparable set sizes."""
    if containment is None or containment <= 0:
        return None
    if containment >= 1:
        return 0.0
    return 1.0 - containment ** (1.0 / k)


def expected_chance_shared(n1, n2, k, canonical=True):
    """Expected number of k-mers shared between two *unrelated* sequences with
    n1, n2 distinct canonical k-mers, given the finite k-mer space size."""
    space = 4.0**k / (2.0 if canonical else 1.0)
    return (n1 * n2) / space


def resolution_z(shared, n1, n2, k, canonical=True):
    """Poisson-ish z-score of the observed shared count above the chance floor.
    None if the expected chance-shared count is 0 (degenerate/empty inputs)."""
    expected = expected_chance_shared(n1, n2, k, canonical=canonical)
    if expected <= 0:
        return None
    return (shared - expected) / math.sqrt(expected)


def is_resolution_limited(shared, n1, n2, k, canonical=True, threshold=RESOLUTION_Z_THRESHOLD):
    z = resolution_z(shared, n1, n2, k, canonical=canonical)
    if z is None:
        return True
    return z < threshold


def summarize_pair_across_k(k_values, per_k_stats, canonical=True, threshold=RESOLUTION_Z_THRESHOLD):
    """
    per_k_stats: {k: (shared, n1, n2)} for every k in k_values.

    Walks k from largest to smallest and reports the distance from the largest k
    that is NOT resolution-limited (best specificity that still carries real
    signal). If every swept k is resolution-limited, the pair is reported
    unresolved rather than emitting a number that's indistinguishable from noise.

    Returns:
      {
        chosen_k: int or None,
        distance: float or None,        # Mash distance at chosen_k, clipped to [0,1]
        raw_distance: float or None,     # same, uncapped
        containment_p: float or None,    # containment-based estimate at chosen_k, QC only
        resolution_limited: bool,        # True if no swept k cleared the noise floor
        z_at_chosen_k: float or None,
        k_consistency_spread: float or None,  # max-min corrected distance among valid k's
        per_k: {k: {...}},
      }
    """
    per_k_out = {}
    valid_distances = {}

    for k in k_values:
        shared, n1, n2 = per_k_stats[k]
        union = n1 + n2 - shared
        jaccard = shared / union if union > 0 else 0.0
        min_n = min(n1, n2)
        containment = shared / min_n if min_n > 0 else 0.0
        z = resolution_z(shared, n1, n2, k, canonical=canonical)
        limited = z is None or z < threshold
        dist = jaccard_to_distance(jaccard, k)
        cont_p = containment_to_p(containment, k)
        per_k_out[k] = dict(
            shared=shared,
            n1=n1,
            n2=n2,
            jaccard=jaccard,
            containment=containment,
            distance=dist,
            containment_p=cont_p,
            z=z,
            resolution_limited=limited,
        )
        if not limited and dist is not None:
            valid_distances[k] = dist

    if not valid_distances:
        # nothing cleared the noise floor at any swept k -- unresolved, not a fabricated number
        largest_k = max(k_values)
        return dict(
            chosen_k=None,
            distance=None,
            raw_distance=None,
            containment_p=None,
            resolution_limited=True,
            z_at_chosen_k=per_k_out[largest_k]["z"],
            k_consistency_spread=None,
            per_k=per_k_out,
        )

    chosen_k = max(valid_distances)  # largest non-resolution-limited k -> best specificity
    raw_distance = valid_distances[chosen_k]
    spread = (
        max(valid_distances.values()) - min(valid_distances.values())
        if len(valid_distances) > 1
        else 0.0
    )

    return dict(
        chosen_k=chosen_k,
        distance=min(max(raw_distance, 0.0), 1.0),
        raw_distance=raw_distance,
        containment_p=per_k_out[chosen_k]["containment_p"],
        resolution_limited=False,
        z_at_chosen_k=per_k_out[chosen_k]["z"],
        k_consistency_spread=spread,
        per_k=per_k_out,
    )
