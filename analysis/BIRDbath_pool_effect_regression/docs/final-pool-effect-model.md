# Final pool-effect model

We evaluated pool-associated differences in BIRD effect estimates using the
log2-scale regression

$$\log_2(\theta^{BIRD}_{vp})=\alpha+\gamma\log_2(\theta^{BIRDbath}_{v})+\beta_p+\epsilon_{vp},\qquad\sum_p\beta_p=0.$$

The common slope $\gamma$ describes the overall relationship between
BIRDbath and BIRD effect magnitudes. The intercept $\alpha$ represents the
expected BIRD value for the average pool when BIRDbath theta equals 1;
retaining it avoids forcing the average-pool relationship through the origin.
The pool coefficient $\beta_p$ represents pool $p$'s additive deviation
from the average pool on the log2 scale, while $\epsilon_{vp}$ represents
variation not explained by the shared BIRDbath relationship or pool identity.
Pool effects were constrained to sum to zero so that they are interpreted
relative to the average pool rather than an arbitrarily selected reference
pool. Because the same variant may occur in multiple pools, confidence
intervals were calculated using variant-clustered standard errors. This model
provides a simple, interpretable screen for systematic pool-associated shifts
conditional on BIRDbath theta; it is not intended as a comprehensive model of
the relationship between the two assays.

