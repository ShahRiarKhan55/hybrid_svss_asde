# Stage 3 similarity variant – defined before any Stage 1/2 result exists (2026-10-21)

Variant "SVSS-cos": identical to the current SVSS (alpha = beta = 0.5, features [mu, sigma^2], joint z-score over candidates and training sites, refit per 1% batch, 3% anchors, same random 1% seed), except the similarity term:

- current: cosine between the mean-centred feature vector of the candidate and the mean-centred training centroid (binary in practice, see Stage 0);
- variant: cosine between the feature vector of the candidate and the training centroid, without mean-centring, `sim = (z_i . c) / (|z_i| |c| + 1e-9)`, then min-max over the candidates as before.

Nothing else changes. It is evaluated in Stage 3 on TUNE only, next to the other Stage 3 candidates, and is adopted only under the Stage 3 rule (better than the current hybrid by more than the seed SD on TUNE).
