import torch
import torch.nn.functional as F
import torch.nn as nn


class _MaskedLoss(torch.nn.Module):
    def forward(self, estimate, output, mask=None):
        feature_mask = mask.expand_as(estimate)
        return self._loss(estimate[feature_mask], output[feature_mask])


class L1Loss(_MaskedLoss):
    def __init__(self):
        super().__init__()
        self._loss = torch.nn.L1Loss()


class L2Loss(_MaskedLoss):
    def __init__(self):
        super().__init__()
        self._loss = torch.nn.MSELoss()


class AELoss(nn.Module):
    def __init__(self, lambda1, lambda2):
        super().__init__()
        self.lambda1 = lambda1
        self.lambda2 = lambda2
        self.l1_loss = nn.L1Loss()
        self.l2_loss = nn.MSELoss()

    def forward(self, outputs, targets):
        l1_loss = self.l1_loss(outputs, targets)
        l2_loss = self.l2_loss(outputs, targets)
        combined_loss = self.lambda1 * l1_loss + self.lambda2 * l2_loss
        return combined_loss


class CustomClipLoss(torch.nn.Module):
    """Modified CLIP contrastive loss with weights for positive and negative samples."""

    def __init__(
        self, linear=None, twin=True, center=False, temp_tau=1.0, negative_weight=10.0
    ):
        super().__init__()
        self.linear = None
        self.center = center
        self.negative_weight = negative_weight
        if linear is not None:
            self.linear_est = torch.nn.LazyLinear(linear)
            if twin:
                self.linear_gt = self.linear_est
            else:
                self.linear_gt = torch.nn.LazyLinear(linear)
        self.temp_tau = nn.Parameter(torch.tensor(temp_tau))

    def get_scores(self, estimates: torch.Tensor, candidates: torch.Tensor):
        """Given estimates that is [B, N] and candidates which is [B', N],
        return a [B, B'] matrix of scores of matching."""
        if self.linear:
            estimates = self.linear_est(estimates)
            candidates = self.linear_gt(candidates)
        if self.center:
            estimates = estimates - estimates.mean(dim=1, keepdim=True)
            candidates = candidates - candidates.mean(dim=1, keepdim=True)

        inv_norms = 1 / (1e-8 + candidates.norm(dim=1, p=2))
        inv_norms_2 = 1 / (1e-8 + estimates.norm(dim=1, p=2))
        # scores = torch.einsum("bn,on,o->bo", estimates, candidates, inv_norms)
        # scores = torch.einsum("bn,bn->b", estimates, candidates)
        scores = torch.einsum(
            "bn,on,b,o -> bo", estimates, candidates, inv_norms_2, inv_norms
        )
        return scores

    def get_probabilities(self, estimates, candidates):
        """Given estimates that is [B, N] and candidates which is [B', N],
        return a [B, B'] matrix of probabilities of matching."""
        scores = self.get_scores(estimates, candidates)
        scores = scores / self.temp_tau
        return F.softmax(scores, dim=1)

    def _weighted_logits(self, logits):
        """Add a hard-negative bias to off-diagonal logits *before* softmax.

        Scaling post-softmax probabilities and feeding them back into
        cross_entropy (the previous implementation) applies softmax twice,
        which is not a valid reweighting. Adding log(negative_weight) to the
        off-diagonal logits is the correct equivalent: after softmax, the
        negative terms' probability mass is scaled by exactly
        `negative_weight` relative to the positive term, so the model is
        pushed harder to separate positives from negatives without breaking
        the softmax normalization.
        """
        n_rows, n_cols = logits.shape
        mask = ~torch.eye(n_rows, n_cols, dtype=torch.bool, device=logits.device)
        bias = torch.zeros_like(logits)
        bias[mask] = torch.log(torch.as_tensor(self.negative_weight, device=logits.device, dtype=logits.dtype))
        return logits + bias

    def forward(self, estimate, candidate):
        """Forward method for ClipLoss.

        Symmetrized like standard CLIP: loss is averaged over the
        estimate->candidate direction and the candidate->estimate direction
        (restricted to the first len(estimate) candidates, the only ones
        with a known positive match), instead of only supervising one
        direction.
        """
        assert estimate.size(0) <= candidate.size(0), (
            "need at least as many targets as estimates"
        )
        n = estimate.size(0)
        logits = self.get_scores(estimate, candidate) / self.temp_tau  # [n, n_candidates]
        weighted_logits = self._weighted_logits(logits)
        target = torch.arange(n, device=estimate.device)

        loss_e2c = F.cross_entropy(weighted_logits, target)
        loss_c2e = F.cross_entropy(weighted_logits[:, :n].t(), target)
        return 0.5 * (loss_e2c + loss_c2e)


def clip_prob(estimates, candidates, temp_tau=1.0):
    inv_norms = 1 / (1e-8 + candidates.norm(dim=1, p=2))
    inv_norms_2 = 1 / (1e-8 + estimates.norm(dim=1, p=2))
    # scores = torch.einsum("bn,on,o->bo", estimates, candidates, inv_norms)
    # scores = torch.einsum("bn,bn->b", estimates, candidates)
    scores = torch.einsum(
        "bn,on,b,o -> bo", estimates, candidates, inv_norms_2, inv_norms
    )
    scores = scores / temp_tau
    return F.softmax(scores, dim=1)


def clip_sim(estimates, candidates):
    inv_norms = 1 / (1e-8 + candidates.norm(dim=1, p=2))
    inv_norms_2 = 1 / (1e-8 + estimates.norm(dim=1, p=2))
    # scores = torch.einsum("bn,on,o->bo", estimates, candidates, inv_norms)
    # scores = torch.einsum("bn,bn->b", estimates, candidates)
    scores = torch.einsum(
        "bn,on,b,o -> bo", estimates, candidates, inv_norms_2, inv_norms
    )
    return scores


def off_diagonal(x):
    n, m = x.shape
    assert n == m
    return x.flatten()[:-1].view(n - 1, n + 1)[:, 1:].flatten()


class Projector(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim, n_hidden_layers, dropout):
        super().__init__()
        self.projector = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(inplace=True),
            nn.Dropout(p=dropout),
            nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x):
        return self.projector(x)
