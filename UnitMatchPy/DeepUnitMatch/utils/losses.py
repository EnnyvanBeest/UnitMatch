import numpy as np
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
    """
    Modified CLIP contrastive loss with weights for positive and negative samples.

    For a batch of B neurons with cosine similarities c_ij between the
    projection of half 1 of neuron i and half 2 of neuron j:
        s_ij  = exp(c_ij / tau) / sum_k exp(c_ik / tau)     (row-wise softmax)
        s'_ij = s_ij * W_ij,  W_ii = 1, W_ij = negative_weight (i != j)
        L     = 1/B sum_i -log( s'_ii / sum_j s'_ij )
    which equals cross-entropy on the logits c_ij / tau + log(W_ij).
    tau is a learned temperature (parameterised as log(tau), so it stays
    positive), initialised at temp_tau.
    """

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
        self.log_temp_tau = nn.Parameter(torch.log(torch.tensor(float(temp_tau))))

    @property
    def temp_tau(self):
        return self.log_temp_tau.exp()

    def _load_from_state_dict(self, state_dict, prefix, *args, **kwargs):
        # checkpoints from before the log parameterisation store "temp_tau" itself
        old_key = prefix + "temp_tau"
        if old_key in state_dict and prefix + "log_temp_tau" not in state_dict:
            state_dict[prefix + "log_temp_tau"] = torch.log(state_dict.pop(old_key))
        super()._load_from_state_dict(state_dict, prefix, *args, **kwargs)

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

    def forward(self, estimate, candidate):
        """
        Weighted contrastive loss (see class docstring).

        Previously the weights multiplied the softmax *probabilities*, which
        were then passed to F.cross_entropy -- that applies a second softmax,
        so the loss was not the one described (-log s'_ii / sum_j s'_ij).
        Adding log(W) to the logits before a single softmax gives exactly that.
        """
        assert estimate.size(0) <= candidate.size(0), (
            "need at least as many targets as estimates"
        )
        logits = self.get_scores(estimate, candidate) / self.temp_tau
        off_diagonal = ~torch.eye(
            logits.size(0), logits.size(1), dtype=torch.bool, device=logits.device
        )
        log_weights = torch.zeros_like(logits)
        log_weights[off_diagonal] = float(np.log(self.negative_weight))
        target = torch.arange(logits.size(0), device=estimate.device)
        return F.cross_entropy(logits + log_weights, target)


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
