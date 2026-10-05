"""End-to-end test of train_paper_models.py (comparison jobs) + compare_training_variants.py.

Run:  python PaperAnalyses/tests/test_train_paper_models.py      (~10-20 min on CPU)

Checks: (1) comparison job list (3 replicates x 3 variants, one AE per replicate);
(2) training on small fake snippet data (1 AE epoch, 2 fine-tuning epochs): AE once, variants
A/B/C fine-tuned from it and published, only C uses channel positions, AE re-fetched from the
shared models folder as on another machine; (3) evaluation of A/B/C through the unified
pipeline on one real location of the previous run (AL032/19011111882/1, READ only);
(4) the comparison report.

All pipeline_config paths are redirected to a local temp folder (printed at the start).
Throwaway experiments in DeepUnitMatch/ModelExp (pytest_drv_ae, cmp_m3_1_*) are removed at
the end. Ends with "ALL CHECKS PASSED".
"""
import tempfile
import json
import os
import shutil
import sys

import numpy as np

PA = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../PaperAnalyses
sys.path.insert(0, PA)


def main():
    import pipeline_config as cfg
    ROOT = os.path.join(tempfile.gettempdir(), "dum_test_train_paper_models")
    print("test output folder:", ROOT)
    shutil.rmtree(ROOT, ignore_errors=True)
    cfg.TRAINING_CACHE = os.path.join(ROOT, "cache")
    cfg.TRAINING_STATE_ROOT = os.path.join(ROOT, "state")
    cfg.MODELS_ROOT = os.path.join(ROOT, "models")
    cfg.ANALYSIS_OUTPUT = os.path.join(ROOT, "out")
    cfg.REPORTS_DIR = os.path.join(ROOT, "reports")
    cfg.LOG_DIR = os.path.join(ROOT, "logs")

    import train_paper_models as tpm
    import run_deepunitmatch_batch_onMerged as onm
    from DeepUnitMatch.utils import param_fun

    tpm.AE_EPOCHS, tpm.FT_EPOCHS, tpm.AE_BATCHSIZE, tpm.FT_BATCHSIZE = 1, 2, 8, 8
    AE = "pytest_drv_ae"
    names = [f"cmp_m3_1_{v}" for v in tpm.VARIANTS]

    def cleanup():
        for d in [os.path.join(tpm.MODELEXP, "AE_experiments", AE), os.path.join(tpm.MODELEXP, "experiments", AE)] + \
                 [os.path.join(tpm.MODELEXP, "experiments", n) for n in names]:
            shutil.rmtree(d, ignore_errors=True)
    cleanup()

    # 1. job list from the real manifest: 3 reps x 3 variants, AE shared per rep
    manifest = tpm.generate_manifest()
    jobs = tpm.comparison_jobs(manifest)
    assert len(jobs) == 9 and len({j["ae"] for j in jobs}) == 3
    assert all(j["train_mice"] == manifest[f"m3_{j['name'].split('_')[2]}"] for j in jobs)
    print("1. comparison jobs:", [j["name"] for j in jobs[:3]], "... AE per rep:", sorted({j["ae"] for j in jobs}))

    # 2. fake training cache (3 'mice', 1 location each, 3 sessions) written by the real get_snippets
    rng = np.random.default_rng(0)
    pos = np.column_stack([np.tile([0.0, 32.0], 192), np.repeat(np.arange(192) * 15.0, 2)])
    groups = {}
    for mouse in ("MA", "MB", "MC"):
        group = f"{mouse}/p/1"
        n = 30
        wf = rng.normal(0, 0.1, (n, 82, 384, 2))
        for u in range(n):
            pk = 20 + 10 * u
            wf[u, 35:45, pk - 4:pk + 4, :] += -5.0
        sid = np.repeat([0, 1, 2], n // 3)
        out = tpm.cache_dir(group)
        param_fun.get_snippets(wf, [pos] * 3, sid, save_path=out, unit_ids=np.arange(n),
                               param=param_fun.get_default_param({"nTime": 82, "nChannels": 384}))
        open(os.path.join(out, ".done"), "w").write("ok")
        groups[group] = None
    real_group = "AL032/19011111882/1"
    groups[real_group] = r"\\znas.cortexlab.net\Lab\Share\UNITMATCHTABLES_ENNY_CELIAN_JULIE\DeepUM_NatMeth2026V2_merged\merged_data_v2\AL032\19011111882\1\DeepUnitMatch"

    test_jobs = [tpm.make_job(n, AE, ["MA", "MB", "MC"], v, evaluate=["AL032"]) for n, v in zip(names, tpm.VARIANTS)]
    try:
        # 3. train: AE once, then A (legacy), B, C; published to MODELS_ROOT
        for job in test_jobs:
            tpm.train_job(job, groups)
        for job in test_jobs:
            assert tpm.read_status(job["name"]).get("trained"), job["name"]
            assert os.path.isfile(os.path.join(tpm.published_dir(job["name"]), "model.pt"))
        assert tpm.read_status(AE).get("trained")
        m = {n: tpm.test.load_trained_model(read_path=os.path.join(tpm.published_dir(n), "model.pt")) for n in names}
        assert m["cmp_m3_1_A_original"].uses_channel_pos is False and m["cmp_m3_1_B_fixed"].uses_channel_pos is False
        assert m["cmp_m3_1_C_v2"].uses_channel_pos is True
        print("3. training OK: one shared AE, A/B/C fine-tuned and published; only C uses channel positions")

        # AE fetched from the share when not local (as on another machine)
        shutil.rmtree(os.path.join(tpm.MODELEXP, "AE_experiments", AE))
        tpm.ensure_local_ae(AE)
        assert os.listdir(os.path.join(tpm.MODELEXP, "AE_experiments", AE, "ckpt")) == ["ckpt_epoch_0"]
        print("   AE re-fetched from the shared models folder")

        # 4. evaluate the three models on one real held-out location (read-only input)
        tpm.evaluate(test_jobs, groups)
        for n in names:
            assert tpm.eval_done(real_group, n), n
            assert os.path.isfile(os.path.join(tpm.eval_dir(real_group, n) + "_AssignUniqueID", "AUC_summary.json"))
        print("4. evaluation OK on", real_group)
        tpm.print_status(test_jobs, groups)

        # 5. summary report
        json.dump({"m3_1": ["MA", "MB", "MC"], "m3_2": [], "m3_3": []}, open(os.path.join(cfg.TRAINING_STATE_ROOT, "manifest.json"), "w"))
        import compare_training_variants as ctv
        ctv.main()
        assert os.listdir(cfg.REPORTS_DIR)
        print("5. summary report OK")
        print("ALL CHECKS PASSED")
    finally:
        cleanup()


if __name__ == '__main__':
    main()
