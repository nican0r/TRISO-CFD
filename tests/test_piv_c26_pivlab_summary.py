"""src.post.piv_c26_pivlab_summary: collecting the MATLAB score.json files."""

import json

import numpy as np

from src.post import piv_c26_pivlab_summary as S


def _score(case, flow, variant, clahe):
    return {
        "run": f"{case}_x", "case": case, "flow_slpm": flow, "variant": variant,
        "window_t0_s": 2.4, "window_t1_s": 3.0, "n_pairs_spout": 300, "pivlab_version": "3.14",
        "pivlab_settings": {"pre": {"clahe": clahe, "binarize": 0}},
        "map_mean_m_s": 0.16, "map_max_m_s": 0.4, "spout_rejected_frac": 0.1,
        "metrics_vs_exp": {"CC": 0.85, "SSI": 0.85, "MSE_m2_s2": 0.0028},
        "metrics_c26_model_vs_exp": {"CC": 0.87, "SSI": 0.74, "MSE_m2_s2": 0.0022},
        "metrics_vs_c26_model": {"CC": 0.9, "SSI": 0.8, "MSE_m2_s2": 0.001},
        "wall_piv_m_s": 1.5e-2, "wall_piv_sd_m_s": 4e-2, "wall_vectors_in_box": 9,
    }


def test_collect_rows_and_targets(tmp_path):
    run = tmp_path / "c02_x"
    for v, clahe in (("pivlab_default", 1), ("pivlab_noclahe", 0)):
        d = run / v
        d.mkdir(parents=True)
        (d / "score.json").write_text(json.dumps(_score("c02", 12.0, v, clahe)))
        np.savetxt(d / "map_fig4grid.csv", np.full((3, 4), 0.1), delimiter=",")
    (run / "reference_native.json").write_text(json.dumps(
        {"wall_piv_m_s": 1.46e-2, "wall_visible_parcels_m_s": 1.53e-2,
         "metrics_vs_exp": {"CC": 0.91, "SSI": 0.88, "MSE_m2_s2": 0.0015}}))
    rows, maps = S.collect(tmp_path, ["c02_x"])
    assert [r["variant"] for r in rows] == ["pivlab_default", "pivlab_noclahe"]
    r = rows[0]
    # C26 Table 3 for c02 (params/sweep_C26.yaml) and Fig. 6d wall targets (params/piv_c26.yaml)
    assert abs(r["CC_table3"] - 0.817) < 1e-9 and abs(r["MSE_m2_s2_table3"] - 0.0042) < 1e-9
    assert abs(r["wall_exp_m_s"] - 1.83e-3) < 1e-12 and abs(r["wall_c26model_m_s"] - 5.05e-3) < 1e-12
    assert abs(r["wall_over_exp"] - 1.5e-2 / 1.83e-3) < 1e-9
    assert r["repo_native_wall_m_s"] == 1.46e-2 and r["clahe"] == 1 and rows[1]["clahe"] == 0
    assert r["wall_parcels_m_s"] == 1.53e-2 and abs(r["wall_parcels_over_exp"] - 1.53e-2 / 1.83e-3) < 1e-9
    assert maps[("c02_x", "pivlab_default")].shape == (3, 4)


def test_write_csv_roundtrip(tmp_path):
    run = tmp_path / "c01_x"
    (run / "pivlab_default").mkdir(parents=True)
    (run / "pivlab_default" / "score.json").write_text(json.dumps(_score("c01", 10.0, "pivlab_default", 1)))
    rows, _ = S.collect(tmp_path, ["c01_x"])
    out = tmp_path / "t.csv"
    S.write_csv(rows, out)
    lines = out.read_text().splitlines()
    assert len(lines) == 2 and lines[0].startswith("run,case,flow_slpm,variant")
    assert "0.6817" in lines[1]            # Table 3 CC for c01
