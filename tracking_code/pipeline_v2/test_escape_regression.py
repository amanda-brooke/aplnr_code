"""Regression test: pipeline_v2/ps_escape reproduces the validated Fig 3 numbers."""
import sys, os; sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, pandas as pd, ps_escape as pe
from scipy.stats import mannwhitneyu
d = pd.read_pickle('data/shaplnr_mps_radial_edges.pkl')
R = d['diameter_um'].iloc[0] / 2
d['dist_from_edge'] = R - d['radial_distance']
edges, ctr = pe.make_bins(0, 300, 30)
res, info = pe.run_escape(d, 'dist_from_edge', edges, ref_window=(42, 48), k=1.5, group_col='Group')
L = res[(res.t_bin >= 45) & (res.t_bin <= 51)].groupby(['File_ID', 'group']).excess.mean().reset_index()
c = L[L.group == 'Control'].excess; k = L[L.group == 'Treatment'].excess
p = mannwhitneyu(c, k).pvalue
with open('pipeline_v2/_regression_out.txt', 'w') as f:
    f.write(f"boundary={info['boundary']:.1f} um (expect ~110)\n")
    f.write(f"true escape late: Control={c.median():.1f}% shAPLNR={k.median():.1f}% p={p:.4g}\n")
    f.write(f"EXPECTED: ~110, Control~17.5, shAPLNR~5.9, p<0.0001\n")
    ok = abs(info['boundary']-110)<8 and abs(c.median()-17.5)<4 and p<0.001
    f.write("PASS\n" if ok else "CHECK\n")
print("done")
