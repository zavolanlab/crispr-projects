import numpy as np
import pandas as pd
import pyBigWig

from lincdisease import encode


def _bigwig(path, value_at):
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr1", 1000)])
    starts = list(range(0, 1000, 100))
    bw.addEntries(["chr1"] * 10, starts, ends=[s + 100 for s in starts], values=[float(value_at(s)) for s in starts])
    bw.close()


def test_condition_signal_uses_feature_strand(tmp_path):
    # plus strand carries 1.0 on 0-500, minus strand 3.0 on 500-1000
    _bigwig(tmp_path / "P.bigWig", lambda s: 1.0 if s < 500 else 0.0)
    _bigwig(tmp_path / "M.bigWig", lambda s: 3.0 if s >= 500 else 0.0)
    regions = pd.DataFrame({"feature_id": ["a", "b", "b"], "chrom": "chr1", "start0": [0, 500, 900],
                            "end1": [200, 600, 1000]})
    sig = encode.condition_signal(regions, pd.Series({"a": "+", "b": "-"}), [("P", "M")], tmp_path)
    assert sig["a"] == 1.0 and sig["b"] == 3.0


def test_log_ratio_masks_unexpressed():
    r = encode.log_ratio(pd.Series([1.0, 0.0, 0.001]), pd.Series([0.25, 0.0, 0.002]))
    assert np.isclose(r[0], np.log2(1.01 / 0.26)) and r[1:].isna().all()
