#!/usr/bin/env python
"""
Write an IGV (desktop) session XML that loads the lincorfz outputs together with the
supporting tracks, using absolute paths (the session works from any location on the cluster).

Track table (--tracks, TSV) columns:
  panel     'data' or 'feature'
  group     overlay name (tracks with the same group are merged into one overlay track) or '-'
  name      display name
  path      absolute path (bigWig, GTF, BED, BAM)
  color     R,G,B
  range     '-' (autoscale) or 'min,max' (fixed data range, e.g. 1,10 for MRL)
  height    track height in pixels
"""
import argparse
import uuid
from xml.sax.saxutils import quoteattr

import pandas as pd


def data_track(r, autoscale_group=None, inner=False):
    fixed = r["range"] not in ("-", "", None) and isinstance(r["range"], str)
    attrs = {
        "altColor": r["color"], "attributeKey": r["name"], "autoScale": "false" if fixed or autoscale_group else "true",
        "color": r["color"], "fontSize": "10" if inner else "12", "id": r["path"], "name": r["name"],
        "renderer": "BAR_CHART", "visible": "true", "windowFunction": "mean"}
    if not inner:
        attrs["clazz"] = "org.broad.igv.track.DataSourceTrack"
        attrs["height"] = str(r["height"])
    if autoscale_group and not fixed:
        attrs["autoscaleGroup"] = str(autoscale_group)
    a = " ".join(f"{k}={quoteattr(str(v))}" for k, v in attrs.items())
    if fixed:
        lo, hi = r["range"].split(",")
        return (f"<Track {a}>\n<DataRange baseline={quoteattr(lo)} drawBaseline=\"true\" flipAxis=\"false\" "
                f"maximum={quoteattr(hi)} minimum={quoteattr(lo)} type=\"LINEAR\"/>\n</Track>")
    return f"<Track {a}/>"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--tracks", required=True)
    ap.add_argument("--locus", default="chr20:49278178-49299600")
    ap.add_argument("--out", default="lincorfz.igv_session.xml")
    args = ap.parse_args()

    t = pd.read_csv(args.tracks, sep="\t", dtype=str).fillna("-")
    res = "\n".join(f'        <Resource path={quoteattr(p)} type={quoteattr(p.rsplit(".", 1)[-1])}/>'
                    for p in t["path"])
    data_xml, feat_xml = [], []
    # reference sequence with three-frame translation (+ strand; click the strand arrow in IGV to switch)
    feat_xml.append('<Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" '
                    'fontSize="10" id="Reference sequence" name="Reference sequence" '
                    'sequenceTranslationStrandValue="POSITIVE" shouldShowTranslation="true" visible="true"/>')
    group_id = 0
    for g, sub in t.loc[t["panel"] == "data"].groupby("group", sort=False):
        if g == "-":
            data_xml += [data_track(r) for _, r in sub.iterrows()]
            continue
        group_id += 1
        r0 = sub.iloc[0]
        inner = "\n".join(data_track(r, autoscale_group=group_id, inner=True) for _, r in sub.iterrows())
        data_xml.append(f'<Track attributeKey="Overlay" autoScale="false" autoscaleGroup="{group_id}" '
                        f'clazz="org.broad.igv.track.MergedTracks" fontSize="12" height="{r0["height"]}" '
                        f'id="{uuid.uuid4()}" name={quoteattr(g)} visible="true">\n{inner}\n</Track>')
    for _, r in t.loc[t["panel"] == "feature"].iterrows():
        feat_xml.append(f'<Track attributeKey={quoteattr(r["name"])} clazz="org.broad.igv.track.FeatureTrack" '
                        f'displayMode="EXPANDED" fontSize="12" groupByStrand="false" id={quoteattr(r["path"])} '
                        f'name={quoteattr(r["name"])} visible="true"/>')
    xml = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome={quoteattr(args.genome)} locus={quoteattr(args.locus)} version="8">
    <Resources>
{res}
    </Resources>
    <Panel name="DataPanel">
{chr(10).join(data_xml)}
    </Panel>
    <Panel name="FeaturePanel">
{chr(10).join(feat_xml)}
    </Panel>
    <PanelLayout dividerFractions="0.65"/>
</Session>
"""
    with open(args.out, "w") as fh:
        fh.write(xml)


if __name__ == "__main__":
    main()
