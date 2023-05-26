import pandas as pd
import svg


def idiogram(karyotype: pd.DataFrame, width=170) -> svg.SVG:
    r"""
    9--10
    |/\|
    1  2
    |  |
    8  3
     \/
     /\
    4  7
    |  |
    5  6
    |\/|
    12-11
    """
    has_centromere = len(karyotype.columns) == 5
    mpx = 3.543307  # "1mm" would be "3.543307px"
    layout_factor = 2.6
    chr_width = width / (layout_factor * len(karyotype)) * mpx  # width(mm)
    left_margin = 20 * mpx

    karyotype["x9"] = karyotype["x1"] = karyotype["x8"] = karyotype["x4"] = karyotype[
        "x5"
    ] = karyotype["x12"] = [
        left_margin + x * layout_factor * chr_width for x in range(len(karyotype))
    ]
    maxchrlen = 150  # the length of the longest chromosome was set to be 150mm
    max_chrom_len = karyotype["End"].max()
    karyotype["y1"] = karyotype["y2"] = (
        25 + maxchrlen * (1 - karyotype["End"] / max_chrom_len)
    ) * mpx + chr_width / 2
    karyotype["x10"] = karyotype["x2"] = karyotype["x3"] = karyotype["x7"] = karyotype[
        "x6"
    ] = karyotype["x11"] = (karyotype["x1"] + chr_width)
    if has_centromere:
        karyotype["y3"] = karyotype["y8"] = (
            25
            + maxchrlen
            * (1 - (karyotype["End"] - karyotype["CE_start"]) / max_chrom_len)
        ) * mpx
        karyotype["y4"] = karyotype["y7"] = (
            25
            + maxchrlen * (1 - (karyotype["End"] - karyotype["CE_end"]) / max_chrom_len)
        ) * mpx
    karyotype["y5"] = karyotype["y6"] = (25 + maxchrlen) * mpx - chr_width / 2
    karyotype["y9"] = karyotype["y10"] = (
        25 + maxchrlen * (1 - karyotype["End"] / max_chrom_len)
    ) * mpx
    karyotype["y11"] = karyotype["y12"] = (25 + maxchrlen) * mpx
    karyotype["path"] = karyotype.apply(
        lambda row: svg.Path(
            fill="none",
            stroke="grey",
            stroke_width=1,
            d=[
                svg.MoveTo(row["x1"], row["y1"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, True, True, row["x2"], row["y2"]
                ),
                svg.LineTo(row["x3"], row["y3"]),
                svg.LineTo(row["x4"], row["y4"]),
                svg.LineTo(row["x5"], row["y5"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, False, False, row["x6"], row["y6"]
                ),
                svg.LineTo(row["x7"], row["y7"]),
                svg.LineTo(row["x8"], row["y8"]),
                svg.ClosePath(),
            ]
            if has_centromere
            else [
                svg.MoveTo(row["x1"], row["y1"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, True, True, row["x2"], row["y2"]
                ),
                svg.LineTo(row["x6"], row["y6"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, True, True, row["x5"], row["y5"]
                ),
                svg.ClosePath(),
            ],
        ),
        axis=1,
    )
    karyotype["hat"] = karyotype.apply(
        lambda row: svg.Path(
            fill="white",
            stroke="white",
            stroke_width=0.75,
            d=[
                svg.MoveTo(row["x1"], row["y1"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, True, True, row["x2"], row["y2"]
                ),
                svg.LineTo(row["x10"], row["y10"]),
                svg.LineTo(row["x9"], row["y9"]),
                svg.ClosePath(),
            ],
        ),
        axis=1,
    )
    karyotype["shoe"] = karyotype.apply(
        lambda row: svg.Path(
            fill="white",
            stroke="white",
            stroke_width=0.75,
            d=[
                svg.MoveTo(row["x5"], row["y5"]),
                svg.Arc(
                    chr_width / 2, chr_width / 2, 0, False, False, row["x6"], row["y6"]
                ),
                svg.LineTo(row["x11"], row["y11"]),
                svg.LineTo(row["x12"], row["y12"]),
                svg.ClosePath(),
            ],
        ),
        axis=1,
    )
    if has_centromere:
        karyotype["bow"] = karyotype.apply(
            lambda row: svg.Path(
                fill="white",
                stroke="white",
                stroke_width=0.75,
                d=[
                    svg.MoveTo(row["x8"], row["y8"]),
                    svg.LineTo(row["x7"], row["y7"]),
                    svg.LineTo(row["x3"], row["y3"]),
                    svg.LineTo(row["x4"], row["y4"]),
                    svg.ClosePath(),
                ],
            ),
            axis=1,
        )
    karyotype["text"] = karyotype.apply(
        lambda row: svg.Text(
            x=(row["x1"] + row["x2"]) / 2 - len(str(row["Chr"])) * 2.2,
            y=(maxchrlen + 25) * mpx + 15,
            font_size=9,
            font_family="Arial",
            fill="black",
            text=row["Chr"],
        ),
        axis=1,
    )

    return svg.SVG(
        width=744.094488189,
        height=1052.36220472,
        elements=[
            *karyotype["hat"],
            *karyotype["shoe"],
            *(karyotype["bow"] if has_centromere else []),
            *karyotype["path"],
            *karyotype["text"],
        ],
    )
