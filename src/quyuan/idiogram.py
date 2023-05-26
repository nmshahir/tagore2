import pandas as pd
import svg


def idiogram(karyotype: pd.DataFrame, total_width=170) -> svg.SVG:
    has_centromere = len(karyotype.columns) == 5
    mpx = 3.543307  # "1mm" would be "3.543307px"
    layout_factor = 2.6
    width = total_width / (layout_factor * len(karyotype)) * mpx  # width(mm)
    left_margin = 20 * mpx
    top_margin = 25 * mpx

    karyotype["x"] = [left_margin + x * layout_factor * width for x in range(len(karyotype))]
    maxchrlen = 150  # the length of the longest chromosome was set to be 150mm
    max_chrom_len = karyotype["End"].max()
    if has_centromere:
        karyotype["y_CE_start"] = maxchrlen * (1 - (karyotype["End"] - karyotype["CE_start"]) / max_chrom_len) * mpx + top_margin
        karyotype["y_CE_end"] = maxchrlen * (1 - (karyotype["End"] - karyotype["CE_end"]) / max_chrom_len) * mpx + top_margin
    karyotype["y_start"] = maxchrlen * (1 - karyotype["End"] / max_chrom_len) * mpx + top_margin
    karyotype["y_end"] = maxchrlen * mpx + top_margin
    karyotype["path"] = karyotype.apply(
        lambda row: svg.Path(
            fill="none",
            stroke="grey",
            stroke_width=1,
            d=[
                svg.MoveTo(row["x"], row["y_start"] + width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, True, True, row["x"] + width, row["y_start"] + width / 2
                ),
                svg.LineTo(row["x"] + width, row["y_CE_start"]),
                svg.LineTo(row["x"], row["y_CE_end"]),
                svg.LineTo(row["x"], row["y_end"] - width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, False, False, row["x"] + width, row["y_end"] - width / 2
                ),
                svg.LineTo(row["x"] + width, row["y_CE_end"]),
                svg.LineTo(row["x"], row["y_CE_start"]),
                svg.ClosePath(),
            ]
            if has_centromere
            else [
                svg.MoveTo(row["x"], row["y_start"] + width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, True, True, row["x"] + width, row["y_start"] + width / 2
                ),
                svg.LineTo(row["x"] + width, row["y_end"] - width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, True, True, row["x"], row["y_end"] - width / 2
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
                svg.MoveTo(row["x"], row["y_start"] + width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, True, True, row["x"] + width, row["y_start"] + width / 2
                ),
                svg.LineTo(row["x"] + width, row["y_start"]),
                svg.LineTo(row["x"], row["y_start"]),
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
                svg.MoveTo(row["x"], row["y_end"] - width / 2),
                svg.Arc(
                    width / 2, width / 2, 0, False, False, row["x"] + width, row["y_end"] - width / 2
                ),
                svg.LineTo(row["x"] + width, row["y_end"]),
                svg.LineTo(row["x"], row["y_end"]),
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
                    svg.MoveTo(row["x"], row["y_CE_start"]),
                    svg.LineTo(row["x"] + width, row["y_CE_end"]),
                    svg.LineTo(row["x"] + width, row["y_CE_start"]),
                    svg.LineTo(row["x"], row["y_CE_end"]),
                    svg.ClosePath(),
                ],
            ),
            axis=1,
        )
    karyotype["text"] = karyotype.apply(
        lambda row: svg.Text(
            x=(row["x"] + row["x"] + width) / 2 - len(str(row["Chr"])) * 2.2,
            y=maxchrlen * mpx + top_margin + 15,
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
