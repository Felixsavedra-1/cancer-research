from __future__ import annotations

_HIGHLIGHT_COLORS = ["magenta", "orange", "lime", "cyan", "yellow", "red"]

_FLEX_STOPS = [(60, 90, 230), (230, 230, 230), (220, 60, 60)]
_PRIORITY_STOPS = [
    (40, 40, 46),
    (110, 92, 52),
    (198, 161, 91),
    (178, 58, 58),
]


def _lerp_color(value: float, stops: list[tuple[int, int, int]]) -> str:
    value = max(0.0, min(1.0, value))
    if len(stops) == 1:
        r, g, b = stops[0]
        return f"0x{r:02x}{g:02x}{b:02x}"
    span = value * (len(stops) - 1)
    i = min(int(span), len(stops) - 2)
    t = span - i
    r = round(stops[i][0] + (stops[i + 1][0] - stops[i][0]) * t)
    g = round(stops[i][1] + (stops[i + 1][1] - stops[i][1]) * t)
    b = round(stops[i][2] + (stops[i + 1][2] - stops[i][2]) * t)
    return f"0x{r:02x}{g:02x}{b:02x}"


def _color_by_values(view, values: dict[int, float], stops, bins: int = 12) -> None:
    buckets: dict[int, list[int]] = {}
    for resi, value in values.items():
        idx = min(bins - 1, max(0, int(value * bins)))
        buckets.setdefault(idx, []).append(resi)
    for idx, residues in buckets.items():
        color = _lerp_color((idx + 0.5) / bins, stops)
        view.addStyle(
            {"resi": [str(r) for r in residues]}, {"cartoon": {"color": color}}
        )


def render_structure(
    pdb_text: str,
    highlights: list[dict] | None = None,
    color_mode: str = "plddt",
    flexibility: dict[int, float] | None = None,
    priority: dict[int, float] | None = None,
    pocket_residues: list[int] | None = None,
    width: int = 900,
    height: int = 600,
) -> str:
    import py3Dmol

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")

    if color_mode == "flexibility" and flexibility:
        view.setStyle({"cartoon": {"color": "0x404046"}})
        _color_by_values(view, flexibility, _FLEX_STOPS)
    elif color_mode == "priority" and priority:
        view.setStyle({"cartoon": {"color": "0x2a2a2e"}})
        _color_by_values(view, priority, _PRIORITY_STOPS)
    elif color_mode == "plddt":
        view.setStyle(
            {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}}
        )
    else:
        view.setStyle({"cartoon": {"color": "spectrum"}})

    if pocket_residues:
        view.addSurface(
            "VDW",
            {"opacity": 0.55, "color": "0xC6A15B"},
            {"resi": [str(r) for r in pocket_residues]},
        )

    for index, hit in enumerate(highlights or []):
        position = str(hit["position"])
        color = hit.get("color") or _HIGHLIGHT_COLORS[index % len(_HIGHLIGHT_COLORS)]
        selection = {"resi": position}
        view.addStyle(selection, {"stick": {"color": color, "radius": 0.3}})
        view.addStyle(selection, {"sphere": {"color": color, "opacity": 0.6}})
        view.addResLabels(
            selection,
            {
                "fontSize": 11,
                "fontColor": "black",
                "backgroundColor": color,
                "backgroundOpacity": 0.75,
            },
        )

    view.zoomTo()
    view.setBackgroundColor("0x0A0A0B")
    return view._make_html()
