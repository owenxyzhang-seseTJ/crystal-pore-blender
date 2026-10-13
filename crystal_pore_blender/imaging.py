"""Composite a transparent render onto white, checking visible bounds."""

from pathlib import Path
from PIL import Image


def composite_white(path):
    path = Path(path)
    with Image.open(path) as source:
        image = source.convert("RGBA")
        if source.mode != "RGBA":
            raise ValueError("Expected an RGBA transparent render.")
    alpha = image.getchannel("A")
    bounds = alpha.getbbox()
    if bounds is None or alpha.getextrema()[0] != 0:
        raise ValueError("Empty render or no transparent background.")
    if (
        bounds[0] <= 1
        or bounds[1] <= 1
        or bounds[2] >= image.width - 1
        or bounds[3] >= image.height - 1
    ):
        raise ValueError(
            "Visible render touches the image border; increase camera margin."
        )
    white = Image.new("RGBA", image.size, "white")
    white.alpha_composite(image)
    target = path.with_name(path.stem.replace("_transparent", "") + "_white.png")
    white.convert("RGB").save(target)
    return target


def main():
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("image", type=Path)
    print(composite_white(p.parse_args().image))
