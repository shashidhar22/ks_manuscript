"""ksplots CLI: ``ksplots figure 4``, ``ksplots all``, ``ksplots list``."""
from __future__ import annotations

from pathlib import Path

import click

from . import figures
from .config import FIGURES_DIR
from .plotting import apply_default_theme


@click.group()
def main() -> None:
    """Render KS manuscript figures."""
    apply_default_theme()


@main.command("list")
def list_cmd() -> None:
    """List available figures."""
    for name in figures.available():
        click.echo(name)


@main.command("figure")
@click.argument("name")
@click.option(
    "--out",
    "out_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=FIGURES_DIR,
    show_default=True,
    help="Output directory.",
)
def figure_cmd(name: str, out_dir: Path) -> None:
    """Render a single figure (e.g. 4, s2)."""
    path = figures.render(name, out_dir)
    click.echo(f"wrote {path}")


@main.command("all")
@click.option(
    "--out",
    "out_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=FIGURES_DIR,
    show_default=True,
)
def all_cmd(out_dir: Path) -> None:
    """Render every figure."""
    for name in figures.available():
        try:
            path = figures.render(name, out_dir)
            click.echo(f"  {name}: {path}")
        except Exception as exc:  # noqa: BLE001
            click.echo(f"  {name}: FAILED — {exc}", err=True)


if __name__ == "__main__":
    main()
