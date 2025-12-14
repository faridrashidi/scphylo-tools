import click

import scphylo as scp


@click.command(short_help="Calculate scores between two trees.")
@click.argument(
    "ground_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
@click.argument(
    "inferred_file",
    required=True,
    type=click.Path(
        exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
    ),
)
def score(ground_file, inferred_file):
    """Scores of Anscenstor-descendant, Different-lineage, MLTD and etc.

    scphylo utils score ground.CFMatrix inferred.CFMatrix
    """
    scp.settings.verbosity = "info"

    df_g = scp.io.read(ground_file)
    df_s = scp.io.read(inferred_file)

    gs = scp.tl.gs(df_g, df_s)
    scp.logg.info(f"GS={gs:0.4f}")
    ad = scp.tl.ad(df_g, df_s)
    scp.logg.info(f"AD={ad:0.4f}")
    dl = scp.tl.dl(df_g, df_s)
    scp.logg.info(f"DL={dl:0.4f}")
    mltd = scp.tl.mltd(df_g, df_s)["normalized_similarity"]
    scp.logg.info(f"MLTSM={mltd:0.4f}")
    tpted = scp.tl.tpted(df_g, df_s)
    scp.logg.info(f"TPTED={tpted:0.4f}")
    rf = scp.tl.rf(df_g, df_s)
    scp.logg.info(f"RF={rf:0.4f}")
    caset = scp.tl.caset(df_g, df_s)
    scp.logg.info(f"CASet={caset:0.4f}")
    disc = scp.tl.disc(df_g, df_s)
    scp.logg.info(f"DISC={disc:0.4f}")
    mp3 = scp.tl.mp3(df_g, df_s)
    scp.logg.info(f"MP3={mp3:0.4f}")

    return None
