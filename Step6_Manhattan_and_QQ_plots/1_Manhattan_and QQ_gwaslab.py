import gwaslab as gl
import matplotlib.font_manager
import pandas as pd

# Change the column names if neccessary
mysumstats = gl.Sumstats("PATH_TO_GWAS_SUMMARY_STATISTICS",
             snpid="MARKER",
             chrom="CHR",
             pos="BP",
             ea="A1",
             nea="A2",
             beta="Beta",
             se="SE",
             p="P",
             build="19")

mysumstats.plot_mqq(
                  mode="mqq",  # Figure layout: manhattan plot on left, QQ plot on right
                  anno="GENENAME", # Annotate gene names
                  sig_level=5e-8,
                  marker_size=(5,10),
                  figargs={"figsize":(15,4),"dpi":300},
                  save=f"Manhattan_QQ.png", 
                  save_args={"dpi":400,"facecolor":"white"})