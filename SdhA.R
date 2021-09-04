library(trackViewer)

UNKNOWN_EFFECT = "#ABABAB"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

feats = read.csv(file = './sdhA_aa_feats.csv', stringsAsFactors=FALSE)

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]
feats = feats[feats$feature!="Turn",]
feats = feats[feats$feature!="Helix",]
feats = feats[feats$feature!="Beta strand",]
feats = feats[feats$feature!="Modified residue",]

muts = read.csv(file = './sdhA_aa_muts.csv', stringsAsFactors=FALSE)

muts$clr <- "black"
muts$clr <- muts$color

mutations <- GRanges("r", IRanges(muts$AA.pos, width=1, names=muts$name))
muts$border_color <- muts$clr
muts$border_color[muts$border_color==UNKNOWN_EFFECT] <- "black"
mutations$score <- muts$mutation.count
mutations$score <- muts$mutation.count
mutations$border <- muts$border_color
muts$text_colors = muts$clr
muts$text_colors[muts$text_colors==UNKNOWN_EFFECT] <- "black"
mutations$label.parameter.gp <- gpar(col=muts$text_colors)
mutations$dashline.col <- muts$clr
mutations$color <- muts$study.color

features <- GRanges("r", IRanges(feats$start, end=feats$end))
feats$height <- rep(0.03, length(features))
features$height = feats$height
# 
feats$layer <- 0
# feats$layer[feats$feature=="Helix"] <- 0
# feats$layer[feats$feature=="Beta strand"] <- 0
feats$layer[feats$feature=="FAD binding domain"] <- 0
feats$layer[feats$feature=="FAD-dependent pyridine nucleotide reductase"] <- 1
feats$layer[feats$feature=="FAD binding site"] <- 2
feats$layer[feats$feature=="Catalytic domain"] <- 1
feats$layer[feats$feature=="Active site"] <- 2
feats$layer[feats$feature=="Substrate binding site"] <- 3
feats$layer[feats$feature=="SdhB interface"] <- 5
features$featureLayerID <- paste(feats$layer)
names(features) <- paste(feats$feature)
features$fill <- feats$color
features$color <- feats$color

legend <- list(
  labels=c(
    "ETC-3",
    "SSW_GLU_XYL",
    "GLU",
    "TOL_hexamethylenediamine",
    "unknown effect",
    "truncation",
    # "deleterious (SIFT < 0.05)",
    # "structurally destabilizing (ΔΔG > 2)",
    "deleterious (SIFT < 0.05) and structural disruption (ΔΔG > 2)"
  ),
  fill=c(
    "#beaed4",
    "#FF9A00",
    "#fdc086",
    "#3875DB",
    "white",
    "white",
    "white"
    # "white",
    # "white"
  ),
  col=c(
    "white",
    "white",
    "white",
    "white",
    # "white",
    "black",
    CODING_SEQ_DISRUPT_COLOR,
    # FUNC_DISRUPT_COLOR,
    # STRUCT_DISRUPT_COLOR,
    FUNC_AND_STRUCT_DISRUPT_COLOR
  )
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$start, aa_chain$end)),
          yaxis = FALSE
)
# grid.text("Mutations across SdhA's amino acid chain", x=.5, y=0.65, just="top", gp=gpar(cex=1.5, fontface="bold"))
