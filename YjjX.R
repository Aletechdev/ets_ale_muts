library(trackViewer)

UNKNOWN_EFFECT = "#ABABAB"
CODING_SEQ_DISRUPT_COLOR = '#CF000F'
FUNC_DISRUPT_COLOR = "#ED7D31"
STRUCT_DISRUPT_COLOR = "#A020F0"
FUNC_AND_STRUCT_DISRUPT_COLOR = "#A52A2A"

feats = read.csv(file = './yjjX_aa_feats.csv', stringsAsFactors=FALSE)

aa_chain = feats[feats$feature=="Chain",]
feats = feats[feats$feature!="Chain",]

muts = read.csv(file = './yjjX_aa_muts.csv', stringsAsFactors=FALSE)

# Renaming experiment label according to that we want to specifically
# use within the manuscript.
muts$study[muts$study=="ndh-cydB-appC"] <- "ETC-4"

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
feats$layer[feats$feature=="Helix"] <- 0
feats$layer[feats$feature=="Beta strand"] <- 0
feats$layer[feats$feature=="Turn"] <- 0
feats$layer[feats$feature=="Manganese or magnesium mental binding site"] <- 1
feats$layer[feats$feature=="Substrate binding region"] <- 2
feats$layer[feats$feature=="YjjX subunit interface"] <- 3
features$featureLayerID <- paste(feats$layer)
names(features) <- paste(feats$feature)
features$fill <- feats$color
features$color <- feats$color

legend <- list(
  labels=c(
    "ETC-4",
    "unknown effect",
    # "deleterious (SIFT < 0.05)",
    # "deleterious (SIFT < 0.05) and structural disruption (ΔΔG > 2)"
    "truncation",
    "structurally destabilizing (ΔΔG > 2)"
  ),
  fill=c(
    "#4878D0",
    "white",
    "white",
    "white"
  ),
  col=c(
    "white",
    "black",
    CODING_SEQ_DISRUPT_COLOR,
    # FUNC_DISRUPT_COLOR,
    STRUCT_DISRUPT_COLOR
    # FUNC_AND_STRUCT_DISRUPT_COLOR
  )
)

lolliplot(mutations,
          features,
          legend = legend,
          ranges = GRanges("r", IRanges(aa_chain$start, aa_chain$end)),
          yaxis = FALSE
)
# grid.text("Mutations across YjjX's amino acid chain", x=.5, y=0.9, just="top", gp=gpar(cex=1.5, fontface="bold"))
