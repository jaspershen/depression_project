library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1-code/100-tools.R')

setwd("3-data_analysis/demographic_data/")

data <- readr::read_csv("SchoolForTheWork-BaselineDemographics labels.csv")

colnames(data)

data <-
  data %>%
  dplyr::rename(subject_id = "Record ID", sex = Sex)

race <-
  1:nrow(data) %>%
  purrr::map(function(idx) {
    temp <-
      c(
        "Black",
        "Asian",
        "Native American",
        "White",
        "Native Hawaiian or other Pacific islander",
        "Other"
      )[which(unlist(as.data.frame(data[idx, 6:11])[1, , drop = TRUE]) == "Checked")]
    if (length(temp) == 0) {
      return(NA)
    }
    temp[1]
  }) %>%
  unlist()

length(race)
race[is.na(race)] <- "Unknown"

data$race <- race

library(circlize)

set.seed(999)
n = 78

df <-
  data.frame(
    factors = data$subject_id,
    x = 1,
    y = 1,
    data,
    stringsAsFactors = TRUE
  )

dim(df)

circos.par(
  "track.height" = 0.15,
  start.degree = 90,
  clock.wise = TRUE,
  gap.after = c(rep(0.02, 77), 20),
  cell.padding = c(0, 0, 0, 0)
)

circos.initialize(factors = df$factors,
                  x = df$x,
                  xlim = c(0.5, 1.5))

## sex
temp_col <- df$sex
temp_col[is.na(temp_col)] <- "grey"
temp_col[temp_col == "Male"] <- "#3B4992FF"
temp_col[temp_col == "Female"] <- "#EE0000FF"

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    circos.text(
      x = mean(xlim),
      y = 1.2,
      labels = name,
      facing = "bending.outside",
      niceFacing = TRUE,
      cex = 0.5,
      adj = aa
    )
    
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)

##bmi
df$Weight
df$Height

height <-
  df$Height %>% 
  sapply(convert2meters) %>% 
  unname()

bmi <- as.numeric(df$Weight) * 0.453592 / height^2

temp_value <- bmi

circos.track(
  factors = df$factors,
  # x = df$x,
  y = temp_value,
  ylim = range(temp_value, na.rm = TRUE),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.2,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.yaxis(
      side = "left",
      at = round(c(min(temp_value, na.rm = TRUE), (
        min(temp_value, na.rm = TRUE) + max(temp_value, na.rm = TRUE)
      ) / 2, max(temp_value, na.rm = TRUE)), 2),
      sector.index = get.all.sector.index()[1],
      labels.cex = 0.4,
      labels.niceFacing = FALSE
    )
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = bmi[i],
      col = ggsci::pal_aaas()(n = 10)[3],
      bg.border = "black"
    )
  }
)

## race
temp_col <- df$race
temp_col[temp_col == "Unknown"] <- "grey"
temp_col[temp_col == "White"] <- ggsci::pal_bmj()(n = 10)[1]
temp_col[temp_col == "Asian"] <- ggsci::pal_bmj()(n = 10)[2]
temp_col[temp_col == "Other"] <- ggsci::pal_bmj()(n = 10)[3]

circos.track(
  factors = df$factors,
  # x = df$x,
  y = df$y,
  ylim = c(0, 1),
  bg.border = "black",
  # bg.col = NA,
  track.height = 0.1,
  panel.fun = function(x, y) {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <-
      ifelse(theta < 90 ||
               theta > 270, "clockwise", "reverse.clockwise")
    aa = c(0.5, 1)
    # if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    #plot country labels
    circos.rect(
      xleft = xlim[1],
      ybottom = ylim[1],
      xright = xlim[2],
      ytop = ylim[2],
      col = temp_col[i],
      bg.border = "black"
    )
  }
)


# ##mother age
# range(df$mother_age)
# temp_value <- df$mother_age
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = temp_value,
#   ylim = range(temp_value),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.08,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.yaxis(
#       side = "left",
#       at = c(min(temp_value), round((
#         min(temp_value) + max(temp_value)
#       ) / 2, 2), max(temp_value)),
#       sector.index = get.all.sector.index()[1],
#       labels.cex = 0.4,
#       labels.niceFacing = FALSE
#     )
#     
#     circos.lines(
#       x = mean(xlim),
#       y =  temp_value[i],
#       pch = 16,
#       cex = 8,
#       type = "h",
#       col = ggsci::pal_aaas()(n = 10)[4],
#       lwd = 2
#     )
#     
#     circos.points(
#       x = mean(xlim),
#       y =  temp_value[i],
#       pch = 16,
#       cex = 1,
#       col = ggsci::pal_aaas()(n = 10)[4]
#     )
#   }
# )
# 
# ##parity
# temp_col <- info$parity
# temp_col[temp_col == 1] <- ggsci::pal_aaas(alpha = 0.3)(n = 10)[4]
# temp_col[temp_col == 2] <- ggsci::pal_aaas(alpha = 0.65)(n = 10)[4]
# temp_col[temp_col == 3] <- ggsci::pal_aaas(alpha = 1)(n = 10)[4]
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_col[i],
#       bg.border = "black"
#     )
#   }
# )
# 
# ##vaginal birth
# temp_col <- info$vaginal_birth
# temp_col[temp_col == 0] <- ggsci::pal_aaas()(n = 10)[1]
# temp_col[temp_col == 1] <- ggsci::pal_aaas()(n = 10)[2]
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#   }
# )
# 
# ##labor onset
# temp_col <- info$labor_onset
# temp_col[temp_col == "natural"] <- ggsci::pal_aaas()(n = 10)[3]
# temp_col[temp_col == "advanced"] <- ggsci::pal_aaas()(n = 10)[6]
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_col[i],
#       bg.border = "black"
#     )
#   }
# )
# 
# ##induce labor
# temp_col <- info$induced_labour
# temp_col[temp_col == 1] <- ggsci::pal_aaas()(n = 10)[2]
# temp_col[temp_col == 0] <- ggsci::pal_aaas()(n = 10)[1]
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_col[i],
#       bg.border = "black"
#     )
#   }
# )
# 
# ##child length
# range(temp_value)
# temp_value <- df$child_length
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = temp_value,
#   ylim = range(temp_value),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.08,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.yaxis(
#       side = "left",
#       at = c(min(temp_value), round((
#         min(temp_value) + max(temp_value)
#       ) / 2, 2), max(temp_value)),
#       sector.index = get.all.sector.index()[1],
#       labels.cex = 0.4,
#       labels.niceFacing = FALSE
#     )
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = temp_value[i],
#       col = ggsci::pal_aaas()(n = 10)[5],
#       bg.border = "black"
#     )
#   }
# )
# 
# 
# ##child weight
# 
# temp_value <- df$child_weight
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = temp_value,
#   ylim = range(temp_value),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.08,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.yaxis(
#       side = "left",
#       at = c(min(temp_value), round((
#         min(temp_value) + max(temp_value)
#       ) / 2, 2), max(temp_value)),
#       sector.index = get.all.sector.index()[1],
#       labels.cex = 0.4,
#       labels.niceFacing = FALSE
#     )
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = temp_value[i],
#       col = ggsci::pal_aaas()(n = 10)[6],
#       bg.border = "black"
#     )
#   }
# )
# 
# ##child sex
# temp_col <- info$child_sex
# temp_col[temp_col == "M"] <- ggsci::pal_aaas()(n = 10)[10]
# temp_col[temp_col == "F"] <- ggsci::pal_aaas()(n = 10)[2]
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_col[i],
#       bg.border = "black"
#     )
#   }
# )
# 
# ###g stage
# temp_value <- df$g_stage
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = temp_value,
#   ylim = range(temp_value),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.08,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.lines(
#       x = mean(xlim),
#       y =  temp_value[i],
#       pch = 16,
#       cex = 8,
#       type = "h",
#       col = ggsci::pal_aaas()(n = 10)[7],
#       lwd = 2
#     )
#     
#     circos.yaxis(
#       side = "left",
#       at = c(round(min(temp_value), 2), round((
#         min(temp_value) + max(temp_value)
#       ) / 2, 2), round(max(temp_value), 2)),
#       sector.index = get.all.sector.index()[1],
#       labels.cex = 0.4,
#       labels.niceFacing = FALSE
#     )
#     
#     circos.points(
#       x = mean(xlim),
#       y =  temp_value[i],
#       pch = 16,
#       cex = 1,
#       col = ggsci::pal_aaas()(n = 10)[7]
#     )
#   }
# )
# 
# ##begin season
# temp_col <- info$begin_season
# temp_col[temp_col == "spring"] <- ggsci::pal_aaas()(n = 10)[3]
# temp_col[temp_col == "summer"] <- ggsci::pal_aaas()(n = 10)[2]
# temp_col[temp_col == "autumm"] <- ggsci::pal_d3()(10)[2]
# temp_col[temp_col == "winter"] <- "white"
# 
# circos.track(
#   factors = df$factors,
#   # x = df$x,
#   y = df$y,
#   ylim = c(0, 1),
#   bg.border = "black",
#   # bg.col = NA,
#   track.height = 0.04,
#   panel.fun = function(x, y) {
#     name = get.cell.meta.data("sector.index")
#     i = get.cell.meta.data("sector.numeric.index")
#     xlim = get.cell.meta.data("xlim")
#     ylim = get.cell.meta.data("ylim")
#     
#     circos.rect(
#       xleft = xlim[1],
#       ybottom = ylim[1],
#       xright = xlim[2],
#       ytop = ylim[2],
#       col = temp_col[i],
#       bg.border = "black"
#     )
#   }
# )
# 
# circos.clear()
# 