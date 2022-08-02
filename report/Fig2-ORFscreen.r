library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

schema = function() {
    schema = grid::rasterGrob(magick::image_read("external/ORFscreen.svg"))
    p = ggplot() + annotation_custom(schema) + theme(panel.background=element_blank())
    wrap_elements(p)
}

volc = function() {
}

sys$run({
})
