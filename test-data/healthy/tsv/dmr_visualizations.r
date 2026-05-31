# =============================================================================
# DMR Visualization Script - TFM Envejecimiento cfDNA
# Regiones Diferencial­mente Metiladas (DMRs): Hyper e Hipometiladas
# Contrastes: MID vs YOUNG | OLD vs YOUNG | OLD vs MID
# =============================================================================
# Requisitos: tidyverse, ggplot2, patchwork, ggrepel
# Instalar si es necesario:
#   install.packages(c("tidyverse", "patchwork", "ggrepel", "scales"))
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)
library(scales)

# -----------------------------------------------------------------------------
# 0. CONFIGURACIÓN
# -----------------------------------------------------------------------------

# ── Rutas a los ficheros TSV ──────────────────────────────────────────────────
# Ajusta esta lista si los ficheros tienen nombres o rutas diferentes
FILES <- list(
  MID_vs_YOUNG_hyper = "C:Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/MID_vs_YOUNG_hyper.tsv",
  MID_vs_YOUNG_hypo  = "C:/Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/MID_vs_YOUNG_hypo.tsv",
  OLD_vs_YOUNG_hyper = "C:/Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/OLD_vs_YOUNG_hyper.tsv",
  OLD_vs_YOUNG_hypo  = "C:/Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/OLD_vs_YOUNG_hypo.tsv",
  OLD_vs_MID_hyper   = "C:/Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/OLD_vs_MID_hyper.tsv",
  OLD_vs_MID_hypo    = "C:/Users/bisite/Documents/GitHub/methyline/test-data/healthy/tsv/OLD_vs_MID_hypo.tsv"
)

# ── Directorio de salida ──────────────────────────────────────────────────────
OUT_DIR <- "figures_DMR"
dir.create(OUT_DIR, showWarnings = FALSE)

# ── Paleta de colores ─────────────────────────────────────────────────────────
COL_HYPER <- "#d94e4e"   # rojo: hipermetilación
COL_HYPO  <- "#4e8fd9"   # azul: hipometilación
COL_DIR   <- c(hyper = COL_HYPER, hypo = COL_HYPO)

COL_COMP  <- c(
  "MID vs YOUNG" = "#2a9d8f",
  "OLD vs YOUNG" = "#e76f51",
  "OLD vs MID"   = "#8338ec"
)

COMP_ORDER <- c("MID vs YOUNG", "OLD vs YOUNG", "OLD vs MID")

# ── Tema base para publicación ────────────────────────────────────────────────
theme_dmr <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(color = "grey40", size = base_size - 1),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1),
      legend.title     = element_text(face = "bold"),
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# ── Tamaño de exportación ─────────────────────────────────────────────────────
save_fig <- function(filename, plot, width = 8, height = 6, dpi = 300) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = dpi,
         bg = "white")
  message("✔  Guardado: ", path)
}


# =============================================================================
# 1. CARGA Y PREPARACIÓN DE DATOS
# =============================================================================

dmr_all <- imap_dfr(FILES, function(path, key) {
  if (!file.exists(path)) {
    warning("Fichero no encontrado: ", path); return(NULL)
  }
  parts      <- str_split(key, "_")[[1]]
  comparison <- paste(parts[1], "vs", parts[3])
  direction  <- parts[4]                       # "hyper" o "hypo"

  read_tsv(path, col_types = cols(), show_col_types = FALSE) |>
    mutate(
      comparison = comparison,
      direction  = direction,
      abs_diff   = abs(meth_diff),
      length_kb  = length / 1000,
      # Etiqueta corta del locus para anotaciones
      locus      = paste0(chr, ":", format(start, big.mark = ","))
    )
}) |>
  mutate(
    comparison = factor(comparison, levels = COMP_ORDER),
    direction  = factor(direction,  levels = c("hyper", "hypo"))
  )

message(
  "Datos cargados: ", nrow(dmr_all), " DMRs en total\n",
  paste(dmr_all |> count(comparison, direction) |>
          mutate(s = paste0(comparison, " [", direction, "]: ", n)) |>
          pull(s), collapse = "\n")
)


# =============================================================================
# 2. FIG 1 – Recuento de DMRs por contraste y dirección
# =============================================================================

counts <- dmr_all |> count(comparison, direction)

fig1 <- ggplot(counts, aes(x = comparison, y = n, fill = direction)) +
  geom_col(position = position_dodge(0.75), width = 0.65, color = "white") +
  geom_text(aes(label = n),
            position = position_dodge(0.75),
            vjust = -0.4, size = 4, fontface = "bold") +
  scale_fill_manual(values = COL_DIR,
                    labels  = c(hyper = "Hipermetiladas", hypo = "Hipometiladas")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "Número de DMRs por grupo",
    x        = NULL,
    y        = "Número de DMRs",
    fill     = NULL
  ) +
  theme_dmr()

save_fig("fig1_dmr_counts.png", fig1, width = 8, height = 5)


# =============================================================================
# 3. FIG 2 – Distribución de |Δβ| (boxplot agrupado)
# =============================================================================

fig2 <- ggplot(dmr_all, aes(x = comparison, y = abs_diff, fill = direction)) +
  geom_boxplot(
    outlier.shape  = 21,
    outlier.size   = 1.5,
    outlier.alpha  = 0.6,
    width          = 0.65,
    position       = position_dodge(0.75),
    color          = "grey30"
  ) +
  scale_fill_manual(values = COL_DIR,
                    labels  = c(hyper = "Hipermetiladas", hypo = "Hipometiladas")) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  labs(
    title    = "Magnitud del cambio de metilación por grupo",
    x        = NULL,
    y        = expression("|" * Delta * beta * "|"),
    fill     = NULL
  ) +
  theme_dmr()

save_fig("fig2_methdiff_boxplot.png", fig2, width = 8, height = 5)


# =============================================================================
# 4. FIG 3 – Scatter β_ref vs β_cond (diagonal de identidad)
# =============================================================================

# Identificar el Top 5 de DMRs por |areaStat| en cada grupo para anotar
top_dmr <- dmr_all |>
  group_by(comparison, direction) |>
  slice_max(abs(areaStat), n = 5) |>
  ungroup()

fig3 <- ggplot(dmr_all, aes(x = meanMethy1, y = meanMethy2,
                              color = comparison, shape = direction)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.6) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_text_repel(
    data        = top_dmr,
    aes(label   = locus),
    size        = 2.5,
    max.overlaps = 10,
    segment.color = "grey60",
    show.legend = FALSE
  ) +
  scale_color_manual(values = COL_COMP) +
  scale_shape_manual(values = c(hyper = 16, hypo = 17),
                     labels  = c(hyper = "Hipermetilada", hypo = "Hipometilada")) +
  scale_x_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.1)) +
  labs(
    title    = "β medio grupo referencia vs grupo condición por DMR",
    x        = "β medio — grupo referencia (meanMethy1)",
    y        = "β medio — grupo condición (meanMethy2)",
    color    = "Contraste",
    shape    = "Dirección"
  ) +
  theme_dmr() +
  theme(legend.position = "right")

save_fig("fig3_scatter_methy.png", fig3, width = 9, height = 7)


# =============================================================================
# 5. FIG 4 – Longitud de DMR (bp) por contraste y dirección (violín)
# =============================================================================

fig4 <- ggplot(dmr_all, aes(x = comparison, y = length, fill = direction)) +
  geom_violin(
    position = position_dodge(0.8),
    width    = 0.7,
    alpha    = 0.6,
    trim     = TRUE,
    color    = "grey30"
  ) +
  geom_boxplot(
    width    = 0.15,
    position = position_dodge(0.8),
    outlier.shape = NA,
    color    = "grey20",
    fill     = "white",
    alpha    = 0.9
  ) +
  scale_fill_manual(values = COL_DIR,
                    labels  = c(hyper = "Hipermetiladas", hypo = "Hipometiladas")) +
  scale_y_log10(labels = label_comma()) +
  labs(
    title    = "Longitud de las DMRs por contraste y dirección",
    x        = NULL,
    y        = "Longitud (bp, escala log₁₀)",
    fill     = NULL
  ) +
  theme_dmr()

save_fig("fig4_dmr_length_violin.png", fig4, width = 8, height = 5)


# =============================================================================
# 6. FIG 5 – nCpG medio en Top 10 DMRs más significativas (areaStat)
# =============================================================================

top10 <- dmr_all |>
  group_by(comparison, direction) |>
  slice_max(abs(areaStat), n = 10) |>
  summarise(
    mean_nCG   = mean(nCG),
    median_nCG = median(nCG),
    n          = n(),
    .groups    = "drop"
  )

fig5 <- ggplot(top10, aes(x = mean_nCG, y = comparison, fill = direction)) +
  geom_col(position = position_dodge(0.7), width = 0.6, color = "white") +
  geom_text(aes(label = round(mean_nCG, 1)),
            position  = position_dodge(0.7),
            hjust     = -0.15,
            size      = 3.8,
            fontface  = "bold") +
  scale_fill_manual(values = COL_DIR,
                    labels  = c(hyper = "Hipermetiladas", hypo = "Hipometiladas")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title    = "nCpGs medio en Top 10 DMRs más significativas",
    x        = "nCpG medio",
    y        = NULL,
    fill     = NULL
  ) +
  theme_dmr()

save_fig("fig5_ncpg_top10.png", fig5, width = 8, height = 4)


# =============================================================================
# 7. FIG 6 – Panel combinado (Fig1 + Fig2) con patchwork
# =============================================================================

fig6 <- (fig1 / fig2) +
  plot_annotation(
    title   = "Resumen de DMRs entre grupos de distinta edad",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 15),
      plot.caption = element_text(color = "grey50", size = 9)
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

save_fig("fig6_panel_combined.png", fig6, width = 9, height = 10)


# =============================================================================
# 8. TABLA RESUMEN (CSV para LaTeX / knitr)
# =============================================================================

summary_table <- dmr_all |>
  group_by(comparison, direction) |>
  summarise(
    n_DMR          = n(),
    median_abs_diff = median(abs_diff) |> round(3),
    mean_abs_diff   = mean(abs_diff)   |> round(3),
    sd_abs_diff     = sd(abs_diff)     |> round(3),
    median_nCG      = median(nCG)      |> round(1),
    median_length   = median(length)   |> round(0),
    .groups = "drop"
  ) |>
  arrange(comparison, direction)

write_csv(summary_table, file.path(OUT_DIR, "summary_DMR_table.csv"))
message("✔  Tabla resumen guardada: ", file.path(OUT_DIR, "summary_DMR_table.csv"))

print(summary_table)

message("\n✅  Script completado. Figuras en: ", OUT_DIR)
