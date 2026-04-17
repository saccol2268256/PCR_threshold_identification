library(tidyverse)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(pROC)
library(lme4)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(scales)
library(glmmTMB)

#' Impostare a TRUE per ESCLUDERE le osservazioni di Breganze 2 (area == 1).
ESCLUDI_BREGANZE2 <- TRUE

#' Impostare a TRUE per ESCLUDERE le osservazioni senza determinazione del sesso.
ESCLUDI_NO_SESSO <- F

#' Soglia QC per qualità estrazione DNA totale (cq_mt16s).
QC_THRESHOLD <- 32
MIN_OBS_PER_TEMPO <- 15

setwd("C:/Users/sacco/Desktop/lavoro/25/estrazione")
df_raw <- read.csv("dataset/estrazione_clean.csv", header = TRUE, sep = ",")

#' Conversione tipi
df_raw$time     <- as.Date(df_raw$time)
df_raw$sesso    <- factor(df_raw$sesso)
df_raw$vigneto  <- factor(df_raw$vigneto)
df_raw$trappola <- factor(df_raw$trappola)

#' Variabile temporale continua
df_raw <- df_raw %>%
  arrange(time) %>%
  mutate(
    times_elapsed = as.integer(time - min(time, na.rm = TRUE)),
    tempo         = times_elapsed + 1
  )

#' Etichette leggibili
df_raw <- df_raw %>%
  mutate(
    vigneto_label = case_when(
      vigneto == 0 & area == 0 ~ "Breganze 1",
      vigneto == 0 & area == 1 ~ "Breganze 2",
      vigneto == 1             ~ "Cimone",
      TRUE                     ~ "NA"
    ),
    vigneto_area = interaction(vigneto, area, sep = "_", drop = TRUE)
  )

#' Flag QC
df_raw <- df_raw %>%
  mutate(
    qc_status = case_when(
      is.na(cq_mt16s) | is.na(cq_sfd) ~ "Missing",
      cq_mt16s > QC_THRESHOLD          ~ "ND (Poor Quality)",
      sfd == 1                         ~ "Positive (QC Pass)",
      sfd == 0                         ~ "Negative (QC Pass)",
      TRUE                             ~ "Unknown"
    ),
    qc_status = factor(qc_status, levels = c("ND (Poor Quality)", 
                                             "Negative (QC Pass)",
                                             "Positive (QC Pass)", 
                                             "Missing"))
  )

#' Applicazione filtri
df <- df_raw

if (ESCLUDI_BREGANZE2) {
  n_prima <- nrow(df)
  df <- df %>% filter(!(vigneto == 0 & area == 1))
  message(sprintf("[FILTRO B2] Rimosse %d osservazioni di Breganze 2. Rimaste: %d (%.1f%%)",
                  n_prima - nrow(df), nrow(df), nrow(df) / n_prima * 100))
}

df_sex <- df %>% filter(!is.na(sesso))

if (ESCLUDI_NO_SESSO) {
  n_prima <- nrow(df)
  df <- df_sex
  message(sprintf("[FILTRO SESSO] Rimosse %d osservazioni senza sesso noto. Rimaste: %d (%.1f%%)",
                  n_prima - nrow(df), nrow(df), nrow(df) / n_prima * 100))
}

#' Diagnostica copertura temporale post-filtro
copertura_tempo <- df %>%
  group_by(time, vigneto_label) %>%
  summarise(n_obs = n(), n_pos = sum(sfd == 1, na.rm = TRUE),
            pct_pos = n_pos / sum(!is.na(sfd)) * 100,
            warn_n = n_obs < MIN_OBS_PER_TEMPO, .groups = "drop")

if (sum(copertura_tempo$warn_n, na.rm = TRUE) > 0) {
  message(sprintf("[WARNING] %d combinazioni data × vigneto con N < %d:",
                  sum(copertura_tempo$warn_n, na.rm = TRUE), MIN_OBS_PER_TEMPO))
  print(copertura_tempo %>% filter(warn_n))
} else {
  message(sprintf("[OK] Tutte le combinazioni data × vigneto hanno almeno %d osservazioni.", 
                  MIN_OBS_PER_TEMPO))
}

#' Grafico copertura temporale
ggplot(copertura_tempo, aes(x = time, y = vigneto_label, fill = n_obs)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = n_obs, color = if_else(warn_n, "red", "white")),
            size = 3, fontface = "bold") +
  scale_fill_gradient(low = "#deebf7", high = "#08306b", name = "N osservazioni") +
  scale_color_identity() +
  labs(title = "Copertura Temporale Post-Filtro",
       subtitle = sprintf("Rosso = N < %d (attenzione)", MIN_OBS_PER_TEMPO),
       x = "Data", y = "Vigneto") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))

#' Dataset QC-pass
df_qc <- df %>% filter(!is.na(cq_mt16s) & cq_mt16s <= QC_THRESHOLD)

# ANALISI DESCRITTIVA E TEST STATISTICI==============================================================================

# Prevalenza SFD
sfd_vig <- df %>% 
  mutate(vig_lab = if_else(vigneto == 1, "Cimone", "Breganze")) %>%
  group_by(vig_lab) %>%
  summarise(n = n(), pos = sum(sfd == 1, na.rm = TRUE),
            pct = pos / sum(!is.na(sfd)) * 100, .groups = "drop")
print(sfd_vig)
#' I due vigneti presentano una differenza statisticamente significativa nella pressione della malattia, con una prevalenza nettamente superiore a Cimone rispetto a Breganze.

sfd_area <- df %>% filter(vigneto == 0) %>%
  group_by(area) %>%
  summarise(n = n(), pos = sum(sfd == 1, na.rm = TRUE),
            pct = pos / sum(!is.na(sfd)) * 100, .groups = "drop") %>%
  mutate(area_label = c("Breganze 1", "Breganze 2")[area + 1])
print(sfd_area)
#' All’interno del vigneto di Breganze si osserva una forte disomogeneità tra le due aree, con Breganze 2 che mostra una prevalenza sensibilmente più elevata.

sfd_ses <- df %>% filter(!is.na(sesso)) %>%
  group_by(sesso) %>%
  summarise(n = n(), pos = sum(sfd == 1, na.rm = TRUE),
            pct = pos / sum(!is.na(sfd)) * 100, .groups = "drop")
print(sfd_ses)
#' A livello globale, il sesso maschile presenta una prevalenza di SFD superiore di circa 2 punti percentuali rispetto al sesso femminile.

# Test Welch su proporzioni
welch_vigneto <- t.test(sfd ~ if_else(vigneto == 1, "Cimone", "Breganze"),
                        data = df %>% filter(!is.na(sfd)), var.equal = FALSE)
welch_sesso   <- t.test(sfd ~ sesso, data = df %>% filter(!is.na(sesso), !is.na(sfd)), var.equal = FALSE)
welch_sesso_B <- t.test(sfd ~ sesso, data = df %>% filter(vigneto == 0, !is.na(sesso), !is.na(sfd)), var.equal = FALSE)
welch_sesso_C <- t.test(sfd ~ sesso, data = df %>% filter(vigneto == 1, !is.na(sesso), !is.na(sfd)), var.equal = FALSE)

tibble(Confronto = c("Cimone vs Breganze", "M vs F (globale)", 
                     "M vs F (Breganze)", "M vs F (Cimone)"),
       t_stat = c(welch_vigneto$statistic, welch_sesso$statistic,
                  welch_sesso_B$statistic, welch_sesso_C$statistic),
       df_welch = c(welch_vigneto$parameter, welch_sesso$parameter,
                    welch_sesso_B$parameter, welch_sesso_C$parameter),
       p_value = c(welch_vigneto$p.value, welch_sesso$p.value,
                   welch_sesso_B$p.value, welch_sesso_C$p.value)) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)),
         sig = case_when(p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
                         p_value < 0.05 ~ "*", TRUE ~ "ns")) %>%
  print()
#' I vigneti differiscono significativamente (p < 0,05) con una pressione della malattia maggiore a Cimone. A livello globale il sesso maschile mostra una prevalenza superiore (p < 0,10). Limitando l’analisi al solo vigneto di Cimone, la differenza tra maschi e femmine diventa significativa al 5%.

# Trend temporale % positivi SFD
sfd_time_vigneto <- df %>%
  group_by(time, vigneto_label) %>%
  summarise(n = n(), pos = sum(sfd == 1, na.rm = TRUE),
            pct = pos / sum(!is.na(sfd)) * 100, .groups = "drop")

ggplot(sfd_time_vigneto, aes(x = time, y = pct, color = vigneto_label, fill = vigneto_label)) +
  geom_line(linewidth = 0.8) + geom_point(aes(size = n), alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.15, linewidth = 1, span = 1.5) +
  scale_color_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1", "Cimone" = "#d73027")) +
  scale_fill_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1", "Cimone" = "#d73027")) +
  scale_size_continuous(range = c(2, 6), name = "N campioni") +
  labs(title = "Trend % Positivi SFD per Vigneto", x = "Data", y = "% Positivi SFD") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
#' Il trend temporale evidenzia una diffusione più rapida e intensa dell’infezione a Cimone, soprattutto nella fase iniziale della stagione, confermando una maggiore pressione della malattia in questo vigneto.

# Boxplot variabilità temporale per trappola
df %>%
  group_by(trappola, vigneto_label, time) %>%
  summarise(pct = sum(sfd == 1, na.rm = TRUE) / sum(!is.na(sfd)) * 100, .groups = "drop") %>%
  ggplot(aes(x = reorder(trappola, pct, FUN = median), y = pct, fill = vigneto_label)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.4, show.legend = FALSE) +
  facet_wrap(~ vigneto_label, scales = "free_x") +
  scale_fill_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1", "Cimone" = "#d73027")) +
  labs(title = "Variabilità Temporale % Positivi per Trappola",
       x = "Trappola (mediana)", y = "% Positivi SFD (per data)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))
#' Il grafico evidenzia una marcata disomogeneità spaziale all’interno di entrambi i vigneti: numerose trappole mostrano percentuali di positivi strutturalmente basse o nulle, confermando l’esistenza di cluster di infezione localizzati.

# Top 10 trappole per % positivi
top_trappole <- df %>%
  group_by(trappola, vigneto_label) %>%
  summarise(n = n(), positivi = sum(sfd == 1, na.rm = TRUE),
            pct = positivi / sum(!is.na(sfd)) * 100, .groups = "drop") %>%
  arrange(desc(pct))

top_trappole %>% slice_head(n = 10) %>%
  ggplot(aes(x = reorder(paste0("T", trappola, "  [", vigneto_label, "]"), pct),
             y = pct, fill = vigneto_label)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(pct, 1), "%  n=", n)), hjust = -0.1, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1", "Cimone" = "#d73027")) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Top 10 Trappole per % Positivi SFD", x = "Trappola [Vigneto]", y = "% Positivi SFD") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
#' I risultati confermano una forte variabilità spaziale: le trappole più positive si concentrano prevalentemente nel vigneto di Cimone, mentre a Breganze prevalgono trappole con bassa o nulla prevalenza.

#' QC Gate: scatter cq_mt16s vs cq_sfd con soglia fissa di laboratorio (classificazione originale).
#' I tre pannelli confrontano: (A) qc_status di laboratorio, (B) classificazione Otsu su delta_cq,
#' (C) classificazione Youden su delta_cq.
p_qc_orig <- df %>%
  filter(!is.na(cq_mt16s), !is.na(cq_sfd)) %>%
  # Rimosso "df" da ggplot() perché arriva già dalla pipe
  ggplot(aes(x = cq_mt16s, y = cq_sfd, color = qc_status, shape = qc_status)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(
    values = c("ND (Poor Quality)" = "#d73027",
               "Negative (QC Pass)" = "#4575b4",
               "Positive (QC Pass)" = "#fee090",
               "Missing" = "#cccccc"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("ND (Poor Quality)" = 4, 
               "Negative (QC Pass)" = 16,
               "Positive (QC Pass)" = 17, 
               "Missing" = 15),
    name = NULL
  ) +
  labs(
    title = "(A) QC Gate – Soglia Laboratorio",
    subtitle = paste0("Soglia CQ_MT16S = ", QC_THRESHOLD),
    x = "CQ_MT16S", 
    y = "CQ_SFD"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )



# INDICI DI CARICA RELATIVA: RAPPORTO_CQ E DELTA_CQ ===========================
#' Analisi approfondita per identificazione soglia dinamica valida.
#' Rapporto_CQ e Delta_CQ come indici di carica relativa di fitoplasma.
#'Utilizziamo dopo aver provato analisi su entrambi perche miglior identificazione di coda destra della distribuzione,
#'parte di modellistica porta a simili risultati, ma data indentificazione di una soglia con ottimo AUC>90 e  miglior interpretazione biologica utilizziamo DELTA_CQ

df_rapporto <- df %>%
  filter(!is.na(rapporto_cq), !is.na(delta_cq), !is.na(vigneto), cq_mt16s <= QC_THRESHOLD)

df_rapporto %>%
  group_by(vigneto_label) %>%
  summarise(n = n(), media_rap = mean(rapporto_cq), sd_rap = sd(rapporto_cq),
            media_delta = mean(delta_cq), sd_delta = sd(delta_cq), .groups = "drop") %>%
  print()

ggplot(df_rapporto, aes(x = vigneto_label, y = delta_cq,
                        fill = vigneto_label, color = vigneto_label)) +
  geom_violin(alpha = 0.45, show.legend = FALSE) +
  geom_boxplot(width = 0.12, alpha = 0.8, color = "black",
               outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1.5, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 5,
               color = "red", show.legend = FALSE) +
  scale_fill_manual(values  = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027")) +
  scale_color_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027")) +
  labs(title = "Delta CQ per Vigneto",
       subtitle = "Valori bassi/negativi = maggiore carico relativo fitoplasma | Rosso = media",
       x = "Vigneto / Area", y = "Delta CQ (Cq_SFD - Cq_MT16S)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# SOGLIA OTTIMALE OTSU SU DELTA_CQ ==============================

otsu_cutoff <- function(x, n_breaks = 500) {
  x   <- x[!is.na(x)]
  rng <- range(x)
  candidati  <- seq(rng[1], rng[2], length.out = n_breaks)
  var_intra  <- sapply(candidati, function(c) {
    g1 <- x[x <= c]; g2 <- x[x > c]
    if (length(g1) < 2 | length(g2) < 2) return(Inf)
    w1 <- length(g1) / length(x);  w2 <- 1 - w1
    w1 * var(g1) + w2 * var(g2)
  })
  c_star <- candidati[which.min(var_intra)]
  list(cutoff = c_star,
       separabilita = 1 - min(var_intra) / var(x),
       n_sotto = sum(x <= c_star),
       n_sopra = sum(x > c_star))
}

#' --- Otsu su delta_cq ---
otsu_delta <- otsu_cutoff(df_qc$delta_cq)
message(sprintf("[OTSU delta_cq] Soglia: %.4f | eta^2: %.3f | n<=: %d, n>: %d",
                otsu_delta$cutoff, otsu_delta$separabilita,
                otsu_delta$n_sotto, otsu_delta$n_sopra))

ggplot(df_qc %>% filter(!is.na(delta_cq)),
       aes(x = delta_cq, fill = vigneto_label)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.5, position = "identity") +
  geom_density(aes(color = vigneto_label), linewidth = 1, alpha = 0) +
  geom_vline(xintercept = otsu_delta$cutoff, linetype = "dashed", color = "black", linewidth = 1.2) +
  annotate("text", x = otsu_delta$cutoff, y = Inf,
           label = sprintf(" Otsu = %.3f\n eta^2 = %.3f", otsu_delta$cutoff, otsu_delta$separabilita),
           hjust = 0, vjust = 1.4, fontface = "bold", size = 3.5) +
  scale_fill_manual(values  = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027"), name = "Vigneto") +
  scale_color_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027"), name = "Vigneto") +
  labs(title    = "Otsu Cut-off su Delta_CQ",
       subtitle = "La linea tratteggiata minimizza la varianza intra-classe",
       x = "Delta CQ (CqsFD - Cqmt16S)", y = "Densita'") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))


#' Concordanza classificazione sfd originale vs soglia Otsu su delta_cq
#' vediamo come le etichette date dalla contro negatività dopo l'analisi PCR venga rispettata anche dalla soglia 
df_qc <- df_qc %>%
  mutate(
    sfd_otsu_delta = if_else(delta_cq <= otsu_delta$cutoff, 1L, 0L)
  )
print(table(sfd_sFD = df_qc$sfd, sfd_otsu_delta = df_qc$sfd_otsu_delta))


# SOGLIA COMBINATA OTSU + ROC (YOUDEN J) =======================================
#' Confronto tra soglia Otsu e soglia ottimale via J di Youden (pROC).
#' J = sensibilità + specificità - 1. Il valore massimo identifica il cut-off
#' che bilancia al meglio vera positività e vera negatività(etichette di partenza).

df_roc_del <- df_qc %>% filter(!is.na(delta_cq),    !is.na(sfd))
roc_del   <- pROC::roc(df_roc_del$sfd, -df_roc_del$delta_cq,    quiet = TRUE)

#' coords restituisce il cut-off su -x, bisogna invertire il segno
best_del  <- pROC::coords(roc_del, "best", ret = c("threshold","sensitivity","specificity"),
                          best.method = "youden", transpose = FALSE)

#' Cut-off ROC sulla scala originale (invertiamo il segno)
youden_cutoff_del <- -as.numeric(best_del$threshold[1])

message(sprintf(
  "[SOGLIA COMBINATA delta_cq]    Otsu = %.4f | Youden = %.4f | sens = %.3f | spec = %.3f",
  otsu_delta$cutoff, youden_cutoff_del,
  as.numeric(best_del$sensitivity[1]), as.numeric(best_del$specificity[1])
))

#' Aggiunta classificazione Youden al dataset
df_qc <- df_qc %>%
  mutate(
    sfd_youden_delta = if_else(delta_cq    <= youden_cutoff_del, 1L, 0L)
  )

#' CONFRONTO VISIVO: QC Gate originale vs classificazioni Otsu e Youden su delta_cq.
df_qc_scatter <- df_qc %>%
  filter(!is.na(cq_mt16s), !is.na(cq_sfd)) %>%
  mutate(
    classe_otsu   = factor(sfd_otsu_delta,   levels = c(0,1), labels = c("Neg (Otsu)",   "Pos (Otsu)")),
    classe_youden = factor(sfd_youden_delta, levels = c(0,1), labels = c("Neg (Youden)", "Pos (Youden)"))
  )

col_neg <- "#4575b4"; col_pos <- "#d73027"

#' Pannello B – Otsu
p_qc_otsu <- ggplot(df_qc_scatter, aes(x = cq_mt16s, y = cq_sfd, color = classe_otsu)) +
  geom_abline(slope = 1, intercept = otsu_delta$cutoff,
              linetype = "dashed", color = "black", linewidth = 1.1) +
  geom_point(size = 2, alpha = 0.6) +
  annotate("text",
           x = min(df_qc_scatter$cq_mt16s, na.rm = TRUE) + 0.5,
           y = min(df_qc_scatter$cq_mt16s, na.rm = TRUE) + 0.5 + otsu_delta$cutoff + 1,
           label = sprintf("delta = %.3f (Otsu)", otsu_delta$cutoff),
           hjust = 0, size = 3, fontface = "italic") +
  scale_color_manual(values = c("Neg (Otsu)" = col_neg, "Pos (Otsu)" = col_pos),
                     name = NULL) +
  labs(title = "(B) Classificazione Otsu su Delta_CQ",
       x = "CQ_MT16S", y = "CQ_SFD") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

#' Pannello C – Youden
p_qc_youden <- ggplot(df_qc_scatter, aes(x = cq_mt16s, y = cq_sfd, color = classe_youden)) +
  geom_abline(slope = 1, intercept = youden_cutoff_del,
              linetype = "dotdash", color = "darkred", linewidth = 1.1) +
  geom_point(size = 2, alpha = 0.6) +
  annotate("text",
           x = min(df_qc_scatter$cq_mt16s, na.rm = TRUE) + 0.5,
           y = min(df_qc_scatter$cq_mt16s, na.rm = TRUE) + 0.5 + youden_cutoff_del + 1,
           label = sprintf("delta = %.3f (Youden)", youden_cutoff_del),
           hjust = 0, size = 3, fontface = "italic", color = "darkred") +
  scale_color_manual(values = c("Neg (Youden)" = col_neg, "Pos (Youden)" = col_pos),
                     name = NULL) +
  labs(title = "(C) Classificazione Youden su Delta_CQ",
       x = "CQ_MT16S", y = "CQ_SFD") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

#' Figura composita: pannello originale (A) + Otsu (B) + Youden (C)
gridExtra::grid.arrange(p_qc_orig, p_qc_otsu, p_qc_youden, ncol = 3)

#' Curva ROC con annotazione del punto Youden
pROC::plot.roc(roc_del, print.auc = TRUE, asp = FALSE,
               main = "ROC Delta_CQ (asse invertito)")
legend("bottomright",
       legend = sprintf("Youden = %.3f", youden_cutoff_del),
       col = "darkgreen", pch = 19, bty = "n")

#'Otteniamo un ottima soglia che va a classificare con un AUC del 92% per i soli campioni con sesso caratterizzato di conformita con le etichette da voi corrisposte,
#'mantenendo però la migliore inclinazione di delta_cq data dalla identificazione della combinazione delle 2 variabili logaritmiche rilevate.

#' Confrontare le distribuzioni con classificazione per valutare la separazione in base alle due soglie.

df_otsu_plot <- df_qc %>%
  filter(!is.na(rapporto_cq), !is.na(delta_cq)) %>%
  mutate(
    classe_otsu_delta = factor(sfd_otsu_delta, levels = c(0, 1),
                               labels = c("Neg (Otsu)", "Pos (Otsu)"))
  )

#' -- Delta_CQ con colore classe Otsu --
ggplot(df_otsu_plot,
       aes(x = vigneto_label, y = delta_cq, fill = vigneto_label)) +
  geom_violin(alpha = 0.3, show.legend = FALSE) +
  geom_boxplot(width = 0.10, alpha = 0.7, color = "black",
               outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(aes(color = classe_otsu_delta, shape = classe_otsu_delta),
              width = 0.12, alpha = 0.55, size = 1.8) +
  geom_hline(yintercept = otsu_delta$cutoff,  linetype = "dashed",
             color = "black",   linewidth = 1, alpha = 0.8) +
  geom_hline(yintercept = youden_cutoff_del,  linetype = "dotdash",
             color = "darkred", linewidth = 1, alpha = 0.8) +
  annotate("text", x = 0.5, y = otsu_delta$cutoff,
           label = sprintf(" Otsu=%.3f",   otsu_delta$cutoff),
           hjust = 0, vjust = -0.4, size = 3, fontface = "italic") +
  annotate("text", x = 0.5, y = youden_cutoff_del,
           label = sprintf(" Youden=%.3f", youden_cutoff_del),
           hjust = 0, vjust = 1.3, size = 3, fontface = "italic", color = "darkred") +
  scale_fill_manual(values  = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027")) +
  scale_color_manual(values = c("Neg (Otsu)" = "#4575b4", "Pos (Otsu)" = "#d73027"),
                     name = "Classe Otsu") +
  scale_shape_manual(values = c("Neg (Otsu)" = 16, "Pos (Otsu)" = 17),
                     name = "Classe Otsu") +
  labs(title    = "Delta_CQ per Vigneto – Categorie Otsu",
       subtitle = "Tratteggio nero = soglia Otsu | Punto-linea rosso = soglia Youden",
       x = "Vigneto / Area", y = "Delta CQ (CqsFD - Cqmt16S)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")



# ANDAMENTO CARICA FITOPLASMATICA NEL TEMPO =====================================
#' Distribuzione di cq_sfd
ggplot(df_qc, aes(x = factor(time), y = cq_sfd,
                  fill = vigneto_label, color = vigneto_label)) +
  # Boxplot per ogni punto temporale
  geom_boxplot(alpha = 0.8, color = "black", outlier.shape = NA, show.legend = FALSE) +
  # Jitter per vedere i singoli punti
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.2, show.legend = FALSE) +
  # Punto rosso per la media
  stat_summary(fun = mean, geom = "point", shape = 18,
               size = 3, color = "red", show.legend = FALSE) +
  # Suddivide il grafico in pannelli per ogni vigneto (es. incolonnati verticalmente)
  facet_wrap(~ vigneto_label, ncol = 1) + 
  scale_fill_manual(values  = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027")) +
  scale_color_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027")) +
  labs(title    = "Andamento Temporale CQ_SFD per Vigneto:Area (positivi QC-pass)",
       subtitle = "Valori bassi = maggiore carica fitoplasma | Rosso = media",
       x = "Data di Campionamento",
       y = "CQ_SFD  <- piu' fitoplasma | meno fitoplasma ->") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        # Ruota le etichette dell'asse X a 45 gradi per una migliore lettura
        axis.text.x = element_text(angle = 45, hjust = 1))

#' Analisi della carica di fitoplasma nel tempo tramite Delta_CQ e Rapporto_CQ.
#' Interpretazione delle pendenze temporali:
#'   beta_tempo > 0  ->  Delta sale nel tempo    ->  carica DIMINUISCE o stabile


df_carica <- df_qc %>% filter(sfd == 1, !is.na(delta_cq), !is.na(rapporto_cq))

#' --- Regressione lineare semplice: indice ~ tempo ---
lm_delta_glob   <- lm(delta_cq   ~ tempo,               data = df_carica)
summary(lm_delta_glob)

ggplot(df_carica, aes(x = tempo, y = delta_cq,
                      color = vigneto_label, fill = vigneto_label)) +
  geom_point(alpha = 0.35, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2, linewidth = 1.2) +
  scale_color_manual(values = c("Breganze 1" = "#2166ac", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027"), name = "Vigneto") +
  scale_fill_manual(values  = c("Breganze 1" = "#2166ac", "Breganze 2" = "#74add1",
                                "Cimone" = "#d73027"), name = "Vigneto") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf(" beta=%.4f  p=%.4f",
                           coef(lm_delta_glob)["tempo"],
                           summary(lm_delta_glob)$coef["tempo","Pr(>|t|)"]),
           hjust = 0, vjust = 1.5, size = 3.5, fontface = "italic") +
  labs(title    = "Delta_CQ nei Positivi nel Tempo",
       subtitle = "beta < 0 -> carica fitoplasmatica crescente nel tempo",
       x = "Giorni dall'inizio stagione", y = "Delta_CQ (CqsFD - Cqmt16S)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")

#' Modelli lineari misti con effetto casuale trappola
#' Intercetta casuale per trappola: cattura differenze di base tra trappole.
lmm_delta_trap <- lmer(delta_cq    ~ tempo +vigneto +(1 | trappola),
                       data = df_carica, REML = TRUE)
summary(lmm_delta_trap)




# MODELLO ZERO-INFLATED (glmmTMB) E GLMM PER PREVALENZA =========================
#' Motivazione zero-inflation: quando molte trappole/date mostrano 0 positivi
#' strutturali (zone sistematicamente indenni), un semplice GLMM binomiale
#' sovrastima la probabilita' di positivo nei cluster negativi.
#' 
#'Verifichiamo se la pendenza della retta cambia da Luglio a Settembre. Se la pendenza cambia, la soglia di classificazione deve essere dinamica perché la biologia dell'insetto è cambiata

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

df_glmm <- df %>%
  filter(!is.na(sfd), !is.na(vigneto), !is.na(trappola)) %>%
  mutate(
    vigneto_f  = factor(vigneto, levels = c(0, 1), labels = c("Breganze", "Cimone")),
    trappola_f = factor(trappola)
  )
df_glmm_sex <- df_glmm %>% filter(!is.na(sesso))

#' Dataset aggregato pool-level
df_pool <- df %>%
  filter(!is.na(sfd)) %>%
  mutate(vigneto_f = factor(vigneto, levels = c(0, 1), labels = c("Breganze", "Cimone"))) %>%
  group_by(trappola, time, tempo, vigneto_f, vigneto_label) %>%
  summarise(n_pos = sum(sfd == 1, na.rm = TRUE),
            n_neg = sum(sfd == 0, na.rm = TRUE),
            n_tot = n(),
            perc_pos = n_pos / n_tot * 100,
            .groups = "drop")


#' perc_positivi per pool (trappola x data): strumento per il campionamento parziale.
df %>%
  group_by(trappola, time, vigneto_label) %>%
  summarise(n_anal = n(), n_pos = sum(sfd == 1, na.rm = TRUE),
            perc_pos = n_pos / n_anal * 100, .groups = "drop") %>%
  ggplot(aes(x = perc_pos, fill = vigneto_label)) +
  geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
  facet_wrap(~ vigneto_label) +
  scale_fill_manual(values = c("Breganze 1" = "#4575b4", "Breganze 2" = "#74add1",
                               "Cimone" = "#d73027")) +
  labs(title = "Distribuzione % Positivi per Pool (trappola x data)",
       x = "% Positivi (su analizzati nel pool)", y = "Frequenza") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

#' Diagnostica zeri strutturali: trappole con 0% positivi in TUTTE le date
prop_zero_trap <- df %>%
  group_by(trappola, vigneto_label) %>%
  summarise(n = n(), pct_pos = sum(sfd == 1, na.rm = TRUE) / sum(!is.na(sfd)) * 100,
            .groups = "drop") %>%
  filter(pct_pos == 0)
n_zi_trap <- nrow(prop_zero_trap)
message(sprintf(
  "[ZI DIAGNOSTICA] %d trappole con 0%% positivi su tutte le date -> candidati zeri strutturali",
  n_zi_trap
))
if (n_zi_trap > 0) print(prop_zero_trap)

#' Flag trappola strutturalmente zero (usato nella zi-formula)
zi_traps <- prop_zero_trap$trappola
df_glmm <- df_glmm %>%
  mutate(zi_trap = factor(if_else(trappola_f %in% zi_traps, 1L, 0L)))
df_glmm_sex <- df_glmm_sex %>%
  mutate(zi_trap = factor(if_else(trappola_f %in% zi_traps, 1L, 0L)))
df_pool <- df_pool %>%
  mutate(zi_trap = factor(if_else(trappola %in% zi_traps, 1L, 0L)))


# ── B: Individuale con sesso (glmer) ─────────────────────────────────────────
#' Manteniamo glmm_B solo per i BLUP e la ROC individuale.
#' Per le curve fitted e il confronto principale si usa ZI-C (pool-level, glmmTMB).
glmm_B <- glmer(
  sfd ~ vigneto_f +sesso+ tempo + (1 | trappola_f),
  data = df_glmm_sex, family = binomial("logit"), control = ctrl
)
summary(glmm_B)
#"La carica batterica scende nel tempo?" Si seppure di pochissimo 
#'Senza sesso 
glmm_C <- glmer(
  sfd ~ vigneto_f + tempo + (1 | trappola_f),
  data = df_glmm_sex, family = binomial("logit"), control = ctrl
)
summary(glmm_C)

# ── ZI-C: Zero-inflato pool-level (glmmTMB) ──────────────────────────────────
zi_formula_pool <- if (n_zi_trap > 0) ~ zi_trap else ~ 0

zi_C <- tryCatch(
  glmmTMB::glmmTMB(
    cbind(n_pos, n_neg) ~ vigneto_f + tempo + (1 | trappola),
    ziformula  = zi_formula_pool,
    family     = binomial("logit"),
    data       = df_pool
  ),
  error = function(e) {
    message("[AVVISO ZI-C] glmmTMB non convergito: ", conditionMessage(e)); NULL
  }
)
if (!is.null(zi_C)) summary(zi_C)

#' Confronto AIC: GLMM B individuale vs ZI-C pool-level zero-inflato
aic_glmm <- tibble(
  Modello = c("B: individuale (glmer)",
              "ZI-C: pool-level (glmmTMB)"),
  AIC = c(
    AIC(glmm_B),
    if (!is.null(zi_C)) AIC(zi_C) else NA_real_
  ),
  BIC = c(
    BIC(glmm_B),
    if (!is.null(zi_C)) BIC(zi_C) else NA_real_
  )
) %>% arrange(AIC) %>% print()

pROC::plot.roc(pROC::roc(df_glmm_sex$sfd, fitted(glmm_C)),
               print.auc = TRUE, asp = FALSE, main = "ROC GLMM C (sesso)")
#identificazione tramite modello con effetti casuali di trappola un modello con buone performance

# CURVE FITTED – TRAIETTORIE PER SESSO E VIGNETO ==================================
tempo_seq  <- seq(min(df_glmm$tempo, na.rm = TRUE),
                  max(df_glmm$tempo, na.rm = TRUE), length.out = 200)

newdata_fx <- expand.grid(tempo = tempo_seq,
                          vigneto_f = factor(c("Breganze", "Cimone")))

newdata_fx$prob_C <- predict(glmm_C, newdata = newdata_fx, type = "response", re.form = NA)

obs_agg <- df_glmm %>%
  group_by(time, tempo, vigneto_f) %>%
  summarise(obs_pct = mean(sfd, na.rm = TRUE), n = n(), .groups = "drop")


ggplot() +
  geom_point(data = df_pool, aes(x = tempo, y = perc_pos / 100, color = vigneto_f, size = n_tot),
             alpha = 0.55) +
  geom_line(data = newdata_fx, aes(x = tempo, y = prob_C, color = vigneto_f), linewidth = 1.4) +
  scale_color_manual(values = c("Breganze" = "#2166ac", "Cimone" = "#d73027")) +
  scale_size_continuous(range = c(2, 7), name = "N pool") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title    = "GLMM_B",
       subtitle = "Una osservazione=(trappola x data)",
       x = "Giorni dall'inizio stagione", y = "P(SFD = 1)", color = "Vigneto") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "bottom")
#'Catturiamo l'effetto di trasmissione maggiore nel vigneto di cimone che piano piano va ad attenuarsi asintoticamente.


# EFFETTI CASUALI PER TRAPPOLA (BLUP) ==========================================
#' BLUP per trappola: intercette > 0 = cluster di viti malate nelle vicinanze.
#'Andiamo a valutare la conformita tra la valutazione di effetti casuali nelle varie trappole e osserviamo come coincida con il grafico delle top 10.
estrai_re <- function(modello, etichetta) {
  ranef(modello)$trappola_f %>%
    as.data.frame() %>%
    rownames_to_column("trappola") %>%
    rename(intercetta = `(Intercept)`) %>%
    mutate(modello = etichetta,
           direzione = if_else(intercetta > 0,
                               "Cluster positivo (> media)", "Sotto la media"))
}

re_B <- tryCatch(estrai_re(glmm_B, "GLMM B - con sesso"), error = function(e) NULL)

bind_rows(re_B) %>%
  left_join(df_glmm %>% distinct(trappola = as.character(trappola_f), vigneto_label),
            by = "trappola") %>%
  ggplot(aes(x = reorder(paste0(trappola, " [", modello, "]"), intercetta),
             y = intercetta, fill = direzione)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(aes(xend = paste0(trappola, " [", modello, "]"), yend = 0),
               linewidth = 0.5, alpha = 0.6) +
  geom_point(size = 3, shape = 21, color = "black") +
  facet_wrap(~ modello, scales = "free_y") +
  scale_fill_manual(values = c("Cluster positivo (> media)" = "#d73027",
                               "Sotto la media" = "#4575b4"), name = NULL) +
  coord_flip() +
  labs(title    = "Effetti Casuali per Trappola (BLUP)",
       subtitle = "Intercette > 0: probabile vicinanza a viti malate",
       x = "Trappola", y = "Intercetta casuale (log-odds)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom", strip.text = element_text(face = "bold"))

# GLMM CON RISPOSTA RICLASSIFICATA DALLA SOGLIA YOUDEN =========================
#' La variabile risposta non e' piu' sfd (soglia fissa di laboratorio) ma
#' sfd_youden_delta

df_youden <- df_qc %>%
  filter(!is.na(sfd_youden_delta), !is.na(vigneto), !is.na(trappola)) %>%
  mutate(
    vigneto_f  = factor(vigneto, levels = c(0, 1), labels = c("Breganze", "Cimone")),
    trappola_f = factor(trappola)
  )

df_youden_sex <- df_youden %>% filter(!is.na(sesso))

ctrl_y <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

#' Modello con sesso (solo campioni con sesso noto)
glmm_youden_sex <- glmer(
  sfd_youden_delta ~ vigneto_f + sesso + tempo + (1 | trappola_f),
  data    = df_youden_sex,
  family  = binomial("logit"),
  control = ctrl_y
)
summary(glmm_youden_sex)
#' Verifichiamo se la pendenza temporale del modello cambia tra luglio e settembre.
#' La pendenza non risulta significativamente diversa, indicando che la dinamica di prevalenza è sostanzialmente costante nell’arco della stagione.
#' Di conseguenza, la soglia di classificazione non necessita di essere resa dinamica nel tempo.

#' "Il decadimento temporale della carica fitoplasmatica è sufficientemente forte da far transitare i campioni dalla classe Positivo a Negativo?"
#' La risposta è negativa: la maggior parte dei campioni mostra fluttuazioni che li mantengono all’interno della stessa classe di classificazione, senza un passaggio statisticamente rilevante tra positivo e negativo.

#' Curva ROC sulla risposta riclassificata con soglia Youden (verifica della discriminazione del modello)
pROC::plot.roc(
  pROC::roc(df_youden_sex$sfd_youden_delta, fitted(glmm_youden_sex), quiet = TRUE),
  print.auc = TRUE, asp = FALSE,
  main = "ROC – GLMM Youden con sesso"
)

#' La soglia Youden ottimizza il cut-off massimizzando la somma di sensibilità e specificità sul "Ground Truth" fornito dalle etichette originali del CREA.
#' Poiché tali etichette sono state assegnate in modo statico (presenza/assenza globale senza tenere conto della degradazione temporale della carica), anche la soglia Youden eredita questa staticità.

#' È importante notare come, nonostante la non significatività delle variabili esplicative (vigneto e tempo), il modello con la soglia Youden ottenga valori di classificazione complessivamente migliori rispetto al modello basato sulle etichette originali fisse.
#' Tuttavia, questa migliore performance comporta una perdita di efficacia temporale nell’identificazione dell’andamento della carica fitoplasmatica, poiché la soglia statica tende a “congelare” le fluttuazioni osservate nel tempo.