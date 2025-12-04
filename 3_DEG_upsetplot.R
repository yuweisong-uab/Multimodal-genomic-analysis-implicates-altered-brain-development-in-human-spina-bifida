library(dplyr)
library(purrr)
library(UpSetR)
library(grid)

deg_dir <- "output/summary/DEG"

celltype_levels <- c(
  "Astrocytes", "OPC", "Oligodendrocytes",
  "Inhibitory neurons", "Excitatory neurons", "Microglia"
)

## cell type colors
my_cols <- c(
  'Astrocytes'          = '#8F7700',
  'OPC'                 = '#3B3B3B',
  'Oligodendrocytes'    = '#EFC001',
  'Microglia'           = '#CD534C',
  'Excitatory neurons'  = '#59753E',
  'Inhibitory neurons'  = '#394E2F'
)

## file types to process
file_types <- c(
  "GO_BP_up_sig",
  "GO_BP_down_sig",
  "KEGG_up_sig",
  "KEGG_down_sig"
)

## pair definitions (orig.ident)
pair_levels <- list(
  Pair_HCPUAB014 = c("HCPUAB014", "control_y10"),
  Pair_HCPUAB023 = c("HCPUAB023", "control_m1")
)

## ============================================================
## LOOP over each pair, then each file type
## ============================================================
for (pair_name in names(pair_levels)) {
  
  message("###############################")
  message("Processing pair: ", pair_name)
  message("###############################")
  
  for (ftype in file_types) {
    
    message("==== Processing file type: ", ftype, " ====")
    
    set_list <- list()   # terms per cell type for this pair + file type
    
    # -------------------------------
    # LOAD FILES FOR EACH CELL TYPE
    # -------------------------------
    for (ct in celltype_levels) {
      base_name <- paste0(pair_name, "_", gsub(" ", "_", ct))
      file_name <- file.path(deg_dir, paste0(base_name, "_", ftype, ".csv"))
      
      if (!file.exists(file_name)) {
        message("[", pair_name, " | ", ftype, "] ", ct, " file not found -> skip")
        next
      }
      
      df <- tryCatch(
        read.csv(file_name, stringsAsFactors = FALSE),
        error = function(e) {
          message("Error reading ", file_name)
          return(NULL)
        }
      )
      if (is.null(df) || nrow(df) == 0) {
        message("[", pair_name, " | ", ftype, "] ", ct, " is empty -> skip")
        next
      }
      
      # Need a Description column
      if (!"Description" %in% colnames(df)) {
        message("[", pair_name, " | ", ftype, "] ", ct, " has no Description col -> skip")
        next
      }
      
      terms <- unique(df$Description)
      if (length(terms) > 0) {
        set_list[[ct]] <- terms
      } else {
        message("[", pair_name, " | ", ftype, "] ", ct, " has 0 terms -> skip")
      }
    } # end for ct
    
    # Must have ≥2 sets to make intersections / UpSet
    if (length(set_list) < 2) {
      message("[", pair_name, " | ", ftype, "] Not enough cell types -> skip")
      next
    }
    
    # -------------------------------
    # BUILD UNIQUE + INTERSECTIONS
    # -------------------------------
    cts <- names(set_list)
    final_results <- list()
    
    ## 1. UNIQUE TERMS per cell type
    for (ct in cts) {
      other_cts  <- setdiff(cts, ct)
      uniq_terms <- setdiff(set_list[[ct]], unlist(set_list[other_cts]))
      
      if (length(uniq_terms) > 0) {
        final_results[[paste0("unique__", ct)]] <-
          data.frame(
            Set_Type   = "unique",
            CellTypes  = ct,
            Description = uniq_terms,
            stringsAsFactors = FALSE
          )
      }
    }
    
    ## 2. ALL intersections of size >= 2
    if (length(cts) >= 2) {
      for (k in 2:length(cts)) {
        combs <- combn(cts, k, simplify = FALSE)
        for (comb in combs) {
          inter_terms <- reduce(set_list[comb], intersect)
          if (length(inter_terms) > 0) {
            final_results[[paste0("intersect__", paste(comb, collapse = "__"))]] <-
              data.frame(
                Set_Type   = "intersect",
                CellTypes  = paste(comb, collapse = "__"),
                Description = inter_terms,
                stringsAsFactors = FALSE
              )
          }
        }
      }
    }
    
    # combine all results (if any) and save CSV
    if (length(final_results) == 0) {
      message("[", pair_name, " | ", ftype, "] No unique/intersection terms -> skip summary CSV")
    } else {
      final_df <- bind_rows(final_results)
      
      out_summary_file <- file.path(
        deg_dir,
        paste0(pair_name, "_", ftype, "_set_summary.csv")
      )
      write.csv(final_df, out_summary_file, row.names = FALSE)
      message("Saved summary: ", out_summary_file)
    }
    
    # -------------------------------
    # MAKE UPSET PLOT (top 20 intersections)
    # -------------------------------
    upset_input <- UpSetR::fromList(set_list)
    
    ## y-axis label: GO vs KEGG
    if (grepl("^GO", ftype)) {
      y_label <- "Number of GO BP terms"
    } else {
      y_label <- "Number of KEGG pathways"
    }
    
    ## patient/control labels by pair
    if (pair_name == "Pair_HCPUAB014") {
      patient_lab <- "10Y patient"
      control_lab <- "10Y control"
    } else if (pair_name == "Pair_HCPUAB023") {
      patient_lab <- "1M patient"
      control_lab <- "1M control"
    } else {
      patient_lab <- "patient"
      control_lab <- "control"
    }
    
    ## Title text by file type
    if (ftype == "GO_BP_up_sig") {
      plot_title <- paste0(
        "Shared and cell-type–specific GO BP terms enriched\n",
        "in genes increased in ", patient_lab, " versus ", control_lab
      )
    } else if (ftype == "GO_BP_down_sig") {
      plot_title <- paste0(
        "Shared and cell-type–specific GO BP terms enriched\n",
        "in genes decreased in ", patient_lab, " versus ", control_lab
      )
    } else if (ftype == "KEGG_up_sig") {
      plot_title <- paste0(
        "Shared and cell-type–specific KEGG pathways enriched\n",
        "in genes increased in ", patient_lab, " versus ", control_lab
      )
    } else if (ftype == "KEGG_down_sig") {
      plot_title <- paste0(
        "Shared and cell-type–specific KEGG pathways enriched\n",
        "in genes decreased in ", patient_lab, " versus ", control_lab
      )
    } else {
      plot_title <- paste0(
        "Shared and cell-type–specific terms\n",
        "in ", patient_lab, " versus ", control_lab
      )
    }
    
    ## color vector for the sets actually present in this plot
    sets_to_plot <- names(set_list)
    set_bar_cols <- my_cols[sets_to_plot]
    set_bar_cols[is.na(set_bar_cols)] <- "grey50"  # fallback if any unexpected name
    
    out_upset_pdf <- file.path(
      deg_dir,
      paste0(pair_name, "_", ftype, "_UpSet.pdf")
    )
    
    pdf(out_upset_pdf, width = 9, height = 4)
    print(
      UpSetR::upset(
        upset_input,
        sets            = sets_to_plot,
        nsets           = length(set_list),
        nintersects     = 20,            # top 20 intersections
        order.by        = "freq",
        sets.x.label    = "Term count",
        mainbar.y.label = y_label,
        sets.bar.color  = set_bar_cols,  # <-- color each cell type bar
        matrix.color    = "black",
        main.bar.color  = "black",
        text.scale      = c(1.2, 1.2, 1, 1, 1.2, 1)
      )
    )
    grid::grid.text(
      plot_title,
      x  = 0.65,
      y  = 0.85,
      gp = grid::gpar(fontsize = 11, fontface = "bold")
    )
    dev.off()
    
    message("Saved UpSet: ", out_upset_pdf)
    
  } # end for ftype
  
} # end for pair_name