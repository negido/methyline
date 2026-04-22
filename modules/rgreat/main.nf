process rgreat {
    tag "$sampleId"
    label 'medium'
    container 'quay.io/biocontainers/bioconductor-rgreat:2.12.2--r45ha27e39d_0'
    input:
    tuple val(sampleId), path(dmr_bed)
    val referenceGenome

    output:
    path "GREAT_enrichment.txt", emit: enrichment
    path "GREAT_summary.html", emit: summary
    path "GO_enrichment.pdf", emit: go_plot
    path "GREAT_volcano_*.pdf", emit: volcano_plots
    path "GREAT_region_gene_associations.pdf", emit: region_gene_plot
    path "GREAT_region_gene_associations.tsv", emit: region_gene_tsv

    script:
    """
    Rscript - <<'EOF'
    library(rGREAT)
    
    write_empty_pdf <- function(path, msg) {
        pdf(path)
        plot.new()
        text(0.5, 0.5, msg)
        dev.off()
    }

    # Read DMR regions
    dmr_gr <- tryCatch(
        rtracklayer::import.bed("${dmr_bed}"),
        error = function(e) GRanges()
    )
    
    # Map organism name
    organism_id <- switch("${referenceGenome}",
        "hg38" = "hg38",
        "hg19" = "hg19",
        "mm10" = "mm10",
        "hg38"
    )
    
    if (length(dmr_gr) == 0) {
        write.table(data.frame(), file = "GREAT_enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE)
        write_empty_pdf("GO_enrichment.pdf", "No regions available for GREAT")
        write_empty_pdf("GREAT_region_gene_associations.pdf", "No regions available for GREAT")
        write_empty_pdf("GREAT_volcano_GO_Biological_Process.pdf", "No regions available for GREAT")
        write_empty_pdf("GREAT_volcano_GO_Cellular_Component.pdf", "No regions available for GREAT")
        write_empty_pdf("GREAT_volcano_GO_Molecular_Function.pdf", "No regions available for GREAT")
        write.table(data.frame(), file = "GREAT_region_gene_associations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
        cat("<html><head><title>GREAT Analysis Summary</title></head><body>", file="GREAT_summary.html")
        cat("<h1>GREAT Enrichment Analysis</h1>", file="GREAT_summary.html", append=TRUE)
        cat("<p>No input regions available.</p>", file="GREAT_summary.html", append=TRUE)
        cat("</body></html>", file="GREAT_summary.html", append=TRUE)
        quit(save = "no", status = 0)
    }

    # Perform GREAT analysis
    job <- tryCatch(
        submitGreatJob(
            gr = dmr_gr,
            genome = organism_id,
            request_interval = 5,
            help = FALSE
        ),
        error = function(e) {
            message(paste0(">> GREAT submit failed: ", conditionMessage(e)))
            NULL
        }
    )

    enrichment_list <- list()
    selected_ontologies <- character(0)
    if (!is.null(job)) {
        # Keep only human-focused ontologies + GO and avoid mouse phenotypes.
        available_onts <- tryCatch(availableOntologies(job), error = function(e) character(0))
        preferred_onts <- c(
            "GO Biological Process",
            "GO Cellular Component",
            "GO Molecular Function",
            "Human Phenotype"
        )
        selected_ontologies <- intersect(preferred_onts, available_onts)

        if (length(selected_ontologies) > 0) {
            enrichment_list <- tryCatch(
                getEnrichmentTables(job, ontology = selected_ontologies, download_by = "tsv"),
                error = function(e) {
                    message(paste0(">> TSV download failed (", conditionMessage(e), "); retrying with JSON."))
                    tryCatch(getEnrichmentTables(job, ontology = selected_ontologies), error = function(e2) list())
                }
            )
        } else {
            message(">> No GO/Human ontologies available for this GREAT job")
        }
    }

    combined_list <- lapply(names(enrichment_list), function(ontology_name) {
        tbl <- enrichment_list[[ontology_name]]
        if (is.null(tbl) || nrow(tbl) == 0) {
            return(NULL)
        }

        df <- as.data.frame(tbl, stringsAsFactors = FALSE)
        if (!"name" %in% colnames(df) && !is.null(rownames(tbl)) && all(rownames(tbl) != "")) {
            df[["name"]] <- rownames(tbl)
        }
        if (!"ID" %in% colnames(df)) {
            df[["ID"]] <- NA_character_
        }
        df[["ontology"]] <- ontology_name
        df
    })

    combined_results <- do.call(rbind, Filter(Negate(is.null), combined_list))

    if (is.null(combined_results) || nrow(combined_results) == 0) {
        write.table(data.frame(), file = "GREAT_enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        if ("Hyper_Raw_PValue" %in% colnames(combined_results)) {
            hp <- suppressWarnings(as.numeric(combined_results[["Hyper_Raw_PValue"]]))
            combined_results <- combined_results[order(hp, na.last = TRUE), , drop = FALSE]
        }
        write.table(combined_results, file = "GREAT_enrichment.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    # Generate plots
    pdf("GO_enrichment.pdf")
    
    plot_df <- NULL
    if (!is.null(combined_results) && nrow(combined_results) > 0) {
        p_col <- if ("Hyper_Raw_PValue" %in% colnames(combined_results)) "Hyper_Raw_PValue" else if ("Binom_Raw_PValue" %in% colnames(combined_results)) "Binom_Raw_PValue" else NULL
        if (!is.null(p_col)) {
            pvals <- suppressWarnings(as.numeric(combined_results[[p_col]]))
            ok <- which(is.finite(pvals) & pvals > 0)
            if (length(ok) > 0) {
                ranked <- ok[order(pvals[ok])]
                top_idx <- head(ranked, 15)
                term_labels <- as.character(combined_results[["name"]][top_idx])
                term_labels[is.na(term_labels) | term_labels == ""] <- as.character(combined_results[["ID"]][top_idx])
                term_labels <- paste0(term_labels, " [", combined_results[["ontology"]][top_idx], "]")
                plot_df <- data.frame(term = term_labels, p_value = -log10(pvals[top_idx]), stringsAsFactors = FALSE)
            }
        }
    }

    if (is.null(plot_df) || nrow(plot_df) == 0) {
        plot.new()
        text(0.5, 0.5, "No enrichment terms to plot")
    } else {
        ordp <- order(plot_df\$p_value)
        barplot(
            plot_df\$p_value[ordp],
            names.arg = plot_df\$term[ordp],
            beside = FALSE,
            main = "Top Enriched Terms in DMRs",
            xlab = "-log10(p-value)",
            las = 2,
            cex.names = 0.7
        )
    }
    dev.off()
    
    # Volcano plots per GO ontology.
    go_onts <- intersect(
        c("GO Biological Process", "GO Cellular Component", "GO Molecular Function"),
        names(enrichment_list)
    )
    for (ont in go_onts) {
        safe_name <- gsub("[^A-Za-z0-9]+", "_", ont)
        out_pdf <- paste0("GREAT_volcano_", safe_name, ".pdf")
        pdf(out_pdf)
        ok <- tryCatch({
            plotVolcano(job, ontology = ont)
            TRUE
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, paste0("Could not generate volcano for ", ont))
            FALSE
        })
        dev.off()
    }
    if (length(go_onts) == 0) {
        write_empty_pdf("GREAT_volcano_none.pdf", "No GO ontologies available for volcano plots")
    }

    # Region-gene association plot.
    pdf("GREAT_region_gene_associations.pdf")
    tryCatch(
        plotRegionGeneAssociations(job),
        error = function(e) {
            plot.new()
            text(0.5, 0.5, "Could not generate region-gene association plot")
        }
    )
    dev.off()

    # Region-gene association table.
    assoc_df <- tryCatch({
        assoc <- getRegionGeneAssociations(job)
        data.frame(
            chr = as.character(GenomicRanges::seqnames(assoc)),
            start = GenomicRanges::start(assoc),
            end = GenomicRanges::end(assoc),
            annotated_genes = sapply(assoc[["annotated_genes"]], function(x) paste(as.character(x), collapse = ";")),
            dist_to_TSS = sapply(assoc[["dist_to_TSS"]], function(x) paste(as.character(x), collapse = ";")),
            stringsAsFactors = FALSE
        )
    }, error = function(e) data.frame())
    write.table(assoc_df, file = "GREAT_region_gene_associations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Generate HTML summary
    cat("<html><head><title>GREAT Analysis Summary</title></head><body>", file="GREAT_summary.html")
    cat("<h1>GREAT Enrichment Analysis</h1>", file="GREAT_summary.html", append=TRUE)
    cat(paste("<p>Total DMRs analyzed: ", length(dmr_gr), "</p>"), 
        file="GREAT_summary.html", append=TRUE)
    cat(paste("<p>Taxon: ${referenceGenome}</p>"), file="GREAT_summary.html", append=TRUE)
    n_categories <- if (!is.null(combined_results) && nrow(combined_results) > 0) length(unique(combined_results[["ontology"]])) else 0
    cat(paste("<p>Categories processed: ", n_categories, "</p>"), file="GREAT_summary.html", append=TRUE)
    cat(paste("<p>Ontologies selected: ", ifelse(length(selected_ontologies) > 0, paste(selected_ontologies, collapse = ", "), "none"), "</p>"), file="GREAT_summary.html", append=TRUE)
    cat("</body></html>", file="GREAT_summary.html", append=TRUE)
    
    cat("GREAT analysis completed\\n")
    EOF
    """

}
