#!/usr/bin/env Rscript

library(qqman)
library(stringr)
library(optparse)
library(data.table)
library(MatrixEQTL)


parse_arguments <- function() {
    parser <- OptionParser(description = "A warpper of MatrixEQTL.")
    parser <- add_option(parser, opt_str = c("-o", "--save-pref"), dest = "save_pref", help = "The prefix for the output.")
    parser <- add_option(parser, opt_str = c("-g", "--gntp-file"), dest = "gntp_file", help = "The file from which read the genotypes.")
    parser <- add_option(parser, opt_str = c("-i", "--gntp-info-file"), dest = "gntp_info_file", help = "The file from which read the genotype information.")
    parser <- add_option(parser, opt_str = c("-p", "--phtp-file"), dest = "phtp_file", help = "The file from which read the phenotypes.")
    parser <- add_option(parser, opt_str = c("-c", "--cvrt-file"), dest = "cvrt_file", help = "The file from which read the coveriates.")

    opts_args <- parse_args2(parser)
    return(opts_args)
}


make_sliced_data <- function(fp, foc = "NA", fss = 2000) {
    # foc: fileOmitCharacters; default: NA
    # fss: fileSliceSize; default: 2000
    scdt <- SlicedData$new()
    scdt$fileOmitCharacters <- foc
    scdt$fileSliceSize <- fss

    if (is.character(fp)) {
        scdt$fileSkipRows <- 1
        scdt$fileSkipColumns <- 1
        scdt$LoadFile(fp)
    } else if (is.matrix(fp))
        scdt$CreateFromMatrix(fp)
    else
        stop("Unsupported type for creating a SliceData object.")

    return(scdt)
}


load_dtfm <- function(phtp_file, gntp_file, cvrt_file, ...) {
    phtp_dtfm <- read.table(phtp_file, header = TRUE, sep = "\t", row.names = 1)
    gntp_dtfm <- read.table(gntp_file, header = TRUE, sep = "\t", row.names = 1)
    cvrt_dtfm <- NULL

    candidate_samples <- intersect(colnames(gntp_dtfm), colnames(phtp_dtfm))
    if (!is.null(cvrt_file)) {
        cvrt_dtfm <- read.table(cvrt_file, header = TRUE, sep = "\t", row.names = 1)
        candidate_samples <- intersect(candidate_samples, colnames(cvrt_dtfm))
        cvrt_dtfm <- cvrt_dtfm[, candidate_samples]
    }

    phtp_dtfm <- phtp_dtfm[, candidate_samples]
    gntp_dtfm <- gntp_dtfm[, candidate_samples]

    dtfm_list <- list(phtp_dtfm = phtp_dtfm, gntp_dtfm = gntp_dtfm, cvrt_dtfm = cvrt_dtfm)

    return(dtfm_list)
}


map_qtls <- function(dtfm_list, save_pref = "./", pv_opt_thrd = 1.0, use_model = MatrixEQTL::modelLINEAR, err_cov = numeric()) {
    phtp_mtrx <- as.matrix(dtfm_list$phtp_dtfm)
    gene <- make_sliced_data(phtp_mtrx)

    gntp_mtrx <- as.matrix(dtfm_list$gntp_dtfm)
    snps <- make_sliced_data(gntp_mtrx)

    cvrt <- NULL
    cvrt_mtrx <- as.matrix(dtfm_list$cvrt_dtfm)
    if (!is.null(cvrt_mtrx))
        cvrt <- make_sliced_data(cvrt_mtrx)

    output_file_name <- file.path(save_pref, "qtls.tsv")
    me <- MatrixEQTL::Matrix_eQTL_engine(
        snps = snps,
        gene = gene,
        cvrt = cvrt,
        output_file_name = output_file_name,
        pvOutputThreshold = pv_opt_thrd,
        useModel = use_model,
        errorCovariance = err_cov,
        verbose = FALSE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = TRUE
    )

    return(me)
}


lambda <- function(stvc, sttp="PVAL", rnd=3) {
    if(sttp == "Z") {
        z = stvc
    } else if (sttp == "CHISQ") {
        z = sqrt(stvc)
    } else if (sttp == "PVAL") {
        z = qnorm(stvc / 2)
    }

    return(round(median(z^2) / qchisq(0.5, 1), rnd))
}

draw_qq_man <- function(meobj, gntp_info_file, info_col_vec, mhtn_fig_p_thrd = 0.01, pv_opt_thrd = 1.0, save_pref = "./") {
    snpid_idx <- info_col_vec[1]
    chrom_idx <- info_col_vec[2]
    position_idx <- info_col_vec[3]

    qtls <- meobj$all$eqtls

    gntp_info_dtfm <- fread(gntp_info_file, data.table = FALSE)
    qtls_info_dtfm <- merge(qtls, gntp_info_dtfm, by.x = "snps", by.y = snpid_idx)

    for (phtp in unique(qtls_info_dtfm$gene)) {
        phtp_qtls <- qtls_info_dtfm[qtls_info_dtfm$gene == phtp, ]

        phtp_qtls_name <- file.path(save_pref, "per_trait", str_glue("{phtp}.tsv.gz"))
        fwrite(phtp_qtls, file = phtp_qtls_name)

        phtp_qtls <- as.data.frame(phtp_qtls[, c("snps", chrom_idx, position_idx, "pvalue")])

        if (nrow(phtp_qtls) > 0) {
            colnames(phtp_qtls) <- c("SNP", "CHR", "BP", "P")

            if (pv_opt_thrd < 1.0) {
                lmb <- "NULL"
            } else {
                lmb <- lambda(phtp_qtls$P, "PVAL")
            }

            cat(str_glue("Inflation factor for {phtp}: {lmb}"), "\n")
            phtp_qtls_mhtn <- phtp_qtls[phtp_qtls$P <= mhtn_fig_p_thrd, ]

            pdf(file.path(save_pref, "manhattan", str_glue("manhattan_plot.{phtp}.pdf")), width = 16, height = 9)
            manhattan(phtp_qtls_mhtn, main = str_glue("Manhattan plot for {phtp}"), ylab = "p-value(-log10)", annotateTop = TRUE)
            dev.off()

            pdf(file.path(save_pref, "qq", str_glue("qq_plot.{phtp}.pdf")), width = 16, height = 16)
            qq(phtp_qtls$P, main = str_glue("Q-Q plot {phtp} (lambda={lmb})"), pch = 18, cex = 1, las = 1)
            dev.off()
        } else {
            warning(str_glue("No QTL is found under p value ({pv_opt_thrd}) for {phtp} (trait or gene)"))
        }

    }
}


main <- function() {
    opts <- parse_arguments()$options

    save_pref <- opts$save_pref
    phtp_file <- opts$phtp_file
    cvrt_file <- opts$cvrt_file
    gntp_file <- opts$gntp_file
    gntp_info_file <- opts$gntp_info_file

    if (is.null(save_pref))
        save_pref <- "./"

    if (is.null(phtp_file))
        stop("-p/--phtp-file is required!!!")

    if (is.null(gntp_file))
        stop("-g/--gntp-file is required!!!")

    if (is.null(gntp_info_file))
        stop("-i/--gntp-info-file is required!!!")

    dtfms <- load_dtfm(phtp_file, gntp_file, cvrt_file)
    meobj <- map_qtls(dtfms, save_pref = save_pref)

    info_col_vec <- c("rsID", "SequenceName", "Position", "EffectAllele", "AlternativeAllele")
    draw_qq_man(meobj, gntp_info_file, info_col_vec, save_pref=save_pref)
}


main()
