function (txt_folder, yaml_obj = list(), report_filenames = NULL) 
{
  DEBUG_PTXQC = FALSE
  time_start = Sys.time()
  if (!any(file.info(txt_folder)$isdir, na.rm = TRUE)) {
    stop(paste0("Argument txt_folder with value '", txt_folder, 
                "' is not a valid directory\n")) #This prompt will show if you pass a file that is not maxquont
  }
  txt_files = list()
  txt_files$param = "parameters.txt"
  txt_files$summary = "summary.txt"
  txt_files$groups = "proteinGroups.txt"
  txt_files$evd = "evidence.txt"
  txt_files$msms = "msms.txt"
  txt_files$msmsScan = "msmsScans.txt"
  txt_files$mqpar = "mqpar.xml"
  txt_files = lapply(txt_files, function(file) file.path(txt_folder, 
                                                         file))
  if (class(yaml_obj) != "list") {
    stop(paste0("Argument 'yaml_obj' is not of type list\n"))
  }
  yc = YAMLClass$new(yaml_obj)
  if (is.null(report_filenames)) {
    use_extended_reportname = yc$getYAML("PTXQC$ReportFilename$extended", 
                                         TRUE)
    rprt_fns = getReportFilenames(txt_folder, use_extended_reportname)
  }
  else {
    rprt_fns = report_filenames
  }
  {
    param_name_PTXQC_UseLocalMQPar = "PTXQC$UseLocalMQPar"
    param_def_PTXQC_UseLocalMQPar = TRUE
    param_useMQPAR = yc$getYAML(param_name_PTXQC_UseLocalMQPar, 
                                param_def_PTXQC_UseLocalMQPar)
    enabled_parameters = yc$getYAML("File$Parameters$enabled", 
                                    TRUE) & file.exists(txt_files$param)
    add_fs_col = yc$getYAML("PTXQC$NameLengthMax_num", 10)
    enabled_summary = yc$getYAML("File$Summary$enabled", 
                                 TRUE) & file.exists(txt_files$summary)
    id_rate_bad = yc$getYAML("File$Summary$IDRate$Thresh_bad_num", 
                             20)
    id_rate_great = yc$getYAML("File$Summary$IDRate$Thresh_great_num", 
                               35)
    GL_name_min_length = 8
    enabled_proteingroups = yc$getYAML("File$ProteinGroups$enabled", 
                                       TRUE) & file.exists(txt_files$groups)
    enabled_pg_ratioLabIncThresh = yc$getYAML("File$ProteinGroups$RatioPlot$LabelIncThresh_num", 
                                              4)
    param_name_PG_intThresh = "File$ProteinGroups$IntensityThreshLog2_num"
    param_def_PG_intThresh = 25
    param_PG_intThresh = yc$getYAML(param_name_PG_intThresh, 
                                    param_def_PG_intThresh)
    if (!is.numeric(param_PG_intThresh) || !(param_PG_intThresh %in% 
                                             1:100)) {
      cat("YAML value for '" %+% param_name_PG_intThresh %+% 
            "' is invalid ('" %+% param_PG_intThresh %+% 
            "'). Using default of " %+% param_def_PG_intThresh %+% 
            ".")
      param_PG_intThresh = param_def_PG_intThresh
    }
    enabled_evidence = yc$getYAML("File$Evidence$enabled", 
                                  TRUE) & file.exists(txt_files$evd)
    param_name_EV_protThresh = "File$Evidence$ProteinCountThresh_num"
    param_def_EV_protThresh = 3500
    param_EV_protThresh = yc$getYAML(param_name_EV_protThresh, 
                                     param_def_EV_protThresh)
    if (!is.numeric(param_EV_protThresh) || !(param_EV_protThresh %in% 
                                              1:1e+05)) {
      cat("YAML value for '" %+% param_name_EV_protThresh %+% 
            "' is invalid ('" %+% param_EV_protThresh %+% 
            "'). Using default of " %+% param_def_EV_protThresh %+% 
            ".")
      param_EV_protThresh = param_def_EV_protThresh
    }
    param_name_EV_intThresh = "File$Evidence$IntensityThreshLog2_num"
    param_def_EV_intThresh = 23
    param_EV_intThresh = yc$getYAML(param_name_EV_intThresh, 
                                    param_def_EV_intThresh)
    if (!is.numeric(param_EV_intThresh) || !(param_EV_intThresh %in% 
                                             1:100)) {
      cat("YAML value for '" %+% param_name_EV_intThresh %+% 
            "' is invalid ('" %+% param_EV_intThresh %+% 
            "'). Using default of " %+% param_def_EV_intThresh %+% 
            ".")
      param_EV_intThresh = param_def_EV_intThresh
    }
    param_name_EV_pepThresh = "File$Evidence$PeptideCountThresh_num"
    param_def_EV_pepThresh = 15000
    param_EV_pepThresh = yc$getYAML(param_name_EV_pepThresh, 
                                    param_def_EV_pepThresh)
    if (!is.numeric(param_EV_pepThresh) || !(param_EV_pepThresh %in% 
                                             1:1e+06)) {
      cat("YAML value for '" %+% param_name_EV_pepThresh %+% 
            "' is invalid ('" %+% param_EV_pepThresh %+% 
            "'). Using default of " %+% param_def_EV_pepThresh %+% 
            ".")
      param_EV_pepThresh = param_def_EV_pepThresh
    }
    contaminant_default = list(cont_MYCO = c(name = "MYCOPLASMA", 
                                             threshold = 1))
    yaml_contaminants = yc$getYAML("File$Evidence$SpecialContaminants", 
                                   contaminant_default)
    param_name_EV_MatchingTolerance = "File$Evidence$MQpar_MatchingTimeWindow_num"
    param_def_EV_MatchingTolerance = 1
    param_EV_MatchingTolerance = yc$getYAML(param_name_EV_MatchingTolerance, 
                                            param_def_EV_MatchingTolerance)
    if (param_useMQPAR) {
      v = getMQPARValue(txt_files$mqpar, "matchingTimeWindow")
      if (!is.null(v)) {
        param_EV_MatchingTolerance = yc$setYAML(param_name_EV_MatchingTolerance, 
                                                as.numeric(v))
      }
    }
    param_name_mbr = "File$Evidence$MatchBetweenRuns_wA"
    param_evd_mbr = yc$getYAML(param_name_mbr, "auto")
    param_name_EV_PrecursorTolPPM = "File$Evidence$MQpar_firstSearchTol_num"
    param_def_EV_PrecursorTolPPM = 20
    param_EV_PrecursorTolPPM = yc$getYAML(param_name_EV_PrecursorTolPPM, 
                                          param_def_EV_PrecursorTolPPM)
    if (param_useMQPAR) {
      v = getMQPARValue(txt_files$mqpar, "firstSearchTol")
      if (!is.null(v)) {
        param_EV_PrecursorTolPPM = yc$setYAML(param_name_EV_PrecursorTolPPM, 
                                              as.numeric(v))
      }
    }
    param_name_EV_PrecursorOutOfCalSD = "File$Evidence$firstSearch_outOfCalWarnSD_num"
    param_def_EV_PrecursorOutOfCalSD = 2
    param_EV_PrecursorOutOfCalSD = yc$getYAML(param_name_EV_PrecursorOutOfCalSD, 
                                              param_def_EV_PrecursorOutOfCalSD)
    param_name_EV_PrecursorTolPPMmainSearch = "File$Evidence$MQpar_mainSearchTol_num"
    param_def_EV_PrecursorTolPPMmainSearch = NA
    param_EV_PrecursorTolPPMmainSearch = yc$getYAML(param_name_EV_PrecursorTolPPMmainSearch, 
                                                    param_def_EV_PrecursorTolPPMmainSearch)
    if (param_useMQPAR) {
      v = getMQPARValue(txt_files$mqpar, "mainSearchTol")
      if (!is.null(v)) {
        param_EV_PrecursorTolPPMmainSearch = yc$setYAML(param_name_EV_PrecursorTolPPMmainSearch, 
                                                        as.numeric(v))
      }
    }
    if (is.na(param_EV_PrecursorTolPPMmainSearch)) {
      warning("PTXQC: Cannot draw borders for calibrated mass error, since neither 'File$Evidence$MQpar_mainSearchTol_num' is set nor a mqpar.xml file is present!", 
              immediate. = TRUE)
    }
    enabled_msms = yc$getYAML("File$MsMs$enabled", TRUE) & 
      file.exists(txt_files$msms)
    enabled_msmsscans = yc$getYAML("File$MsMsScans$enabled", 
                                   TRUE) & file.exists(txt_files$msmsScan)
    param_name_MSMSScans_ionInjThresh = "File$MsMsScans$IonInjectionThresh_num"
    param_def_MSMSScans_ionInjThresh = 10
    param_MSMSScans_ionInjThresh = yc$getYAML(param_name_MSMSScans_ionInjThresh, 
                                              param_def_MSMSScans_ionInjThresh)
    if (!is.numeric(param_MSMSScans_ionInjThresh)) {
      cat("YAML value for '" %+% param_name_MSMSScans_ionInjThresh %+% 
            "' is invalid ('" %+% param_MSMSScans_ionInjThresh %+% 
            "'). Using default of " %+% param_def_MSMSScans_ionInjThresh %+% 
            ".")
      param_MSMSScans_ionInjThresh = param_def_MSMSScans_ionInjThresh
    }
    param_name_PTXQC_OutputFormats = "PTXQC$OutputFormats"
    out_formats_supported = c("html", "plainPDF")
    param_def_PTXQC_OutputFormats = out_formats_supported
    param_OutputFormats = yc$getYAML(param_name_PTXQC_OutputFormats, 
                                     param_def_PTXQC_OutputFormats)
    param_name_PTXQC_PageNumbers = "PTXQC$PlainPDF$AddPageNumbers"
    param_def_PTXQC_PageNumbers = "on"
    param_PageNumbers = yc$getYAML(param_name_PTXQC_PageNumbers, 
                                   param_def_PTXQC_PageNumbers)
  }
  mq = MQDataReader$new()
  mq$readMappingFile(rprt_fns$filename_sorting)
  lst_qcMetrics = getMetricsObjects(DEBUG_PTXQC)
  df.meta = getMetaData(lst_qcMetrics = lst_qcMetrics)
  df.meta
  lst_qcMetrics_ord = lst_qcMetrics[df.meta$.id]
  i = 1
  for (i in 1:nrow(df.meta)) {
    pname = paste0("order$", df.meta$.id[i])
    pval = df.meta$order[i]
    param = yc$getYAML(pname, pval)
    if (is.numeric(param)) {
      lst_qcMetrics_ord[[i]]$orderNr = param
    }
    else {
      stop("YAML param '", pname, "' is not numeric (", 
           param, "). Please fix the YAML configuration!")
    }
  }
  df.meta = getMetaData(lst_qcMetrics = lst_qcMetrics)
  lst_qcMetrics_ord = lst_qcMetrics[df.meta$.id]
  yc$writeYAML(rprt_fns$yaml_file)
  if (enabled_parameters) {
    d_parAll = mq$readMQ(txt_files$param, type = "par")
    lst_qcMetrics[["qcMetric_PAR"]]$setData(d_parAll)
  }
  if (enabled_summary) {
    d_smy = mq$readMQ(txt_files$summary, type = "sm", add_fs_col = add_fs_col)
    lst_qcMetrics[["qcMetric_SM_MSMSIdRate"]]$setData(d_smy$raw, 
                                                      id_rate_bad, id_rate_great)
  }
  if (enabled_proteingroups) {
    df_pg = mq$readMQ(txt_files$groups, type = "pg", col_subset = NA, 
                      filter = "R")
    clusterCols = list()
    colsSIL = grepv("^intensity\\.[hlm](\\.|$)", colnames(df_pg))
    colsLF = grepv("^intensity\\..", colnames(df_pg))
    colsOneCond = "intensity"
    if (length(colsSIL)) {
      plain_channel = grepv("^intensity\\.[hlm]$", colnames(df_pg))
      if (all(plain_channel == colsSIL)) 
        colsW = colsSIL
      else colsW = setdiff(colsSIL, plain_channel)
    }
    else if (length(colsLF)) {
      colsW = colsLF
    }
    else {
      colsW = colsOneCond
    }
    MAP_pg_groups = data.frame(long = colsW)
    MAP_pg_groups$short = shortenStrings(simplifyNames(delLCP(MAP_pg_groups$long, 
                                                              min_out_length = GL_name_min_length, add_dots = TRUE), 
                                                       min_out_length = GL_name_min_length))
    lst_qcMetrics[["qcMetric_PG_Cont"]]$setData(df_pg, colsW, 
                                                MAP_pg_groups)
    clusterCols$raw.intensity = colsW
    lst_qcMetrics[["qcMetric_PG_RawInt"]]$setData(df_pg, 
                                                  colsW, MAP_pg_groups, param_PG_intThresh)
    colsSIL = grepv("^lfq.intensity\\.[hlm](\\.|$)", colnames(df_pg))
    colsLF = grepv("^lfq.intensity\\..", colnames(df_pg))
    MAP_pg_groups_LFQ = NA
    if (length(c(colsSIL, colsLF)) > 0) {
      if (length(colsSIL)) {
        colsW = colsSIL
      }
      else colsW = colsLF
      MAP_pg_groups_LFQ = data.frame(long = colsW)
      MAP_pg_groups_LFQ$short = shortenStrings(simplifyNames(delLCP(MAP_pg_groups_LFQ$long, 
                                                                    min_out_length = GL_name_min_length, add_dots = TRUE), 
                                                             min_out_length = GL_name_min_length))
      clusterCols$lfq.intensity = colsW
      lst_qcMetrics[["qcMetric_PG_LFQInt"]]$setData(df_pg, 
                                                    colsW, MAP_pg_groups_LFQ, param_PG_intThresh)
    }
    colsITRAQ = grepv("^reporter.intensity.[0-9]", colnames(df_pg))
    MAP_pg_groups_ITRAQ = NA
    if (length(colsITRAQ) > 0) {
      MAP_pg_groups_ITRAQ = data.frame(long = c(colsITRAQ))
      MAP_pg_groups_ITRAQ$short = shortenStrings(simplifyNames(delLCP(MAP_pg_groups_ITRAQ$long, 
                                                                      min_out_length = GL_name_min_length, add_dots = TRUE), 
                                                               min_out_length = GL_name_min_length))
      clusterCols$reporter.intensity = colsITRAQ
      lst_qcMetrics[["qcMetric_PG_ReporterInt"]]$setData(df_pg, 
                                                         colsITRAQ, MAP_pg_groups_ITRAQ, param_PG_intThresh)
    }
    MAP_pg_groups_ALL = rbind(MAP_pg_groups, MAP_pg_groups_LFQ, 
                              MAP_pg_groups_ITRAQ)
    lst_qcMetrics[["qcMetric_PG_PCA"]]$setData(df_pg, clusterCols, 
                                               MAP_pg_groups_ALL)
    ratio_cols = grepv("^ratio\\.[hm]\\.l", colnames(df_pg))
    ratio_cols = grepv("^ratio.[hm].l.normalized", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.count", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.variability", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.significance.a", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.significance.b", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.iso.count", ratio_cols, 
                       invert = TRUE)
    ratio_cols = grepv("^ratio.[hm].l.type", ratio_cols, 
                       invert = TRUE)
    ratio_cols
    if (length(ratio_cols) > 0) {
      lst_qcMetrics[["qcMetric_PG_Ratio"]]$setData(df_pg = df_pg, 
                                                   ratio_cols = ratio_cols, thresh_LabelIncorp = enabled_pg_ratioLabIncThresh, 
                                                   GL_name_min_length = GL_name_min_length)
    }
  }
  if (enabled_evidence) {
    df_evd = mq$readMQ(txt_files$evd, type = "ev", filter = "R", 
                       col_subset = c("proteins", numeric = "Retention.Length", 
                                      numeric = "retention.time.calibration", numeric = "Retention.time$", 
                                      numeric = "Match.Time.Difference", numeric = "^intensity$", 
                                      "^Type$", numeric = "Mass\\.Error", numeric = "^uncalibrated...calibrated.", 
                                      numeric = "^m.z$", numeric = "^score$", numeric = "^fraction$", 
                                      "Raw.file", "^Protein.Group.IDs$", "Contaminant", 
                                      numeric = "[RK]\\.Count", numeric = "^Charge$", 
                                      "modified.sequence", numeric = "^Mass$", "^protein.names$", 
                                      numeric = "^ms.ms.count$", numeric = "^reporter.intensity."))
    if (class(yaml_contaminants) == "list") {
      if (exists("df_pg")) {
        lst_qcMetrics[["qcMetric_EVD_UserContaminant"]]$setData(df_evd, 
                                                                df_pg, yaml_contaminants)
      }
      else {
        lst_qcMetrics[["qcMetric_EVD_UserContaminant"]]$setData(df_evd, 
                                                                NULL, yaml_contaminants)
      }
    }
    lst_qcMetrics[["qcMetric_EVD_PeptideInt"]]$setData(df_evd, 
                                                       param_EV_intThresh)
    if (length(grep("^reporter.intensity.", colnames(df_evd))) > 
        0) {
      lst_qcMetrics[["qcMetric_EVD_ReporterInt"]]$setData(df_evd)
    }
    df_evd$hasMTD = !is.na(df_evd$match.time.difference)
    reportMTD = any(df_evd$hasMTD)
    lst_qcMetrics[["qcMetric_EVD_ProteinCount"]]$setData(df_evd, 
                                                         param_EV_protThresh)
    lst_qcMetrics[["qcMetric_EVD_PeptideCount"]]$setData(df_evd, 
                                                         param_EV_pepThresh)
    if ("retention.length" %in% colnames(df_evd)) {
      lst_qcMetrics[["qcMetric_EVD_RTPeakWidth"]]$setData(df_evd)
    }
    if (("retention.time.calibration" %in% colnames(df_evd))) {
      MBR_HAS_DATA = (sum(df_evd$type == "MULTI-MATCH") > 
                        0)
      if ((param_evd_mbr == FALSE) || (MBR_HAS_DATA == 
                                       FALSE)) {
      }
      else {
        lst_qcMetrics[["qcMetric_EVD_MBRAlign"]]$setData(df_evd, 
                                                         param_EV_MatchingTolerance, mq$raw_file_mapping)
        avg_peak_width = lst_qcMetrics[["qcMetric_EVD_RTPeakWidth"]]$outData[["avg_peak_width"]]
        if (is.null(avg_peak_width)) {
          stop("RT peak width module did not run, but is required for MBR metrics. Enable it and try again or switch off MBR metrics!")
        }
        lst_qcMetrics[["qcMetric_EVD_MBRIdTransfer"]]$setData(df_evd, 
                                                              avg_peak_width)
        lst_qcMetrics[["qcMetric_EVD_MBRaux"]]$setData(df_evd)
      }
    }
    lst_qcMetrics[["qcMetric_EVD_Charge"]]$setData(df_evd)
    lst_qcMetrics[["qcMetric_EVD_IDoverRT"]]$setData(df_evd)
    if (enabled_summary) {
      df_idrate = d_smy$raw[, c("fc.raw.file", "ms.ms.identified....")]
    }
    else {
      df_idrate = NULL
    }
    lst_qcMetrics[["qcMetric_EVD_PreCal"]]$setData(df_evd, 
                                                   df_idrate, param_EV_PrecursorTolPPM, param_EV_PrecursorOutOfCalSD)
    lst_qcMetrics[["qcMetric_EVD_PostCal"]]$setData(df_evd, 
                                                    df_idrate, param_EV_PrecursorTolPPM, param_EV_PrecursorOutOfCalSD, 
                                                    param_EV_PrecursorTolPPMmainSearch)
    lst_qcMetrics[["qcMetric_EVD_Top5Cont"]]$setData(df_evd)
    lst_qcMetrics[["qcMetric_EVD_MS2OverSampling"]]$setData(df_evd)
    lst_qcMetrics[["qcMetric_EVD_MissingValues"]]$setData(df_evd)
    if (!DEBUG_PTXQC) 
      df_evd = df_evd[, c("id", "contaminant")]
  }
  if (enabled_msms) {
    df_msms = mq$readMQ(txt_files$msms, type = "msms", filter = "", 
                        col_subset = c(numeric = "Missed\\.cleavages", "^Raw.file$", 
                                       "^mass.deviations", "^masses$", "^mass.analyzer$", 
                                       "fragmentation", "reverse", numeric = "^evidence.id$"), 
                        check_invalid_lines = FALSE)
    lst_qcMetrics[["qcMetric_MSMS_MSMSDecal"]]$setData(df_msms = df_msms, 
                                                       fc_raw_files = mq$raw_file_mapping$to)
    if (exists("df_evd")) {
      lst_qcMetrics[["qcMetric_MSMS_MissedCleavages"]]$setData(df_msms, 
                                                               df_evd)
    }
    else {
      lst_qcMetrics[["qcMetric_MSMS_MissedCleavages"]]$setData(df_msms)
    }
    rm(df_msms)
    if (!DEBUG_PTXQC) 
      rm(df_evd)
  }
  if (enabled_msmsscans) {
    df_msmsScan = mq$readMQ(txt_files$msmsScan, type = "msms", 
                            filter = "", col_subset = c(numeric = "^ion.injection.time", 
                                                        numeric = "^retention.time$", "^Identified", 
                                                        "^Scan.event.number", "^total.ion.current", "^base.?peak.intensity", 
                                                        "^Raw.file", "^dp.aa$", "^dp.modification$"), 
                            check_invalid_lines = FALSE)
    if (ncol(df_msmsScan) > 3) {
      df_msmsScan$rRT = round(df_msmsScan$retention.time/2) * 
        2
      lst_qcMetrics[["qcMetric_MSMSScans_TopNoverRT"]]$setData(df_msmsScan)
      lst_qcMetrics[["qcMetric_MSMSScans_IonInjTime"]]$setData(df_msmsScan, 
                                                               param_MSMSScans_ionInjThresh)
      lst_qcMetrics[["qcMetric_MSMSScans_MSMSIntensity"]]$setData(df_msmsScan)
      lst_qcMetrics[["qcMetric_MSMSScans_TopN"]]$setData(df_msmsScan)
      lst_qcMetrics[["qcMetric_MSMSScans_TopNID"]]$setData(df_msmsScan)
      if ("dp.modification" %in% colnames(df_msmsScan)) {
        lst_qcMetrics[["qcMetric_MSMSScans_DepPep"]]$setData(df_msmsScan)
      }
    }
    rm(df_msmsScan)
  }
  print("#Metrics: ")
  print(length(lst_qcMetrics))
  hm = getQCHeatMap(lst_qcMetrics, raw_file_mapping = mq$raw_file_mapping)
  write.table(hm[["table"]], file = rprt_fns$heatmap_values_file, 
              quote = TRUE, sep = "\t", row.names = FALSE)
  pl_nameMapping = mq$plotNameMapping()
  cat("Creating Report file ...")
  out_formats = unlist(strsplit(param_OutputFormats, "[ ,]+"))
  out_formats
  out_format_requested = out_formats_supported[match(out_formats, 
                                                     out_formats_supported)]
  if (any(is.na(out_format_requested))) {
    stop("Output format(s) not supported: '", paste(out_formats[is.na(out_format_requested)], 
                                                    collapse = "', '"), "'")
  }
  if ("html" %in% out_format_requested) {
    if (pandoc_available()) {
      if (DEBUG_PTXQC) {
        html_template = "Z:/projects/QC/PTXQC/package/inst/reportTemplate/PTXQC_report_template.Rmd"
      }
      else {
        html_template = system.file("./reportTemplate/PTXQC_report_template.Rmd", 
                                    package = "PTXQC")
      }
      cat(paste0("HTML TEMPLATE: ", html_template, "\n"))
      out_dir = dirname(rprt_fns$report_file_HTML)
      file.copy(html_template, out_dir, overwrite = TRUE)
      out_template = file.path(out_dir, basename(html_template))
      render(out_template, output_file = rprt_fns$report_file_HTML)
    }
    else {
      warning("The 'Pandoc' converter is not installed on your system or you do not have read-access to it!\n", 
              "Pandoc is required for HTML reports.\n", "Please install Pandoc <http://pandoc.org/installing.html> or make sure you have access to pandoc(.exe).\n", 
              "Restart your R-session afterwards.", immediate. = TRUE)
    }
  }
  if ("plainPDF" %in% out_format_requested) {
    report_file_PDF = rprt_fns$report_file_PDF
    if (!wait_for_writable(report_file_PDF)) {
      stop("Target file not writable")
    }
    if (param_PageNumbers == "on") {
      printWithPage = function(gg_obj, page_nr, filename = report_file_PDF) {
        filename = basename(filename)
        printWithFooter(gg_obj, bottom_left = filename, 
                        bottom_right = page_nr)
      }
    }
    else {
      printWithPage = function(gg_obj, page_nr, filename = report_file_PDF) {
        print(gg_obj)
      }
    }
    pdf(report_file_PDF)
    printWithPage(hm[["plot"]], "p. 1")
    printWithPage(pl_nameMapping$plots, "p. 2")
    pc = 3
    for (qcm in lst_qcMetrics_ord) {
      for (p in qcm$plots) {
        printWithPage(p, paste("p.", pc))
        pc = pc + 1
      }
    }
    dev.off()
    cat(" done\n")
  }
  mq$writeMappingFile(rprt_fns$filename_sorting)
  cat(paste("Report file created at\n\n    ", rprt_fns$report_file_prefix, 
            ".*\n\n", sep = ""))
  cat(paste0("\n\nTime elapsed: ", round(as.double(Sys.time() - 
                                                     time_start, units = "mins"), 1), " min\n\n"))
  return(rprt_fns)
}
