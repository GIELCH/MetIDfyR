###############################################################
# Project : MetIDfyR                                          
# Author  : Awa Bousso GUEYE                                  
# Date    : 13/08/2026                                      
# Version : v2.0 (corrigee + logging structure)
#
# CORRECTIONS DE BUGS (cf version precedente MetIDfyR_corrige.R) :
#  FIX 1 : filename_tsv recalcule a chaque pic (evite l'ecrasement
#          de table_<formule>_<index>_..._.tsv par un pic suivant).
#  FIX 2 : condition d'ecriture basee sur nrow(data_ms2$match) > 0
#          plutot que length(data_ms2) > 0 (presque toujours vrai).
#  FIX 3 : table vide ecrite meme si ms2_ref est NULL.
#  FIX 4 : transfo$name <- transfo$Transformation (compareMS2()
#          cherchait `transfo$name`, absent -> mass_change = NA pour
#          toute transformation nommee).
#  FIX 5 : out_current_mlc normalise en chemin absolu avant tout
#          makeCluster()/clusterExport().
#  FIX 6 : bloc de cluster mort (sans foreach) supprime.
#  FIX 7 : "save_structure_svg" retire de p_load() sur les workers
#          (c'est une fonction, pas un package).
#
# AMELIORATIONS "PUBLICATION-READY" AJOUTEES DANS CETTE VERSION :
#  LOG 1 : remplacement de cat()/print()/message()/warning() par le
#          package `logger` (niveaux INFO/WARN/ERROR), sortie
#          simultanee console + fichier horodate unique par run.
#  LOG 2 : logging en environnement parallele -> chaque worker ecrit
#          dans son propre fichier (nomme par PID), fusionnes dans le
#          log principal apres le foreach.
#  GLOB 1 : clusterExport(cl_big, ls()) remplace par une liste EXPLICITE
#          des objets reellement necessaires aux workers (plus sur,
#          plus lisible, evite d'exporter des objets non pertinents/
#          non serialisables).
#  ERR 1  : les etapes a risque (BioTransformer, ecriture fichiers)
#          sont enveloppees dans tryCatch() avec log_error() avant de
#          passer a la molecule suivante, au lieu de laisser l'erreur
#          remonter silencieusement ou faire planter tout le run.
#  PATH 1 : tous les chemins utilises viennent de opt$output / opt$config /
#          objets deja parametres -> aucun chemin absolu code en dur
#          (deja le cas dans cette version du script, verifie et documente).
#
# NOTE : les points restants de la checklist "publication" (tests
# testthat, package R structure, renv/Docker, CI, DOI Zenodo) sont
# des livrables complementaires (nouveaux fichiers / structure de
# projet), pas des corrections internes a ce script -> a construire
# a part si besoin.
###############################################################


#### DEPENDENCIES ####

cat(" Working directory dans MetIDfyR.R : ", getwd(), "\n")
if(!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load("BiocManager", "optparse", "logger")

# Install dependencies from Bioconductor
if(!"MSnbase" %in% installed.packages()) BiocManager::install("MSnbase")
if(!"Rdisop" %in% installed.packages()) BiocManager::install("Rdisop")
if(!"reticulate" %in% installed.packages()) BiocManager::install("reticulate")

#### PARAMETERS ####

option_list = list(
  make_option(c("-i", "--input"), type="character",
              help="input informations (see TEMPLATE_start_mlc.tsv)", 
              metavar="character"),
  make_option(c("-o", "--output"), type="character", default = "out_MetIDfyR",
              help="output directory", 
              metavar="character"),
  make_option(c("-c", "--config"), type="character", default="input/config.R",
              help="configuration file (see TEMPLATE_config.R)", 
              metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Test option
if (is.null(opt$input )){
  print_help(opt_parser)
  stop("Input tsv required.n", call.=FALSE)
}
if (is.null(opt$output )){
  print_help(opt_parser)
  stop("Output directory required.n", call.=FALSE)
}

dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

#### LOG 1 : LOGGER SETUP (avant tout autre traitement) ####
# Chemin absolu + horodate : un seul fichier de log principal par run,
# facile a joindre en annexe/supplement d'un article.
log_dir  <- normalizePath(opt$output, mustWork = FALSE)
log_file <- file.path(log_dir, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))

log_appender(appender_tee(file = log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_threshold(INFO)

log_info("=== Demarrage MetIDfyR v2 ===")
log_info("Input   : {opt$input}")
log_info("Output  : {opt$output}")
log_info("Config  : {opt$config}")

### LOG 2 : petit helper pour les logs cote worker ###
# Les workers %dopar% sont des processus separes : ils n'ecrivent pas
# dans le meme flux que le processus principal. Chaque worker log dans
# son propre fichier (par PID) ; ces fichiers sont fusionnes dans le
# log principal une fois le foreach termine (voir plus bas).
worker_log <- function(msg, log_dir, level = "INFO"){
  cat(sprintf("%s [%s] [PID:%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), level, Sys.getpid(), msg),
      file = file.path(log_dir, paste0("worker_", Sys.getpid(), ".log")),
      append = TRUE)
}

#### CREATION D'UN ENVIRONNEMENT VIRTUEL ####
library("reticulate")

use_virtualenv("rdkit-env", required = TRUE) 
import("rdkit")
py_module_available("rdkit")


#### SOURCES ####
source_python("diff.py")
source(opt$config) #load config informations
source("util_smiles-copie.R") #load functions file
source("plot_functions.R") #load plot functions file


.libPaths( c( .libPaths(), lib_perso))
suppressPackageStartupMessages(pacman::p_load("dplyr", "plyr", "doParallel", "stringr", "MSnbase", "Rdisop", 
                                              "readr", "tibble", "arrangements", "snow", "svglite"))

options(digits=10) #allow more digits
E_MASS = c(minus = 5.489E-4, plus = -5.489E-4) 


#### Initialisation : Data loading and preparation ####

log_info("### Initialisation ###")
data_tsv = read_tsv(opt$input, col_types = cols()) #table containing the molecules
data_tsv = as_tibble(cbind(data_tsv, fill_table(data_tsv)))

#### VERIFICATION DE LA COLONNE SMILES ####
if(!"smiles" %in% tolower(names(data_tsv))){
  log_error("Colonne 'smiles' absente du fichier d'entree {opt$input}")
  stop("Error: the 'smiles' column is missing. Cannot proceed without chemical structure (SMILES).")
}

## Gestion 
sanitize <- function(x) gsub("[^A-Za-z0-9_-]", "_",x)
bt_dir <- file.path(opt$output, "BioTransformer")

log_info("### Configuration ###")
log_info(paste(readLines(opt$config), collapse = "\n"))
log_info("{nrow(data_tsv)} molecule(s) a traiter")

#### === LOOP OVER MOLECULES === ####

for(row in 1:nrow(data_tsv)){
  
  start = Sys.time()
  
  current_target = data_tsv[row, ] ; do = c(); ms_file = c(); adduct = c() 
  
  log_info("--- Molecule {row}/{nrow(data_tsv)} : {current_target$name} ---")
  
  do$plus = !is.na(current_target$adduct_plus) ; do$minus = !is.na(current_target$adduct_minus)
  
  if(!do$plus) current_target$adduct_plus = FALSE ; if(!do$minus) current_target$adduct_minus = FALSE
  if(is.na(current_target$ms2_reference_tsv)) current_target$ms2_reference_tsv = FALSE
  
  #### === PREDICTION BIOTRANSFORMER === ####
  
  if (bool_phase_1 & !bool_phase_2) {
    bt_mode <- "cyp450"
  } else if (!bool_phase_1 & bool_phase_2) {
    bt_mode <- "phaseII"
  } else if (bool_phase_1 & bool_phase_2) {
    bt_mode <- "allHuman"
  } else {
    log_error("Aucune phase (1 ou 2) selectionnee dans la configuration")
    stop("Veuillez selectionner au moins une phase (1 ou 2).")
  }
  
  log_info("Mode BioTransformer selectionne : {bt_mode}")
  
  ### ERR 1 : BioTransformer peut echouer (jar manquant, timeout, SMILES
  ### invalide...) -> on capture l'erreur, on logue, et on passe a la
  ### molecule suivante au lieu de faire planter tout le run.
  bt_df <- tryCatch({
    run_biotransformer(
      smiles = current_target$smiles,
      name = current_target$name,
      bt_mode = bt_mode,
      bt_dir = bt_dir,
      bt_jar = "biotransformer-1.1.5.jar"
    )
  }, error = function(e){
    log_error("BioTransformer a echoue pour {current_target$name} : {conditionMessage(e)}")
    return(NULL)
  })
  
  if (is.null(bt_df)) {
    log_warn("Aucun resultat BioTransformer pour {current_target$name}, molecule ignoree")
    next
  }
  
  results <- tryCatch({
    process_biotransformer_results(bt_df, current_target, bt_dir, do_smiles)
  }, error = function(e){
    log_error("Echec du traitement des resultats BioTransformer pour {current_target$name} : {conditionMessage(e)}")
    return(NULL)
  })
  
  if (is.null(results) || is.null(results$info_all_combi) || nrow(results$info_all_combi) == 0) {
    log_warn("Aucune transformation pour {current_target$name}, passage a la molecule suivante")
    next
  }
  
  info_all_combi <- results$info_all_combi
  transfo_ref <- unique(info_all_combi[, c("Transformation", "Diff", 'Formula')])
  
  transfo_ref <- calculate_mass_changes(transfo_ref, current_target)
  
  ### FIX 4 ###
  # compareMS2() (plot_functions.R) fait match(x, transfo$name) : on ajoute
  # l'alias `name` attendu, sans casser le reste du script (`Transformation`).
  transfo <- transfo_ref
  transfo$name <- transfo$Transformation
  
  log_info("{nrow(transfo)} transformation(s) unique(s) identifiee(s) pour {current_target$name}")
  
  #### Chromatogram : check signal for each metabolites ####
  
  log_info("### Debut chromatogramme pour {current_target$name} ###")
  
  out_current_mlc <- file.path(
    opt$output,
    paste0("out_", sanitize(current_target$name), "_", 
           sanitize(sub('\\..*$', '', basename(current_target$ms_file))), "_", 
           format(Sys.time(), "%d%m%y_%H%M"))
  )
  
  ### FIX 5 ###
  # Chemin absolu AVANT tout makeCluster()/clusterExport(), pour que les
  # workers ecrivent (donnees ET logs) au meme endroit que le process principal.
  out_current_mlc <- normalizePath(out_current_mlc, mustWork = FALSE)
  log_info("Repertoire de sortie : {out_current_mlc}")
  
  lapply(c("", "/POS", "/NEG", "/POSNEG", "/input"), function(newdir) {
    dir.create(file.path(out_current_mlc, newdir), showWarnings = FALSE, recursive = TRUE)
  })
  
  ### FIX 6 : bloc de cluster mort (makeCluster + clusterExport + stopCluster
  ### sans foreach entre les deux) retire.
  
  log_info("### Lecture {current_target$ms_file} ###")
  log_info("### Recuperation de la reference MS2 ###")
  
  ms_file = c() ; adduct = c()
  
  if(do$minus){
    adduct$minus = list(formula = current_target$adduct_minus,
                        mz = ifelse(current_target$adduct_minus=="H",
                                    -getMolecule(current_target$adduct_minus)$isotopes[[1]][1,1],
                                    getMolecule(current_target$adduct_minus)$isotopes[[1]][1,1]))
    ms_file$minus = tryCatch(
      readMSData(current_target$ms_file, mode="onDisk") %>% filterPol(polarity = 0),
      error = function(e){
        log_error("Lecture MS impossible (polarite negative) pour {current_target$name} : {conditionMessage(e)}")
        NULL
      }
    )
  }
  
  if(do$plus){
    adduct$plus = list(formula = current_target$adduct_plus,
                       mz = getMolecule(current_target$adduct_plus)$isotopes[[1]][1,1])
    ms_file$plus = tryCatch(
      readMSData(current_target$ms_file, mode="onDisk") %>% filterPol(polarity = 1),
      error = function(e){
        log_error("Lecture MS impossible (polarite positive) pour {current_target$name} : {conditionMessage(e)}")
        NULL
      }
    )
  }
  
  if((do$minus && is.null(ms_file$minus)) || (do$plus && is.null(ms_file$plus))){
    log_warn("Fichier MS illisible pour {current_target$name}, molecule ignoree")
    next
  }
  
  if(current_target$ms2_reference_tsv != FALSE ) {
    if(file.exists(current_target$ms2_reference_tsv)){
      ref_ms2_tsv = read_tsv(current_target$ms2_reference_tsv, col_types = cols())
      ref_ms2_tsv$polarity = ifelse(ref_ms2_tsv$intensity > 0, "plus", "minus")
    }else{
      log_warn("Le fichier MS2 de reference fourni n'existe pas : {current_target$ms2_reference_tsv}")
    }
  }
  
  ref_ms2 = getMS2Reference(current_target, names(do)[which(do==T)])
  
  # Save input for the current molecule
  file.copy(opt$config, file.path(out_current_mlc, "input"))
  write_tsv(data_tsv[row, c(1:6)], file.path(out_current_mlc, "input", "molecule.tsv"))
  if(exists("ref_ms2_tsv")) file.copy(current_target$ms2_reference_tsv, file.path(out_current_mlc, "input"))
  
  
  #### PARALLEL LOOP ####
  
  freeCores <- parallel::detectCores()
  if(!exists("cores") || is.na(cores)){
    log_warn("`cores` non defini dans la config, utilisation de freeCores - 1")
    cores <- freeCores - 1
  }
  cores <- ifelse(cores >= freeCores | is.na(cores), freeCores - 1, cores)
  log_info("Utilisation de {cores} coeur(s) sur {freeCores} disponible(s)")
  
  cl_big = makeCluster(cores)
  clusterApply(cl_big, 1:cores, opt, fun = function(x, opt){
    library(reticulate)
    use_virtualenv("rdkit-env", required = TRUE)
    import("rdkit")
    py_module_available("rdkit")
    source(opt$config)
    source("util_smiles-copie.R")
    source("plot_functions.R")
    library(pacman)
    
    ### FIX 7 ###
    # "save_structure_svg" est une fonction (plot_functions.R, deja
    # source() ci-dessus), pas un package -> retiree de p_load().
    p_load("tibble", "MSnbase", "plyr", "ggplot2", "ggpubr", "Rdisop", "dplyr", "ggrepel",
           "stringr", "htmlwidgets", "readr")
  })
  
  ### GLOB 1 ###
  # clusterExport(cl_big, ls()) exportait TOUT l'environnement local
  # (y compris des objets volumineux ou non pertinents pour le worker).
  # On exporte desormais explicitement seulement ce qui est utilise dans
  # le corps du foreach ci-dessous. Les FONCTIONS utilisees (getChromato,
  # doMS2, compareMS2, close_match, save_structure_svg, ...) sont deja
  # disponibles sur chaque worker via les source() executes dans
  # clusterApply() ci-dessus -> pas besoin de les exporter.
  vars_to_export <- c(
    "info_all_combi", "current_target", "ms_file", "adduct", "do", "ref_ms2",
    "E_MASS", "transfo", "mz_ppm", "min_peak_intensity", "nb_scan", "minimum_mz",
    "rt_windows", "wdw_mz_ms2", "out_current_mlc", "worker_log"
  )
  missing_vars <- setdiff(vars_to_export, ls())
  if(length(missing_vars) > 0){
    log_warn("Variables attendues absentes avant clusterExport : {paste(missing_vars, collapse=', ')}")
  }
  clusterExport(cl_big, intersect(vars_to_export, ls()))
  
  registerDoParallel(cl_big)
  
  info_all_combi$metabolite_id <- paste(info_all_combi$SMILES, info_all_combi$Transformation, sep = "|")
  
  log_info("Lancement du traitement parallele sur {length(unique(info_all_combi$metabolite_id))} candidat(s)")
  
  BIG_TABLE = foreach::foreach(current_id = unique(info_all_combi$metabolite_id), .combine = rbind) %dopar% {
    
    parts <- strsplit(current_id, "\\|")[[1]]
    current_smiles <- parts[1]
    current_transfo <- parts[2]
    
    current_formula <- info_all_combi$Formula[info_all_combi$SMILES == current_smiles & info_all_combi$Transformation == current_transfo][1]
    current_diff <- unique(info_all_combi$Diff[info_all_combi$SMILES == current_smiles & info_all_combi$Transformation == current_transfo])
    nb_transfo <- info_all_combi$Nb_Transfo[info_all_combi$SMILES == current_smiles & info_all_combi$Transformation == current_transfo][1]
    
    worker_log(sprintf("Debut candidat : formule=%s transfo=%s", current_formula, current_transfo), out_current_mlc)
    
    current_mlc = getMolecule(current_formula)        
    if(getValid(current_mlc)=="Valid"){
      
      plot_subtitle = paste("Shorter transformation :", current_transfo,
                            ifelse(current_diff != "", paste("\nFormula differences :", current_diff), ""))
      
      combine_plot_chromato = data.frame() ; combine_ms_data = data.frame()
      
      for(pol in c("minus", "plus")){
        if(do[[pol]]){
          plot_chromato = tryCatch(
            getChromato(ms_file[[pol]], current_mlc, adduct[[pol]], E_MASS[[pol]], min_intensity = min_peak_intensity,
                        mz_precision = mz_ppm, min_scan = nb_scan, min_mz = minimum_mz),
            error = function(e){
              worker_log(sprintf("getChromato() a echoue pour %s (%s) : %s", current_formula, pol, conditionMessage(e)),
                         out_current_mlc, level = "ERROR")
              data.frame()
            }
          )
          
          if(!empty(plot_chromato)){
            ms_data = tryCatch(
              getMassSpectrum(ms_file[[pol]], current_mlc, adduct[[pol]], plot_chromato, E_MASS[[pol]], peak_wdw = rt_windows,
                              min_intensity = min_peak_intensity),
              error = function(e){
                worker_log(sprintf("getMassSpectrum() a echoue pour %s (%s) : %s", current_formula, pol, conditionMessage(e)),
                           out_current_mlc, level = "ERROR")
                data.frame()
              }
            )
            
            if(!empty(ms_data)){
              plot_chromato$index = 0 ; plot_chromato$polarity = pol
              plot_chromato$index[which((plot_chromato$rt/60) %in% unique(ms_data$rtime) &
                                          plot_chromato$isotope == 1)] = unique(ms_data$index)
              ms_data$polarity = pol
              combine_plot_chromato = rbind(combine_plot_chromato, plot_chromato)
              combine_ms_data = rbind(combine_ms_data, ms_data)
            }
          }
        }
      }
      
      
      ### BUILD AND SAVE FIGURES ###
      
      if(!empty(combine_plot_chromato) && "polarity" %in% colnames(combine_plot_chromato)){
        polarities <- unique(combine_plot_chromato$polarity)
        if(length(polarities) == 2){
          out_polarity = "POSNEG"
        }else if(polarities == "minus"){
          out_polarity = "NEG"
        }else if(polarities == "plus") {
          out_polarity = "POS"
        }else{
          out_polarity <- "UNKNOWN"
        }
        
        plot_path = paste0(out_current_mlc, "/", out_polarity)
        image_path = paste0(out_current_mlc, "/images/", out_polarity)
        dir.create(image_path, showWarnings = FALSE, recursive = TRUE)
        
        tmp_table = data.frame()
        
        for(pol in unique(combine_plot_chromato$polarity)){
          perc_common_peak = c() ; dotp_ms2 = c() ; mono_ppm = c()
          
          ms2_ref = ref_ms2[[pol]][2]  ;  mz_ref = ref_ms2[[pol]][1]
          
          plot_chromato = combine_plot_chromato[combine_plot_chromato$polarity == pol, ]
          ms_data = combine_ms_data[which(combine_ms_data$polarity == pol), ]
          
          group_ms = split(ms_data, ms_data$rtime)
          
          filenames = c()
          
          for(current_ms in group_ms) {
            
            ggp = doChromato(plot_chromato[which(plot_chromato$rt/60 < unique(current_ms$rtime) +2 &
                                                   plot_chromato$rt/60 > unique(current_ms$rtime) -2), ]) +
              labs(title = "Metabolite chromatogram",
                   subtitle = paste(plot_subtitle, "\nm/z : ", round(min(plot_chromato$mz),5)))
            
            ggp_ms = doMassSpectrum(current_ms)
            mono_ppm = c(mono_ppm, current_ms$exp_ppm[!is.na(current_ms$exp_ppm)][1])
            
            ggparr_ms = ggarrange(ggp, ggp_ms)
            
            ### FIX 1 + 2 ###
            filename_tsv = paste0("table_", current_formula, "_", unique(current_ms$index), "_",
                                  ifelse(out_polarity=="POSNEG", pol, ""), "_", out_polarity, ".tsv")
            
            if( 2 %in% unique(msLevel(ms_file[[pol]])) ){
              if(!is.null(ms2_ref)){
                data_ms2 = tryCatch(
                  compareMS2(ms_file[[pol]], current_transfo, ms2_ref[[1]], mz_ref[[1]], unique(plot_chromato$mz),
                             unique(current_ms$rtime), wdw_mz = wdw_mz_ms2, tol_mz = mz_ppm),
                  error = function(e){
                    worker_log(sprintf("compareMS2() a echoue pour %s : %s", current_formula, conditionMessage(e)),
                               out_current_mlc, level = "ERROR")
                    list()
                  }
                )
                ggp_ms2 = doMS2(data_ms2)
                dotp_ms2 = c(dotp_ms2, ifelse(is.null(data_ms2$dotp_ms2), 0, data_ms2$dotp_ms2))
                
                ### FIX 2 ###
                if(length(data_ms2) > 1 && !is.null(data_ms2$match) && nrow(data_ms2$match) > 0){
                  write_tsv(data_ms2$match, file.path(plot_path, filename_tsv))
                } else {
                  empty_table <- data.frame(note = "No MS2 reference available")
                  write_tsv(empty_table, file.path(plot_path, filename_tsv))
                }
              }else{
                data_ms2 = compareMS2(ms_file[[pol]], current_transfo, mz_exp = unique(plot_chromato$mz)[1],
                                      rt_exp = unique(current_ms$rtime), exists_ref_ms2 = F, wdw_mz = wdw_mz_ms2)
                ggp_ms2 = doMS2(data_ms2, exists_ref = F)
                dotp_ms2 = c(dotp_ms2, 0)
                
                ### FIX 3 ###
                empty_table <- data.frame(note = "No MS2 reference available for this polarity")
                write_tsv(empty_table, file.path(plot_path, filename_tsv))
              }
              ggp_tot = ggarrange(ggparr_ms, ggp_ms2, nrow=2)
            }else{
              data_ms2 = c()
              dotp_ms2 = c(dotp_ms2, 0)
              ggp_tot = ggparr_ms
            }
            
            perc_peak = 0
            if(length(data_ms2) > 1 && nrow(data_ms2$match) > 0) {
              nb_ref_peaks = length(grep("ref", data_ms2$data$type))
              perc_peak = round(max(data_ms2$match$index_peak) / nb_ref_peaks*100,2)
              
              ggp_tot = annotate_figure(ggp_tot, bottom = text_grob(
                paste0("Percentage of common peaks : ", perc_peak, "%", " / Dotp : ", round(dotp_ms2[length(dotp_ms2)],3)), size = 8 ))
            }
            
            ggp_tot = annotate_figure(ggp_tot, top = text_grob(
              paste(current_target$name, ":", current_formula,
                    ifelse(pol == "plus", "+", "-"), 
                    adduct[[pol]]$formula,
                    "/ peak number", unique(current_ms$index))
              , size = 9))
            
            perc_common_peak = c(perc_common_peak, perc_peak)
            
            filename = paste0(current_formula, "_", unique(current_ms$index),"_",
                              ifelse(out_polarity=="POSNEG", pol, ""), "_", out_polarity,".svg")
            ggsave(file.path(plot_path, filename), ggp_tot)
            
            filenames = append(filenames, filename)
            
            if(!is.na(current_smiles) && current_smiles != "") {
              structure_filename = paste0(current_smiles, "_structure.svg")
              structure_path = file.path(image_path, structure_filename)
              save_structure_svg(current_smiles, structure_path, width=400, height=400)
            }
            
          }  # end for retention time 
          peak_info = unique(ms_data[, c("rtime", "dotp", "rscore", "abscore", "peak_intensity", "index")])
          
          tmp_table = rbind.data.frame(tmp_table,
                                       cbind.data.frame(name = current_target$name, formula = current_formula,
                                                        polarity = pol,
                                                        smiles=current_smiles,
                                                        adduct = unique(plot_chromato$adduct), mz = round(min(plot_chromato$mz),5),
                                                        transfo = current_transfo, diff = current_diff, rt = round(peak_info$rtime, 3),
                                                        abscore = round(peak_info$abscore, 3), dotp_ms2 = round(dotp_ms2, 3), 
                                                        common_ms2_peak = perc_common_peak, mono_ppm = mono_ppm,
                                                        intensity = peak_info$peak_intensity, index_peak = peak_info$index,
                                                        nb_transfo = nb_transfo,
                                                        filepath = plot_path, filename = filenames))
          
          
        } #end for polarity
        
        tmp_table
        
      }
    } else {
      worker_log(sprintf("Molecule invalide (formule=%s)", current_formula), out_current_mlc, level = "WARN")
    }
    
  } #end foreach
  
  
  stopCluster(cl_big)
  closeAllConnections()
  
  ### LOG 2 : fusion des logs workers dans le log principal ###
  worker_logs <- list.files(out_current_mlc, pattern = "^worker_.*\\.log$", full.names = TRUE)
  if(length(worker_logs) > 0){
    file.append(log_file, worker_logs)
    file.remove(worker_logs)
    log_info("{length(worker_logs)} fichier(s) de log worker fusionne(s) pour {current_target$name}")
  }
  
  
  if(!empty(BIG_TABLE)){
    log_info("{nrow(BIG_TABLE)} ligne(s) de resultats bruts pour {current_target$name}")
    BIG_TABLE$nb_transfo = as.numeric(as.character(BIG_TABLE$nb_transfo))
    
    ## Compute score
    
    parent_intensity = list(plus="", minus="")
    for(pol in c("plus", "minus")){
      if(do[[pol]]){
        parent = which(BIG_TABLE$nb_transfo==0 & BIG_TABLE$polarity==pol)
        if(length(parent) > 1){
          if( length(ref_ms2[[pol]]) == 3 && !is.na(ref_ms2[[pol]][3]) ){
            good_parent = which.min( abs( BIG_TABLE$rt[parent] - ref_ms2[[pol]][[3]] ) )
            parent_intensity[[pol]] = BIG_TABLE$intensity[parent[good_parent]]
          }else{
            good_parent = which.max(BIG_TABLE$intensity[parent])
            parent_intensity[[pol]] = BIG_TABLE$intensity[parent[good_parent]]
          }
        }else if(length(parent) == 0){
          parent_intensity[[pol]] = 0
          log_warn("Aucun signal parent trouve pour {current_target$name} en {pol}")
        }else{
          parent_intensity[[pol]] = BIG_TABLE$intensity[parent]
        }
      }
    }
    
    ratio_intensity = ifelse(BIG_TABLE$polarity == "plus",
                             ifelse(BIG_TABLE$intensity/parent_intensity$plus == Inf,
                                    0, BIG_TABLE$intensity/parent_intensity$plus),
                             ifelse(BIG_TABLE$intensity/parent_intensity$minus == Inf,
                                    0, BIG_TABLE$intensity/parent_intensity$minus))
    ratio_intensity[ratio_intensity>1] = 1
    
    penalty = ifelse(BIG_TABLE$nb_transfo == 0, 1, 1/sqrt(BIG_TABLE$nb_transfo))
    BIG_TABLE$rintensity = ratio_intensity
    BIG_TABLE$score = round(penalty * ( 1/2 * BIG_TABLE$common_ms2_peak / 100 * BIG_TABLE$dotp_ms2 +
                                          1/2 * BIG_TABLE$abscore * ratio_intensity )
                            , 3 )
    
    BIG_TABLE$score[BIG_TABLE$score == Inf] = 1
    BIG_TABLE = arrange(BIG_TABLE, desc(score))
    
    # Remove isobare (close M/Z)
    for(tab_row in 1:nrow(BIG_TABLE)){
      isobare = close_match(x= BIG_TABLE$mz, target=BIG_TABLE$mz[tab_row], tolerance = BIG_TABLE$mz[tab_row]/1e6*mz_ppm)
      true_isobare = which(BIG_TABLE$formula[isobare] != BIG_TABLE$formula[tab_row] & 
                             BIG_TABLE$polarity[isobare] == BIG_TABLE$polarity[tab_row])
      if(length(true_isobare) > 0) BIG_TABLE = BIG_TABLE[-isobare[true_isobare],]
    }
    
    BIG_TABLE[,c("formula", "filepath", "filename")] = apply(BIG_TABLE[,c("formula", "filepath", "filename")], 2, as.character)
    # Group by RT and remove isotopes false positive
    BIG_TABLE_FINAL = split(BIG_TABLE, BIG_TABLE$rt) %>%
      lapply(function(current_rt){
        
        index = 1
        while (index <= nrow(current_rt)) {
          mol = getMolecule(current_rt$formula[index])
          mz_iso = mol$isotopes[[1]][1, 2:6] + adduct[[current_rt$polarity[index]]]$mz + E_MASS[[current_rt$polarity[index]]]
          false_pos = unlist( lapply(mz_iso, function(iso) which(between(current_rt$mz, iso-iso/1e6*mz_ppm, iso+iso/1e6*mz_ppm))) )
          
          names = unlist(strsplit(current_rt$filename[false_pos], ".svg"))
          file.remove(
            list.files(path = unique(current_rt$filepath[false_pos]), pattern = paste0(names, collapse="|"), full.names = T)
          )
          if(length(false_pos) > 0){
            current_rt = current_rt[-false_pos,]
          }
          index = index +1
        }
        
        return(current_rt[,!names(current_rt) %in% c("filepath", "filename")])
        
      }) 
    
    BIG_TABLE_FINAL = do.call(rbind, BIG_TABLE_FINAL)
    
    BIG_TABLE_FINAL$polarity = ifelse(BIG_TABLE_FINAL$polarity == "plus", "Positive mode", "Negative mode")
    
    out_tsv <- file.path(out_current_mlc, paste0("out_", sub('\\..*$', '', basename(current_target$ms_file)), ".tsv"))
    BIG_TABLE_FINAL %>% write_tsv(out_tsv)
    log_info("Resultats finaux ecrits : {out_tsv} ({nrow(BIG_TABLE_FINAL)} ligne(s))")
  } else {
    log_warn("Aucun resultat pour {current_target$name} (BIG_TABLE vide)")
  }
  
  stop_time = Sys.time()
  log_info("### Temps d'execution : {current_target$name} - {round(difftime(stop_time, start, units='mins'),2)} min ###")
  
} #end for loop

log_info("=== Fin MetIDfyR v2 ===")