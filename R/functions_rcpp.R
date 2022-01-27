# functions_rcpp.R
# main functions for running the DynaMICE model based on Rcpp files
# update: 2022/01/20

# ------------------------------------------------------------------------------
## Execute the Rcpp measles model for a country and a PSA run
# ------------------------------------------------------------------------------
runCountry_rcpp <- function (
  #variables specific for loop
  iso3,
  years,

  #infection dynamic variables
  vaccination,
  using_sia,
  parms_rcpp,

  #input data
  c_coverage_routine,
  c_coverage_sia,
  c_timeliness,
  contact_mat,
  c_contact,
  c_rnought,
  c_population,

  # dynaMICE model options
  tstep,
  save_scenario,
  foldername,
  log_name,

  # PSA options
  r,
  runs,
  psa,
  psa_var  # used only when psa > 0

) {

  # adjusted filename
  if (iso3 == "XK") {iso3 <- "XKX"}

  # stochastic parameters
  if (psa > 0) {
    # vaccine efficacy
    ve1_intcp <- as.numeric (psa_var [r, "take1_input"])
    ve1_slope <- as.numeric (psa_var [r, "take2_input"])
    ve2plus   <- as.numeric (psa_var [r, "take3_input"])
    age_ve1   <- ve1_intcp + ve1_slope * 12 * c(1:(3*52)/52, 4:101)  # based on age in months
    age_ve1   <- ifelse (age_ve1 >= ve2plus, ve2plus, age_ve1)
    parms_rcpp$ve1     <- age_ve1
    parms_rcpp$ve2plus <- ve2plus
    # R0
    c_rnought <- as.numeric (psa_var [r, "R0_input"])
    # distribution of SIA doses (portnoy = 1+0, nreach = 1+1)
    if (using_sia > 0) {
      using_sia <- 1 + as.numeric (psa_var [r, "siadose_input"] == "nreach")
    }
  }

  # country specific timeliness curve
  country_timeliness <- c_timeliness [!is.na(age), timeliness]
  timeliness_ages    <- c_timeliness [!is.na(age), age]

  # expand 0-2 years old to weekly age strata
  s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
  jt        <- 3  # how many ages to expand to s (or to weeks)
  beta_full <- matrix (0, ncol = 254, nrow = 254)

  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
    A = c_contact [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
    expand_rows =  s, expand_cols =  s,
    rescale_rows = FALSE, rescale_cols = FALSE)

  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
    A = c_contact [1:jt,(jt+1):ncol(c_contact)],
    expand_rows = s, expand_cols = 1,
    rescale_rows = F, rescale_cols = F)

  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
    A = c_contact [(jt+1):nrow(c_contact),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
    expand_rows = 1, expand_cols = s,
    rescale_rows = F, rescale_cols = F)

  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
    c_contact [(jt+1):nrow(c_contact),(jt+1):ncol(c_contact)]

  beta_full_unadj <- beta_full

  # infection rate under target R0
  beta_tstep  <- c_rnought*parms_rcpp$gamma

  # setup inputs for DynaMICE Rcpp model
  n_years   <- length(years)

  y_out       <- array(0, c(254, 14, n_years))   # numbers of age groups, compartments, years
  #case_out    <- array(0, c(254, n_years))
  case0d_out  <- array(0, c(254, n_years))
  case1d_out  <- array(0, c(254, n_years))
  case2d_out  <- array(0, c(254, n_years))
  pop_out     <- array(0, c(254, n_years))
  pop0d_out   <- array(0, c(254, n_years))
  popSus_out  <- array(0, c(254, n_years))
  dose_out    <- array(0, c(254, n_years))       # number of doses administrated
  reach0_out  <- array(0, c(254, n_years))       # number of zero-dose population reached
  fvp_out     <- array(0, c(254, n_years))       # number of fully vaccinated population reached

  init_Comp     <- matrix(0, 254, 14)
  init_Comp[,2] <- 0.95
  init_Comp[,3] <- 0.05

  t_spinup <- 1:1e5 # assume a fixed period for equilibrium period

  writelog (log_name, paste0(iso3, "; Run ",r,"/",runs,"; Generated input data & started model run"))

  # run model by yearly input
  for (y in years) {

    #old script groups those aged 70-80, but division is by actual population size
    pop.vector <- c_population[year == y, value]

    # first expand polymod matrix (contact_tstep) and population vector and
    # then divide by population sizes, otherwise it doesn't work.
    pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])

    # change zero values to 1 to avoid division by 0
    pop.vector_full[pop.vector_full==0] <- 1

    # generate contact matrices for DynaMICE using uniform age structure
    # adjust contact reciprocity by population structure of each year
    if (contact_mat == "unimix") {

      beta_full <- beta_tstep * beta_full_unadj

    } else if (contact_mat == "prpmix") {

      beta_full <- beta_tstep * matrix (rep(pop.vector_full/sum(pop.vector_full), each = 254), nrow = 254)

    } else {

      beta_full <- matrix(0, 254, 254)
      beta_full_R0 <- matrix(0, 254, 254)
      for (i in 1:254){
        for (j in 1:254) {
          beta_full [i, j] <- (beta_full_unadj[i, j] * pop.vector_full[i] +
                                 beta_full_unadj[j, i] * pop.vector_full[j])/(2*pop.vector_full[i])

          # calculate infection rate of the contact matrix based on the R0 definition
          # transform from "contactees per contactor" to "contactors per contactee"
          # This step can be skipped, since the largest eigenvalue remains unchanged.
          beta_full_R0 [i, j] <-  beta_full [i, j] * (pop.vector_full[i] / pop.vector_full[j])
        }
      }
      # make sure the contact matrix to represent target R0
      beta_full <- (beta_tstep / Re(eigen(beta_full_R0, only.values=T)$values[1])) * beta_full
    }

    # run spin-up period
    if (y == years[1]){
      out_Comp <- rcpp_spinup (init_Comp, parms_rcpp, beta_full, pop.vector_full, length(t_spinup))
      # print ('Spin-up period finished')
    }

    if (vaccination >= 1) {

      # Maximum coverage can (obviously) only be 100%
      # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
      # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
      # In essence, this becomes the inverse of the cumulative timeliness curve
      cycov <- c_coverage_routine [year == y & vaccine == "MCV1", coverage] /
        c_timeliness [is.na(age), prop_final_cov]

      if (length(cycov) == 0) {   # check if vaccine is not yet introduced and thereby, coverage value missing for this year
        cycov <- 0
      } else if (is.na(cycov)) {
        cycov <- 0
      }

      country_year_timeliness_mcv1 <- 1 - min(cycov,1) * country_timeliness

      country_year_timeliness_mcv1 <- -diff(country_year_timeliness_mcv1) /
        (country_year_timeliness_mcv1 [1:(length(country_year_timeliness_mcv1)-1)])

      # Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
      country_year_timeliness_mcv1 [is.na(country_year_timeliness_mcv1)] <- 0
      country_year_timeliness_mcv1 [is.nan(country_year_timeliness_mcv1)] <- 0
      country_year_timeliness_mcv1_allages <- rep(0, 254)
      country_year_timeliness_mcv1_allages [round(timeliness_ages)] <- country_year_timeliness_mcv1

    } else {
      country_year_timeliness_mcv1_allages <- rep(0, 254)
    }

    if(vaccination == 2){
      country_year_mcv2 <- c_coverage_routine [year == y & vaccine == "MCV2", coverage]
    } else {
      country_year_mcv2 <- 0
    }

    # if ( length (country_year_mcv2) == 0 || is.na(country_year_mcv2) ) {
    #   country_year_mcv2 <- 0
    # }


    # set up SIA inputs
    setorder(c_coverage_sia, year, start_m)
    if ((using_sia == 1) && (dim(c_coverage_sia [year == y])[1] > 0)) {
      # set up timesteps based on mid-month
      sia_mn <- c_coverage_sia [year == y, start_m]
      t_sia_m <- c()
      for (im in unique(sia_mn)){
        t_sia_m <- c(t_sia_m,
                     round((tstep/24) + (im-1)*(tstep/12)) + (1:sum(sia_mn == im)))
      }

      sia_input <- list (
        sia_implement = 2, # use 7.7% never-reach assumption
        a0 = c_coverage_sia [year == y, a0],
        a1 = c_coverage_sia [year == y, a1],
        sia_cov = c_coverage_sia [year == y, coverage],
        sia_mt = c(t_sia_m, 2000) # add a larger-than-timestep number for Rcpp model run
      )
    } else {
      sia_input <- list (
        sia_implement = as.integer(0),
        a0 = as.integer(0),
        a1 = as.integer(0),
        sia_cov = as.numeric(0),
        sia_mt = as.integer(0)
      )
    }

    t_start <- length(t_spinup) + (y-years[1])*tstep + 1
    outp <- rcpp_vaccine_oney (out_Comp,
                               parms_rcpp,
                               sia_input,
                               beta_full,
                               pop.vector_full,
                               country_year_timeliness_mcv1_allages,
                               country_year_mcv2,
                               t_start)
    out_Comp <- outp$out_Comp

    y_out    [, , (y-years[1])+1] <- out_Comp
    #case_out   [, (y-years[1])+1] <- outp$cases*pop.vector_full
    case0d_out [, (y-years[1])+1] <- outp$cases_0d*pop.vector_full     # new cases among 0-dose
    case1d_out [, (y-years[1])+1] <- outp$cases_1d*pop.vector_full     # new cases among 1-dose
    case2d_out [, (y-years[1])+1] <- outp$cases_2d*pop.vector_full     # new cases among >=2-dose
    pop_out    [, (y-years[1])+1] <- rowSums(out_Comp[, 1:13])*pop.vector_full         # all compartments
    pop0d_out  [, (y-years[1])+1] <- rowSums(out_Comp[, 1:4])*pop.vector_full          # sum of M, S, I, R
    popSus_out [, (y-years[1])+1] <- rowSums(out_Comp[, c(2,5,8,11)])*pop.vector_full  # sum of S, V1S, V2S, V3S
    dose_out   [, (y-years[1])+1] <- outp$doses*pop.vector_full
    reach0_out [, (y-years[1])+1] <- outp$reach_d0*pop.vector_full
    fvp_out    [, (y-years[1])+1] <- outp$fvps*pop.vector_full


    #if(y %% 20 == 0) {print (paste0 ('year ', y, ' finished'))}
  }

  saveRDS (list(#cases   = rbind (colSums(  case_out[1:52,]), colSums(  case_out[53:104,]), colSums(  case_out[105:156,]),
                cases0d  = rbind (colSums(case0d_out[1:52,]), colSums(case0d_out[53:104,]), colSums(case0d_out[105:156,]), case0d_out[157:254,]),
                cases1d  = rbind (colSums(case1d_out[1:52,]), colSums(case1d_out[53:104,]), colSums(case1d_out[105:156,]), case1d_out[157:254,]),
                cases2d  = rbind (colSums(case2d_out[1:52,]), colSums(case2d_out[53:104,]), colSums(case2d_out[105:156,]), case2d_out[157:254,]),
                pops     = rbind (colSums(   pop_out[1:52,]), colSums(   pop_out[53:104,]), colSums(   pop_out[105:156,]),    pop_out[157:254,]),
                pops0d   = rbind (colSums( pop0d_out[1:52,]), colSums( pop0d_out[53:104,]), colSums( pop0d_out[105:156,]),  pop0d_out[157:254,]),
                popsSus  = rbind (colSums(popSus_out[1:52,]), colSums(popSus_out[53:104,]), colSums(popSus_out[105:156,]), popSus_out[157:254,]),
                doses    = rbind (colSums(  dose_out[1:52,]), colSums(  dose_out[53:104,]), colSums(  dose_out[105:156,]),   dose_out[157:254,]),
                reachs0d = rbind (colSums(reach0_out[1:52,]), colSums(reach0_out[53:104,]), colSums(reach0_out[105:156,]), reach0_out[157:254,]),
                fvps     = rbind (colSums(   fvp_out[1:52,]), colSums(   fvp_out[53:104,]), colSums(   fvp_out[105:156,]),    fvp_out[157:254,])),
           file = ifelse (psa > 0,
                          paste0 ("outcome/", save_scenario, "/", foldername, "/", iso3, "_", r, ".RDS"),
                          paste0 ("outcome/", save_scenario, "/", foldername, "/", iso3, ".RDS"))
  )

  writelog (log_name, paste0 (iso3, "; Run ",r,"/",runs,"; Finished model run & saved outputs"))
}
# end of function -- runCountry_rcpp
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
## Run the Rcpp measles model for a selected vaccination strategy
# ------------------------------------------------------------------------------
runScenario_rcpp <- function (vaccine_coverage_folder    = "",
                              vaccine_coverage_subfolder = "",
                              coverage_prefix            = "",
                              antigen                    = "",
                              scenario_name,
                              save_scenario,
                              burden_estimate_folder,                  # burden estimate folder
                              group_name,                              # modelling group name
                              log_name,
                              countries                  = "all",
                              cluster_cores              = 1,
                              psa                        = 0,          # psa runs; 0 for single run
                              vaccination,  # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
                              using_sia,    # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA  (Portnoy), 2: with SIA (7.7%)
                              contact_mat,                             # contact matrix: "prpmix","polymod","syn"
                              sim_years                  = 1980:2100   # calendar years for simulation
) {
  # --------------------------------------------------------------------------
  # define global parameters
  # --------------------------------------------------------------------------
  ages 		    <- c (0:100)	              # Numeric vector with age-strata that are reported (Model ALWAYS models 101 age-groups; 0-2yo in weekly age-strata, and 3-100yo in annual age-strata)
  psa 		    <- psa		                  # Number of PSAs to run (0 to use deterministic).
  dinf		    <- 14				                # duration of infection (days)
  tstep			  <- 1000				              # number of timesteps in a year

  # age-dependent vaccine efficacy for first dose, based on a linear model (Hughes et al. 2020)
  ve1_intcp   <- 0.64598                  # intercept of the linear model
  ve1_slope   <- 0.01485                  # slope of the linear model, per month of age
  ve2plus     <- 0.98                     # vaccine efficacy for two and more doses

  age_ve1 <- ve1_intcp + ve1_slope * 12 * c(1:(3*52)/52, 4:101)  # based on age in months
  age_ve1 <- ifelse (age_ve1 >= ve2plus, ve2plus, age_ve1)

  # parameters for Rcpp functions
  parms_rcpp <- list (gamma         = 1 / (dinf * tstep/365),    # rate of losing infectivity
                      tstep         = tstep,
                      amp           = 0.05,		                   # amplitude for seasonality
                      ve1           = age_ve1,
                      ve2plus       = ve2plus)

  # number of clusters to use
  # if larger than 1, country-specific model runs are distributed over specified number of clusters
  # note that model uses a lot of memory, so might not want to max out all clusters
  use_cluster  <- cluster_cores

  # --------------------------------------------------------------------------
  # set up folders and parallel computation
  # --------------------------------------------------------------------------
  # load correct libraries when on cluster (using open MPI)
  if ("cluster_using_openmpi" %in% commandArgs()){
    using_openmpi <- TRUE
    # require(doMPI)
    require(parallel)
    cl <- doMPI::startMPIcluster()
    doMPI::registerDoMPI(cl)
    print("using openMPI")
  } else {
    using_openmpi <- FALSE
    print("no openMPI")
  }

  # determine the folder name: deterministic or stochastic?
  if (psa == 0) {
    det_stoch <- "deter" # deterministic
    runs      <- 1

  } else {
    det_stoch <- "stoch" # stochastic
    runs      <- psa
  }

  # create folders with correct name for in- and output data if not yet exists
  # typically foldername should not exist - but may come in handy when only processing results
  if ( !exists("foldername_analysis") ) {

    foldername <- paste0 (
      format(Sys.time(),format="%Y%m%d"),
      "_v",
      vaccination,
      "_s",
      using_sia,
      "_",
      det_stoch
    )

    dir.create(
      file.path(
        paste0(
          getwd(),
          "/outcome/", save_scenario, "/",
          foldername
        )
      ), recursive = T
    )

    if (psa > 0) {

      for(p in 1:psa){
        if(p<10){
          p <- paste0("00",p)
        } else if(p<100){
          p <- paste0("0",p)
        }
      }
    }
  } else {
    foldername <- foldername_analysis
  }

  # log
  writelog (log_name, paste0("Main; started"))

  if (using_openmpi) {
    writelog (log_name, paste0 ("Main; OpenMPI enabled"))
  } else {
    writelog (log_name, paste0 ("Main; OpenMPI disabled"))
  }

  # --------------------------------------------------------------------------
  # prepare input data
  # --------------------------------------------------------------------------
  # define filename of coverage data
  data_coverage_routine <- paste0 (vaccine_coverage_folder,
                                   vaccine_coverage_subfolder,
                                   "routine_",
                                   scenario_name,
                                   ".csv")

  data_coverage_sia <- paste0 (vaccine_coverage_folder,
                               vaccine_coverage_subfolder,
                               "sia_",
                               scenario_name,
                               ".csv")

  # import/read data
  coverage_routine	<- copy (fread (data_coverage_routine))[year %in% sim_years]
  coverage_sia		  <- copy (fread (data_coverage_sia))[year %in% sim_years]
  timeliness  		  <- setDT (data_timeliness)
  rnought	    		  <- setDT (data_r0)
  population  		  <- setDT (data_pop)
  template    		  <- setDT (data_template)

  # if psa variables file does not exist, then create psa variables file
  if (psa > 0) {

    # data file name for PSA variables using create_PSA_data()
    data_psa 		<- "psa_variables.csv"

    if (file.exists (data_psa)) {

      # read csv if file already exists
      psa_var <- fread (data_psa)

      # check if the number of PSA runs in the parameter file corresponds with psa
      if(nrow (psa_var) != psa){
        stop(paste0 ("Number of runs in PSA file is not the same as those specified! Variables used in the same run in each scenario should be similar. Delete or rename ",
                     data_psa,
                     " if a new file needs to be created, or change the number of PSA runs to ",
                     nrow(psa_var),
                     " if the same file should be used."))
      }
    } else {
      stop (paste0 ("There is no files for inputting PSA variables. Use create_PSA_data() to generate the file."))
    }
  }

  # if countries are specified to all, then set countries to all countries in coverage file
  if (countries[[1]] == "all") {
    countries	<- as.character (unique (coverage_routine [, country_code] ) )
  }

  # Process contact matrices
  contact_list <- switch (contact_mat,
                          "syn"     = sapply (countries,
                                              function(cty){data_contact_syn[[cty]]},
                                              simplify = FALSE, USE.NAMES = TRUE),
                          "polymod" = sapply (countries,
                                              function(cty){data_contact_polymod[[cty]]},
                                              simplify = FALSE, USE.NAMES = TRUE),
                          "prpmix"  = sapply (countries,
                                              function(x = NULL){matrix (1/101, nrow = 101, ncol = 101)},
                                              simplify = FALSE, USE.NAMES = TRUE), # not actually used but to fit in the model structure
                          "unimix"  = sapply (countries,
                                              function(x = NULL){matrix (1/101, nrow = 101, ncol = 101)},
                                              simplify = FALSE, USE.NAMES = TRUE)
  )


  # ----------------------------------------------------------------------------
  # Run model
  # ----------------------------------------------------------------------------
  writelog (log_name, paste0 ("Main; Foldername: ", foldername))

  if (using_openmpi | use_cluster > 1){
    if (!using_openmpi){

      if (use_cluster != parallel::detectCores()) {
        warning (paste0(parallel::detectCores(), " cores detected but ", use_cluster, " specified."))
      }

      cl <- parallel::makeCluster (use_cluster)
      doParallel::registerDoParallel (cl)

    } else {
      print (paste0 ("Clustersize: ", doMPI::clusterSize(cl)))
    }
  }

  # foreach will run countries and PSA runs in parallel if a parallel backend
  # is registered, and sequentially otherwise   #, "rcpp_spinup","rcpp_vaccine_oney"
  # require(foreach)
  combine <-
    foreach (
      ii = 1:length(countries),
      .packages = c("data.table", "Rcpp"), #, "dynamiceRcpp"
      .errorhandling = "stop", # "pass
      .export = c("runCountry_rcpp", "rcpp_spinup", "rcpp_vaccine_oney",
                  "writelog", "expandMatrix")) %:%   # "updateProgress"
    foreach (
      r = 1:runs,
      .errorhandling = "stop") %dopar% {
        iso3    <- countries[ii]
        out_run <- runCountry_rcpp (iso3               = iso3,
                                    years              = as.numeric (sim_years),
                                    vaccination        = vaccination,
                                    using_sia          = using_sia,
                                    parms_rcpp         = parms_rcpp,
                                    c_coverage_routine = coverage_routine[country_code == iso3,],
                                    c_coverage_sia     = coverage_sia[country_code == iso3 & coverage != 0,],
                                    c_timeliness       = timeliness[country_code == iso3,],
                                    contact_mat        = contact_mat,
                                    c_contact          = contact_list[[iso3]],
                                    c_rnought          = rnought[country_code == iso3, r0],
                                    c_population       = population[country_code == iso3,],
                                    tstep              = tstep,
                                    save_scenario      = save_scenario,
                                    foldername         = foldername,
                                    log_name           = log_name,
                                    r                  = r,
                                    runs               = runs,
                                    psa                = psa,
                                    psa_var            = psa_var
        )
        return(out_run)
      }

  if(using_openmpi | use_cluster > 1){
    if (!using_openmpi){
      parallel::stopCluster(cl)
    } else {
      doMPI::closeCluster(cl)
    }
  }

  # check for errors
  errorcount <- 0
  for (i in 1:length(combine)){
    if ("error" %in% class(combine[[i]])){
      errormessage <- paste0 ("Error in task ", i, ": ", combine[[i]])
      warning(errormessage)
      writelog (log_name, errormessage)
      errorcount <- errorcount + 1
      #remove from data
      combine[[i]] <- NULL
    }
  }
  if(errorcount > 0){
  stop(paste0("There were ",errorcount," errors."))
  }

  if (using_openmpi) {
    Rmpi::mpi.quit()
  }
  return()

} # end of function -- runScenario_rcpp
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
## Merge the output of cases
# ------------------------------------------------------------------------------
merge_case_csv <- function (vaccine_coverage_folder,
                            vaccine_coverage_subfolder,
                            antigen,
                            scenario_name,
                            save_scenario,
                            burden_estimate_folder,
                            group_name,
                            log_name,
                            psa = 0,
                            vaccination,
                            using_sia,
                            folder_date,
                            sim_years,
                            chunksize
) {
  # ----------------------------------------------------------------------------
  # merge and process results
  # ----------------------------------------------------------------------------
  # define folder name
  foldername <- paste0 (
    folder_date,
    "_v", vaccination,
    "_s", using_sia,
    "_", ifelse (psa == 0, "deter", "stoch"))

  # merge results
  output_files <- list.files (path = paste0 ("outcome/", save_scenario, "/",foldername,"/"),
                              recursive = T, full.names = T)
  nfile        <- length(output_files)
  chunckdist   <- split (1:nfile, rep(1:chunksize, each = nfile/chunksize))

  # set up format of output file
  years          <- as.numeric (sim_years)
  ages           <- c(0:100)
  template    	 <- copy (setDT (data_template)) [year %in% years]   # drop rows if simulation period is shorter
  report_years   <- sort (unique (template$year))
  country_names  <- unique (subset (template, select = c("country", "country_name")))
  c_names        <- country_names$country_name
  names(c_names) <- country_names$country

  # file name for burden estimates
  burden_estimate_file <- paste0 (ifelse (psa == 0, "central", "stochastic"),
                                  "_burden_estimate_",
                                  antigen,
                                  group_name,
                                  scenario_name)

  # coverage file
  coverage_routine	<- copy (fread (paste0 (vaccine_coverage_folder,
                                           vaccine_coverage_subfolder,
                                           "routine_",
                                           scenario_name,
                                           ".csv"))
  )[year %in% sim_years]

  for (iprocess in 1:chunksize) {
    # read RDS files
    all_runs <- rbindlist (lapply (output_files[chunckdist[[iprocess]]], function (filename, ...) {
      res <- withCallingHandlers (
        readRDS (filename),
        warning = function(w) {warning(w, filename);}
      )
      filefinal <- stringr::str_extract (filename, "[^/]+$")      # remove path but keep filename
      res2 <- data.table (country     = rep (stringr::str_sub(filefinal, 1, 3), length(ages)*length(years)), # 101 ages * 121 years
                          year        = rep (years, each = length(ages)),
                          age         = rep (ages, length(years)),
                          #cases       = as.vector(res$cases),
                          cases0d     = as.vector(res$cases0d),
                          cases1d     = as.vector(res$cases1d),
                          cases2d     = as.vector(res$cases2d),
                          pops        = as.vector(res$pops),
                          pops0d      = as.vector(res$pops0d),
                          popsSus     = as.vector(res$popsSus),
                          doses       = as.vector(res$doses),
                          reachs0d    = as.vector(res$reachs0d),
                          fvps        = as.vector(res$fvps)
      )

      res2 [country == "XKX", country := "XK"]
      if (psa > 0) {
        res2 [, run_id := stringr::str_extract (filefinal, "\\d+")]
      }
      return (res2)
    }))

    # add country names and disease (Measles) to match template file
    all_runs[, c("country_name", "disease") := list (c_names[country], "Measles")]

    # select output years
    all_runs <- subset (all_runs, year %in% report_years)

    # --------------------------------------------------------------------------
    # add columns for remaining life expectancy & MCV1
    # --------------------------------------------------------------------------
    sel_countries <- unique (all_runs$country)

    # remaining life expectancy
    lexp_remain <- tailor_data_lexp_remain (sel_countries)

    if (psa > 0) {
      all_runs <- lexp_remain [all_runs,
                               .(i.country, run_id, year, age, cases0d, cases1d, cases2d,
                                 pops, pops0d, popsSus, doses, reachs0d, fvps,
                                 country_name, disease, value),
                               on = .(country_code = country,
                                      age          = age,
                                      year         = year)]
    } else {
      all_runs <- lexp_remain [all_runs,
                               .(i.country, year, age, cases0d, cases1d, cases2d,
                                 pops, pops0d, popsSus,  doses, reachs0d, fvps,
                                 country_name, disease, value),
                               on = .(country_code = country,
                                      age          = age,
                                      year         = year)]
    }

    # rename column names for output
    setnames (all_runs,
              old = c("i.country", "value"      ),
              new = c("country"  , "remain_lexp"))


    # MCV1 coverage
    coverage_routine_MCV1 <- coverage_routine [(vaccine == "MCV1") & (country_code %in% sel_countries)]

    if (psa > 0) {
      all_runs <- coverage_routine_MCV1 [all_runs,
                                         .(i.country, run_id, i.year, age,
                                           cases0d, cases1d, cases2d,
                                           pops, pops0d, popsSus,
                                           doses, reachs0d, fvps, country_name,
                                           disease, coverage, remain_lexp),
                                         on = .(country_code = country,
                                                year         = year) ]
    } else {
      all_runs <- coverage_routine_MCV1 [all_runs,
                                         .(i.country, i.year, age,
                                           cases0d, cases1d, cases2d,
                                           pops, pops0d, popsSus,
                                           doses, reachs0d, fvps, country_name,
                                           disease, coverage, remain_lexp),
                                         on = .(country_code = country,
                                                year         = year) ]
    }

    # rename column "coverage" to "MCV1"
    setnames (x = all_runs,
              old = c("i.country", "i.year", "coverage"),
              new = c("country"  , "year"  , "MCV1"    ))

    # save burden estimates to file
    suffix <- ifelse (chunksize == 1, ".csv", paste0 ("_", iprocess, ".csv"))
    fwrite (x     = all_runs [order (country, year, age)],
            file  = paste0 (burden_estimate_folder,  burden_estimate_file, suffix))

    # calculate mean estimates by country
    if (psa > 0){

      if ((dim(all_runs)[1] / (length(sel_countries)*length(ages)*length(report_years))) != psa) {
        stop ("Full PSA runs for the selected countries are not completely included. Check the RDS files.")
      } else {
        mean_runs <- all_runs [ , lapply(.SD, mean),
                                by = c("disease", "year", "age", "country", "country_name", "MCV1", "remain_lexp"),
                                .SDcols = pops:fvps]

        fwrite (x     = mean_runs [order (country, year, age)],
                file  = paste0 (burden_estimate_folder,
                                "mean_burden_estimate_",
                                antigen,
                                group_name,
                                scenario_name,
                                ".csv"),
                append = TRUE)
      }
    }

    # release memory for the next round
    remove (list = c("all_runs", "sel_countries"))

    writelog (log_name, paste0 ("case csv generated; ", iprocess, "/", chunksize))
  }

  return ()
} # end of function -- merge_case_csv
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## Estimate deaths and DALYs of measles - 2021 Portnoy's CFRs
# ------------------------------------------------------------------------------
estimate_deaths_dalys_21 <- function (antigen,
                                      group_name,
                                      scenario_name,
                                      log_name,
                                      burden_estimate_folder,
                                      psa = 0,
                                      chunksize
) {
  burden_estimate_file <- paste0 (ifelse (psa == 0, "central", "stochastic"),
                                  "_burden_estimate_",
                                  antigen,
                                  group_name,
                                  scenario_name)

  for (iprocess in 1:chunksize){
    # read case estimates
    suffix <- ifelse (chunksize == 1, ".csv",  paste0("_", iprocess, ".csv"))
    burden <- fread (paste0 (burden_estimate_folder, burden_estimate_file, suffix))

    data_cfr_21 <- setDT (copy (data_cfr_portnoy_21))

    # The data contains CFR estimates between 1981 and 2020, for age between 0 to 99.
    # Extrapolate CFRs for year 1980 and 2021-2100 and for age 100.

    min_year = min (burden [, year])
    max_year = max (burden [, year])
    chunk_countries = unique (burden$country)
    data_cfr_21 <- data_cfr_21 [year %in% min_year:max_year & country %in% chunk_countries]

    if (min_year < 1981) {
      data_cfr_add <- rbindlist (lapply (min_year:1981, function(i) copy (data_cfr_21 [year == 1981, ])[, year := i]))
      data_cfr_21  <- rbind     (data_cfr_add, data_cfr_21, use.names = TRUE)
    }

    if (max_year > 2020) {
      data_cfr_add <- rbindlist (lapply (2021:max_year, function(i) copy (data_cfr_21 [year == 2020, ])[, year := i]))
      data_cfr_21  <- rbind     (data_cfr_21, data_cfr_add, use.names = TRUE)
    }

    data_cfr_21  <- rbind (data_cfr_21,
                           copy (data_cfr_21 [age == 99, ])[, age := 100],
                           use.names = TRUE)
    setorder (data_cfr_21, country, year, age)

    if (psa > 0) {
      burden <- data_cfr_21 [burden,
                             .(disease, run_id, year, age, country, country_name,
                               pops, pops0d, popsSus, cases0d, cases1d, cases2d,
                               doses, reachs0d, fvps, cfr, remain_lexp),
                             on = .(country = country,
                                    year    = year,
                                    age     = age)]

      psa_var	<- fread ("psa_variables.csv")
      burden <- burden [psa_var [ , .(cfrpast_input, cfrproj_input, run_id)],
                        on = .(run_id = run_id)]
      burden [year <= 2015 , cfr := cfr + cfrpast_input]
      burden [year >  2015 , cfr := cfr + cfrproj_input]
      burden [cfr < 0, cfr := 0]

    } else {
      burden <- data_cfr_21 [burden,
                             .(disease, year, age, country, country_name,
                               pops, pops0d, popsSus,
                               cases0d, cases1d, cases2d, doses,
                               reachs0d, fvps, cfr, remain_lexp),
                             on = .(country = country,
                                    year    = year,
                                    age     = age)]
    }

    # estimate deaths
    burden [, `:=` (deaths0d = cases0d * cfr,
                    deaths1d = cases1d * cfr,
                    deaths2d = cases2d * cfr)]

    # calculate DALYs = (YLDs) + (YLLs)
    burden [, dalys := (((cases0d+cases1d+cases2d) - (deaths0d+deaths1d+deaths2d)) * 0.002) +
              ((deaths0d+deaths1d+deaths2d) * remain_lexp)]

    # adjust columns for output
    select.cols <- c("disease", "country", "country_name", "year", "age",
                     "pops", "pops0d", "popsSus",
                     "cases0d", "cases1d", "cases2d", "deaths0d", "deaths1d", "deaths2d",
                     "dalys", "doses", "reachs0d", "fvps")
    if (psa > 0) {
      save.cols <- c("run_id", select.cols)
    } else {
      save.cols <- select.cols
    }
    burden <- subset (burden, select = save.cols)

    # save updated burden estimate file (cases + deaths) to file
    # cfr_option is also the name of the subfolder
    fwrite (x      = burden,
            file   = paste0 (burden_estimate_folder,
                             "Portnoy/",
                             burden_estimate_file,
                             ".csv"),
            append = (iprocess != 1)
    )

    writelog (log_name, paste0 ("deaths csv generated; ",
                                iprocess, "/", chunksize))

    # calculate mean estimates by country
    if (psa > 0){
      central_burden <- burden [ , lapply(.SD, mean),
                                 by = select.cols[1:5],
                                 .SDcols = select.cols[6:16]]

      fwrite (x     = central_burden [order (country, year, age)],
              file  = paste0 (burden_estimate_folder,
                              "central_burden_estimate_",
                              antigen,
                              group_name,
                              scenario_name,
                              ".csv"),
              append = TRUE)
    }

    # release memory for the next round
    remove (list = c("burden", "data_cfr_21"))

  }
  return ()
} # end of function -- estimate_deaths_dalys_21
# ------------------------------------------------------------------------------
