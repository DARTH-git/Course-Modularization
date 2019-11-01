#### Wrap decision tree in a function ####
calculate_ce_out <- function(l_params_all, n_wtp=500){
  with(
    as.list(l_params_all),
    {
      #### Create vector of weights for each strategy ####
      # weight per path for the ondasentron arm
      v_w_on <-  c(     p_sem_on    *   (1 - p_sem_ade_on),                                                                            # emesis with no ADEs
                        p_sem_on    *        p_sem_ade_on        *       p_sem_ade_trt_on       *       p_sem_ade_trt_res_on,          # ADE after emesis that was treated and resolved
                        p_sem_on    *        p_sem_ade_on        *       p_sem_ade_trt_on       *  (1 - p_sem_ade_trt_res_on),         # ADE after emesis that was treated and did not resolve
                        p_sem_on    *        p_sem_ade_on        *  (1 - p_sem_ade_trt_on)      *       p_sem_ade_no_trt_res_on,       # ADE after emesis that was not treated and resolved
                        p_sem_on    *        p_sem_ade_on        *  (1 - p_sem_ade_trt_on)      *  (1 - p_sem_ade_no_trt_res_on),      # ADE after emesis that was not treated and did not resolve
                   (1 - p_sem_on)   *   (1 - p_no_sem_ade_on),                                                                         # no emesis with no ADEs
                   (1 - p_sem_on)   *        p_no_sem_ade_on     *       p_no_sem_ade_trt_on    *       p_no_sem_ade_trt_res_on,       # ADE after no emesis that was treated and resolved
                   (1 - p_sem_on)   *        p_no_sem_ade_on     *       p_no_sem_ade_trt_on    *  (1 - p_no_sem_ade_trt_res_on),      # ADE after no emesis that was treated and did not resolve
                   (1 - p_sem_on)   *        p_no_sem_ade_on     *  (1 - p_no_sem_ade_trt_on)   *       p_no_sem_ade_no_trt_res_on,    # ADE after no emesis that was not treated and resolved
                   (1 - p_sem_on)   *        p_no_sem_ade_on     *  (1 - p_no_sem_ade_trt_on)   *  (1 - p_no_sem_ade_no_trt_res_on))   # ADE after no emesis that was not treated and did not resolve
      
      # weight per path for the metoclopramide arm
      v_w_met <- c(     p_sem_met   *   (1 - p_sem_ade_met),                                                                           # emesis with no ADEs
                        p_sem_met   *        p_sem_ade_met       *       p_sem_ade_trt_met      *       p_sem_ade_trt_res_met,         # ADE after emesis that was treated and resolved
                        p_sem_met   *        p_sem_ade_met       *       p_sem_ade_trt_met      *  (1 - p_sem_ade_trt_res_met),        # ADE after emesis that was treated and did not resolve
                        p_sem_met   *        p_sem_ade_met       *  (1 - p_sem_ade_trt_met)     *       p_sem_ade_no_trt_res_met,      # ADE after emesis that was not treated and resolved
                        p_sem_met   *        p_sem_ade_met       *  (1 - p_sem_ade_trt_met)     *  (1 - p_sem_ade_no_trt_res_met),     # ADE after emesis that was not treated and did not resolve
                   (1 - p_sem_met)  *   (1 - p_no_sem_ade_met),                                                                        # no emesis with no ADEs
                   (1 - p_sem_met)  *        p_no_sem_ade_met    *       p_no_sem_ade_trt_met   *       p_no_sem_ade_trt_res_met,      # ADE after no emesis that was treated and resolved
                   (1 - p_sem_met)  *        p_no_sem_ade_met    *       p_no_sem_ade_trt_met   *  (1 - p_no_sem_ade_trt_res_met),     # ADE after no emesis that was treated and did not resolve
                   (1 - p_sem_met)  *        p_no_sem_ade_met    *  (1 - p_no_sem_ade_trt_met)  *       p_no_sem_ade_no_trt_res_met,   # ADE after no emesis that was not treated and resolved
                   (1 - p_sem_met)  *        p_no_sem_ade_met    *  (1 - p_no_sem_ade_trt_met)  *  (1 - p_no_sem_ade_no_trt_res_met))  # ADE after no emesis that was not treated and did not resolve
      
      
      #### Create vector of costs for each strategy ####
      # Estimating cost per path for the ondasentron arm
      v_c_on <-                 
        c(c_trt_on + c_sem_on,
          c_trt_on + c_sem_on + c_ade_on + c_ade_trt_on,
          c_trt_on + c_sem_on + c_ade_on + c_ade_trt_on, 
          c_trt_on + c_sem_on + c_ade_on,
          c_trt_on + c_sem_on + c_ade_on,
          c_trt_on,
          c_trt_on + c_ade_on + c_ade_trt_on,
          c_trt_on + c_ade_on + c_ade_trt_on, 
          c_trt_on + c_ade_on,
          c_trt_on + c_ade_on)  
      
      # Estimating cost per path for the metoclopramide arm
      v_c_met <-              
        c(c_trt_met + c_sem_met,
          c_trt_met + c_sem_met + c_ade_met + c_ade_trt_met,
          c_trt_met + c_sem_met + c_ade_met + c_ade_trt_met, 
          c_trt_met + c_sem_met + c_ade_met,
          c_trt_met + c_sem_met + c_ade_met,
          c_trt_met,
          c_trt_met + c_ade_met + c_ade_trt_met,
          c_trt_met + c_ade_met + c_ade_trt_met, 
          c_trt_met + c_ade_met,
          c_trt_met + c_ade_met) 
      
      #### Create vector of utilities for each strategy ####
      # vector of health outcomes for both therapies for each path
      v_e_on  <- v_e_met <- c(0, 
                              0,
                              0,
                              0, 
                              0, 
                              1,
                              0,
                              0, 
                              0,
                              0) 
      
      #### Calculate total costs for each strategy ####
      tc_on   <- v_w_on   %*%  v_c_on      
      tc_met  <- v_w_met  %*%  v_c_met
      
      #### Calculate total utilities for each strategy ####
      te_on   <- v_w_on   %*%  v_e_on      
      te_met  <- v_w_met  %*%  v_e_met
      
      v_tc    <- c(tc_on, tc_met)    # vector of total costs
      v_te    <- c(te_on, te_met)    # vector of total life years
      v_nmb   <- v_te * n_wtp - v_tc # calculate vector of nmb
      
      # Name outcomes
      names(v_tc)  <- v_names_str   # names for the elements of the  tc  vector
      names(v_te)  <- v_names_str   # names for the elements of the  te vector
      names(v_nmb) <- v_names_str   # names for the elements of the nmb vector
      
      df_output  <- data.frame(Strategy        =  v_names_str,
                               Cost            =  v_tc,
                               Effectiveness   =  v_te,
                               NMB             =  v_nmb)
     return(df_output)})
}

