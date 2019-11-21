#### Wrap decision tree in a function ####
calculate_ce_out <- function(l_params_all, n_wtp = 100000){
  with(
    as.list(l_params_all),
    {
      #### Create vector of weights for each strategy ####
      
      v_w_no_tx  <- c(  p_HVE  *    p_HVE_comp     ,  # HVE, complications
                        p_HVE  * (1-p_HVE_comp)    ,  # HVE, no complications
                     (1-p_HVE) *    p_OVE_comp     ,  # OVE, complications
                     (1-p_HVE) * (1-p_OVE_comp))      # OVE, no complications
      
      v_w_tx     <- c(  1                          ,  # On treatment
                        p_HVE  *    p_HVE_comp_tx  ,  # HVE w/tx, complications
                        p_HVE  * (1-p_HVE_comp_tx) ,  # HVE w/tx, no complications
                     (1-p_HVE) *    p_OVE_comp_tx  ,  # OVE w/tx, complications
                     (1-p_HVE) * (1-p_OVE_comp_tx))   # OVE w/tx, no complications
      
      v_w_biopsy <- c(  1                                              ,  # Undergo biopsy
                        p_biopsy_comp                                  ,  # biopsy complications
                     (1-p_biopsy_comp) *    p_HVE  *    p_HVE_comp_tx  ,  # no biopsy comp., HVE w/tx, complications 
                     (1-p_biopsy_comp) *    p_HVE  * (1-p_HVE_comp_tx) ,  # no biopsy comp., HVE w/tx, no complications
                     (1-p_biopsy_comp) * (1-p_HVE) *    p_OVE_comp     ,  # no biopsy comp., OVE, complications
                     (1-p_biopsy_comp) * (1-p_HVE) * (1-p_OVE_comp))      # no biopsy comp., OVE, no complications
      
      #### Create vector of outcomes (QALYs) for each strategy ####
      
      v_qaly_no_tx  <- c(q_VE_comp ,  # HVE, complications
                         q_VE      ,  # HVE, no complications
                         q_VE_comp ,  # OVE, complications
                         q_VE)        # OVE, no complications
      
      v_qaly_tx     <- c(0         ,  # treatment does not directly add any QALYs 
                         q_VE_comp ,  # HVE, complications
                         q_VE      ,  # HVE, no complications
                         q_VE_comp ,  # OVE, complications
                         q_VE)        # OVE, no complications
      
      
      v_qaly_biopsy <- c(q_loss_biopsy  ,  # loss due to biopsy
                         q_VE_comp      ,  # biopsy complications
                         q_VE_comp      ,  # no biopsy comp., HVE w/tx, complications 
                         q_VE           ,  # no biopsy comp., HVE w/tx, no complications
                         q_VE_comp      ,  # no biopsy comp., OVE, complications
                         q_VE)             # no biopsy comp., OVE, no complications
      
      #### Create vector of costs for each strategy ####
      
      v_cost_no_tx  <- c(c_VE_comp ,  # HVE, complications
                         c_VE      ,  # HVE, no complications
                         c_VE_comp ,  # OVE, complications
                         c_VE)        # OVE, no complications
      
      v_cost_tx     <- c(c_tx      ,  # cost of treatment
                         c_VE_comp ,  # HVE, complications
                         c_VE      ,  # HVE, no complications
                         c_VE_comp ,  # OVE, complications
                         c_VE)        # OVE, no complications
      
      
      v_cost_biopsy <- c(c_biopsy         ,  # cost of biopsy procedure
                         c_VE_comp        ,  # biopsy complications
                         c_VE_comp + c_tx ,  # no biopsy comp., HVE w/tx, complications 
                         c_VE + c_tx      ,  # no biopsy comp., HVE w/tx, no complications
                         c_VE_comp        ,  # no biopsy comp., OVE, complications
                         c_VE)               # no biopsy comp., OVE, no complications
      
      #### Calculate total utilities for each strategy ####
      total_qaly_no_tx  <- v_w_no_tx  %*%  v_qaly_no_tx      
      total_qaly_tx     <- v_w_tx     %*%  v_qaly_tx
      total_qaly_biopsy <- v_w_biopsy %*%  v_qaly_biopsy
      
      #### Calculate total costs for each strategy####
      total_cost_no_tx  <- v_w_no_tx  %*%  v_cost_no_tx    
      total_cost_tx     <- v_w_tx     %*%  v_cost_tx
      total_cost_biopsy <- v_w_biopsy %*%  v_cost_biopsy
      
      v_total_qaly <- c(total_qaly_no_tx, total_qaly_tx, total_qaly_biopsy)  # vector of total QALYs
      v_total_cost <- c(total_cost_no_tx, total_cost_tx, total_cost_biopsy)  # vector of total costs
      v_nmb        <- v_total_qaly * n_wtp - v_total_cost                    # calculate vector of nmb
      
      # Name outcomes
      names(v_total_qaly) <- v_names_str  # names for the elements of the total QALYs vector
      names(v_total_cost) <- v_names_str  # names for the elements of the total cost vector
      names(v_nmb)        <- v_names_str  # names for the elements of the nmb vector
      
      df_output <- data.frame(Strategy =  v_names_str,
                              Cost     =  v_total_cost,
                              Effect   =  v_total_qaly,
                              NMB      =  v_nmb)
     return(df_output)})
}

