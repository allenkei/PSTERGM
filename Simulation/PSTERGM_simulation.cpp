#include <RcppArmadillo.h>
#include <math.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

//////////////////////////
// Network Construction //
//////////////////////////

//[[Rcpp::export]]
List aug_dim_network(arma::mat network_before, arma::mat network_after){
  int n = network_before.n_rows;
  arma::mat network_augmentation = network_before;
  arma::mat network_diminution = network_before;
  List output;
  
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      network_augmentation(i,j) = max(network_after(i,j),network_before(i,j));
      network_diminution(i,j) = min(network_after(i,j),network_before(i,j));
    }
  }
  
  output.push_back(network_augmentation);
  output.push_back(network_diminution);
  
  return output;
  
}

///////////////////////////
// Proposal Distribution //
///////////////////////////

//[[Rcpp::export]]
int rZIpoi(double lambda, double pi0){
  double u;
  int y;
  
  u = unif_rand();
  if(u < pi0){
    y = 0;
  }else{
    y = R::rpois(lambda+0.5);
  }
  
  return y;
}

//[[Rcpp::export]]
double dZIpoi(int y, double lambda, double pi0){
  double prob;
  
  if(y == 0){
    prob = pi0 + (1-pi0) * exp(-(lambda+0.5));
  }else{
    prob = (1-pi0) * R::dpois(y, lambda+0.5, false);
  }
  
  return prob;
}

///////////////////////////////
// Generate Network Features //
///////////////////////////////

//[[Rcpp::export]]
arma::vec gen_feature_valued(arma::mat y, arma::vec node_attr){
  int n = y.n_rows;
  arma::vec holder(n);
  
  double network_sum = 0; 
  double num_zeros = 0;
  double mutual_geom = 0;
  double triad = 0;
  
  for(int i = 0; i < n; ++i){
    for(int j = 0; j < n; ++j){
      if(i != j){
        
        network_sum += y(i,j);
        
        if(y(i,j) == 0){num_zeros += 1;}
        
        if(i < j){mutual_geom += sqrt(y(i,j)*y(j,i));}
        
        holder.zeros(n);
        for(int k = 0; k < n; ++k){
          holder.at(k) = min(y(i,k),y(k,j));
        }
        triad += min(y(i,j),max(holder));
        
      }
    }
  }
  
  arma::vec output = {network_sum, num_zeros, mutual_geom, triad};
  
  return output;
}

//[[Rcpp::export]]
arma::vec change_stats(arma::mat y, arma::vec node_attr, double yij_new, int i, int j){
  int n = y.n_rows;
  
  arma::vec holder(n);
  arma::vec holder_new(n);
  arma::mat y_new = y;
  y_new(i,j) = yij_new;
  
  double d_network_sum = 0; 
  double d_num_zeros = 0;
  double d_mutual_geom = 0;
  double d_triad = 0;
  
  d_network_sum += yij_new - y(i,j);
  
  if(yij_new == 0 & y(i,j) > 0){d_num_zeros = 1;}
  if(yij_new > 0 & y(i,j) == 0){d_num_zeros = -1;}
  
  d_mutual_geom += sqrt(yij_new*y(j,i)) - sqrt(y(i,j)*y(j,i));
  
  holder.zeros(n);
  holder_new.zeros(n);
  for(int k = 0; k < n; ++k){
    holder.at(k) = min(y(i,k),y(k,j));
    holder_new.at(k) = min(y_new(i,k),y_new(k,j));
  }
  d_triad -= min(y_new(i,j),max(holder_new)) - min(y(i,j),max(holder)); // subtract
  
  for(int a = 0; a < n; ++a){ // fixed j
    holder.zeros(n);
    holder_new.zeros(n);
    for(int k = 0; k < n; ++k){
      holder.at(k) = min(y(a,k),y(k,j));
      holder_new.at(k) = min(y_new(a,k),y_new(k,j));
    }
    d_triad += min(y_new(a,j),max(holder_new)) - min(y(a,j),max(holder)); // add
  }
  
  for(int b = 0; b < n; ++b){ // fixed i
    holder.zeros(n);
    holder_new.zeros(n);
    for(int k = 0; k < n; ++k){
      holder.at(k) = min(y(i,k),y(k,b));
      holder_new.at(k) = min(y_new(i,k),y_new(k,b));
    }
    d_triad += min(y_new(i,b),max(holder_new)) - min(y(i,b),max(holder)); // add
  }
  
  arma::vec output = {d_network_sum, d_num_zeros, d_mutual_geom, d_triad};
  
  return output;
}

///////////////////////////////
// Generate MCMC sample (CD) //
///////////////////////////////

//[[Rcpp::export]]
List gen_MCMC_valued_CD(int nsim, int n_CD, arma::mat y_condition, arma::mat y_obs, double pi0, int m, 
                        arma::vec eta_atleast, arma::vec eta_atmost, arma::vec node_attr){
  int flip_counter;
  int yij_new;
  double q, r, u;
  int x, y;
  arma::mat y_syn, y_aug, y_dim, y_aug_temp, y_dim_temp;
  arma::vec delta_g_atleast, delta_g_atmost;
  List holder, aug_dim_holder;
  
  int n = y_condition.n_rows;
  
  for(int MCMC_iter = 0; MCMC_iter < nsim; ++MCMC_iter){
    
    flip_counter = 0;
    y_syn = y_obs; // initialize as obs
    
    aug_dim_holder = aug_dim_network(y_condition, y_syn);
    y_aug = as<arma::mat>(aug_dim_holder[0]);
    y_dim = as<arma::mat>(aug_dim_holder[1]);
    
    while(flip_counter < n_CD){
      
      x = n * unif_rand();y = n * unif_rand();
      while(y == x){x = n * unif_rand();y = n * unif_rand();}
      
      yij_new = rZIpoi(y_syn(x,y), pi0);
      q = dZIpoi(y_syn(x,y), yij_new, pi0) / dZIpoi(yij_new, y_syn(x,y), pi0);
      
      if(yij_new <= y_condition(x,y)){

        delta_g_atleast = change_stats(y_aug, node_attr, y_condition(x,y), x, y);
        delta_g_atmost = change_stats(y_dim, node_attr, yij_new, x, y);

        r = q * tgamma(y_aug(x,y)+1)/tgamma(y_condition(x,y)+1) * exp(dot(eta_atleast,delta_g_atleast)) *
          Rf_choose(m,yij_new)/Rf_choose(m,y_dim(x,y)) * exp(dot(eta_atmost,delta_g_atmost));

      }else{

        delta_g_atleast = change_stats(y_aug, node_attr, yij_new, x, y);
        delta_g_atmost = change_stats(y_dim, node_attr, y_condition(x,y), x, y);

        r = q * tgamma(y_aug(x,y)+1)/tgamma(yij_new+1) * exp(dot(eta_atleast,delta_g_atleast)) *
          Rf_choose(m,y_condition(x,y))/Rf_choose(m,y_dim(x,y)) * exp(dot(eta_atmost,delta_g_atmost));

      }
      
      u = unif_rand();
      if(u<r){
        // not symmetry
        y_syn(x,y) = yij_new;
        
        if(yij_new <= y_condition(x,y)){
          y_aug(x,y) = y_condition(x,y);
          y_dim(x,y) = yij_new;
        }else{
          y_aug(x,y) = yij_new;
          y_dim(x,y) = y_condition(x,y);
        }
        
      }
      
      flip_counter += 1;
      
    }
    
    holder.push_back(y_syn);
    
  }
  
  return holder;
  
}

/////////////////////////////
// Generate sample feature //
/////////////////////////////

//[[Rcpp::export]]
arma::mat gen_features_MCMC_valued(List MCMC_list, List y_list, int iter, arma::vec node_attr, int num_feature){
  int n = MCMC_list.size();
  List aug_dim_holder;
  arma::mat output; output.zeros(n,num_feature);
  arma::mat yt_1 = y_list[iter-1]; // previous network from the list of observed networks
  
  for(int gen_f_iter = 0; gen_f_iter < n; ++gen_f_iter){
    aug_dim_holder = aug_dim_network(yt_1, MCMC_list[gen_f_iter]); // MCMC_list a list of simulated y_t
    output.submat(gen_f_iter,0,gen_f_iter,num_feature/2-1) = gen_feature_valued(aug_dim_holder[0], node_attr).t();
    output.submat(gen_f_iter,num_feature/2,gen_f_iter,num_feature-1) = gen_feature_valued(aug_dim_holder[1], node_attr).t();
  }
  
  return output;
}

/////////////////////////////////
// Parameter Learning Temporal //
/////////////////////////////////

//[[Rcpp::export]]
arma::vec partial_stepping_temporal(double learning_iter, int MCMC_size, int n_CD, List y_list, arma::vec eta, arma::vec node_attr, double pi0, int m){
  
  List MCMC_list, aug_dim_holder;
  arma::vec eta_atleast, eta_atmost, g_obs_sum, mu_sum;
  arma::mat cov_sum;
  arma::mat xi_hat, MCMC_list_features;
  double par_gamma;
  int num_feature = eta.size();
  
  //features of observed networks
  g_obs_sum.zeros(num_feature);
  for(int iter = 1; iter < y_list.size(); ++iter){
    aug_dim_holder = aug_dim_network(y_list[iter-1], y_list[iter]);
    g_obs_sum.subvec(0, num_feature/2-1) += gen_feature_valued(aug_dim_holder[0],node_attr);
    g_obs_sum.subvec(num_feature/2,num_feature-1) += gen_feature_valued(aug_dim_holder[1],node_attr);
  }
  
  for(int l_iter = 0; l_iter < learning_iter; ++l_iter){
    
    //features of sampled networks
    mu_sum.zeros(num_feature);
    cov_sum.zeros(num_feature,num_feature);
    for(int iter = 1; iter < y_list.size(); ++iter){
      eta_atleast = eta.subvec(0, num_feature/2-1);
      eta_atmost = eta.subvec(num_feature/2,num_feature-1);
      MCMC_list = gen_MCMC_valued_CD(MCMC_size, n_CD, y_list[iter-1], y_list[iter], pi0, m, eta_atleast, eta_atmost, node_attr);
      MCMC_list_features = gen_features_MCMC_valued(MCMC_list, y_list, iter, node_attr, num_feature);
      mu_sum += arma::mean(MCMC_list_features,0).t();
      cov_sum += cov(MCMC_list_features);
    }
    
    par_gamma = (l_iter+1)/learning_iter;
    xi_hat = par_gamma*g_obs_sum + (1-par_gamma)*mu_sum;
    eta = eta + arma::solve(cov_sum,arma::eye<arma::mat>(num_feature,num_feature)) * (xi_hat-mu_sum);
    
    Rcout << "learning iter = " << l_iter+1 << ", eta (partial stepping): \n" << eta << "\n";
  }
  
  
  return eta;
  
}

//[[Rcpp::export]]
arma::vec newton_raphson_temporal(int learning_iter, int MCMC_size, int n_CD, List y_list, arma::vec eta, arma::vec node_attr, double pi0, int m){
  
  List MCMC_list, aug_dim_holder;
  arma::vec g_obs_sum, mu_sum, S, eta_atleast, eta_atmost;
  arma::mat cov_sum, H;
  arma::mat MCMC_list_features;
  int num_feature = eta.size();
  
  //features of observed networks
  g_obs_sum.zeros(num_feature);
  for(int iter = 1; iter < y_list.size(); ++iter){
    aug_dim_holder = aug_dim_network(y_list[iter-1], y_list[iter]);
    g_obs_sum.subvec(0, num_feature/2-1) += gen_feature_valued(aug_dim_holder[0],node_attr);
    g_obs_sum.subvec(num_feature/2,num_feature-1) += gen_feature_valued(aug_dim_holder[1],node_attr);
  }
  
  Rcout << "Observed network statistics calculated\n";
  
  for(int l_iter = 0; l_iter < learning_iter; ++l_iter){
    
    mu_sum.zeros(num_feature);
    cov_sum.zeros(num_feature,num_feature);
    for(int iter = 1; iter < y_list.size(); ++iter){
      eta_atleast = eta.subvec(0, num_feature/2-1);
      eta_atmost = eta.subvec(num_feature/2,num_feature-1);
      MCMC_list = gen_MCMC_valued_CD(MCMC_size, n_CD, y_list[iter-1], y_list[iter], pi0, m, eta_atleast, eta_atmost, node_attr);
      MCMC_list_features = gen_features_MCMC_valued(MCMC_list, y_list, iter, node_attr, num_feature);
      mu_sum += arma::mean(MCMC_list_features,0).t();
      cov_sum += cov(MCMC_list_features);
    }
    
    S = g_obs_sum - mu_sum;
    H = -cov_sum;
    eta = eta - inv(H) * S;
    
    Rcout << "learning iter = " << l_iter+1 << "\n";
    Rcout << "gradient =" << S.t();
    Rcout << "eta (newton raphson):\n" << eta << "\n";
    
  }
  
  return eta;
  
}

//[[Rcpp::export]]
arma::vec SE_temporal(int MCMC_size, int n_CD, List y_list, arma::vec eta, arma::vec node_attr, double pi0, int m){
  
  List MCMC_list;
  arma::vec eta_atleast, eta_atmost;
  arma::mat cov_sum, H;
  arma::mat MCMC_list_features;
  int num_feature = eta.size();
  
  arma::vec output;
  eta_atleast = eta.subvec(0, num_feature/2-1);
  eta_atmost = eta.subvec(num_feature/2,num_feature-1);
  
  cov_sum.zeros(num_feature,num_feature);
  for(int iter = 1; iter < y_list.size(); ++iter){
    MCMC_list = gen_MCMC_valued_CD(MCMC_size, n_CD, y_list[iter-1], y_list[iter], pi0, m, eta_atleast, eta_atmost, node_attr);
    MCMC_list_features = gen_features_MCMC_valued(MCMC_list, y_list, iter, node_attr, num_feature);
    cov_sum += cov(MCMC_list_features);
  }
  
  H = -cov_sum;
  output = sqrt(diagvec(-inv(H)));
  
  return output;
  
}


