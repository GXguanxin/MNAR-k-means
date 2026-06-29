// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


// [[Rcpp::export]]
List kmpp_init_centers_cpp(const arma::mat xmat, const int k){
  int n = xmat.n_rows, p = xmat.n_cols;
  double tmp_dist = 0;
  arma::mat init_cmat = arma::zeros(k,p);
  arma::vec all_index = arma::regspace(0,n-1), tmp_dist_vec = arma::ones(n), init_index(1);
  //arma::vec init_centers_index_vec(k);
  for(int kk=0; kk<k; ++kk){
    //deal with NA in tmp_dist_vec for sample probability
    //arma::vec tmp_dist_vec_copy = tmp_dist_vec;
    //tmp_dist_vec_copy.elem(arma::find_nonfinite(tmp_dist_vec_copy)).zeros();
    //if(arma::accu(tmp_dist_vec_copy) == 0){tmp_dist_vec_copy=arma::ones(n);}
    
    //sample kk-th center with probability=tmp_dist_vec_copy
    init_index = RcppArmadillo::sample(all_index, 1, false, tmp_dist_vec);
    init_cmat(kk,arma::span::all) = xmat(init_index(0),arma::span::all);
    
    for(int i = 0; i < n; ++i){//update weights for sampling
      tmp_dist = sum( pow(xmat(i,arma::span::all) - init_cmat(kk,arma::span::all), 2) );
      if(kk==0){
        tmp_dist_vec(i) = tmp_dist;
      }else if(tmp_dist < tmp_dist_vec(i)){
        tmp_dist_vec(i) = tmp_dist;
      }
    }
  }
  
  List res;
  res["centers"] = init_cmat;
  
  return res;
}



// [[Rcpp::export]]
List pkpod_once_1218cpp(const arma::mat x, const int k, const double lambda,
                const arma::mat init_centers, const int itermax){
  int n = x.n_rows, p = x.n_cols;
  arma::mat missing = arma::conv_to<arma::mat>::from(x != x); // (NAN != NAN) = TRUE，(NAN == NAN) = FALSE
  arma::mat observed = 1.0 - missing;
  arma::mat coefmat = observed + lambda * missing;
  
  double obj_old = 1e-8, obj_new = 0, term1 = 0, term2 = 0; 
  int iter_end = 0;
  arma::uvec cluster(n);
  arma::rowvec masked_center(p), numerator(p), denominator(p);
  arma::vec obj_list(itermax); 
  arma::mat cmat= init_centers, dist_mat = arma::zeros(n,k), Q_umat = arma::zeros(n,k), ab_mat = arma::zeros(n,p);
  arma::umat label_list(n, itermax, arma::fill::zeros);
  arma::cube cmat_list = arma::zeros(k,p,itermax);
  
  arma::mat Pomega_x = x;
  Pomega_x.replace(arma::datum::nan, 0);
  
  for (int iter = 0; iter < itermax; ++iter){
    
    //given tmp_cmat, estimate tmp_cluster
    dist_mat = arma::zeros(n,k);
    for (int i = 0; i < n; ++i){
      for (int kk = 0; kk < k; ++kk){
        //partial distance between x_i and mu_kk on observed dimensions
        masked_center = cmat(kk,arma::span::all)%observed(i,arma::span::all);
        dist_mat(i,kk) = sum( pow( Pomega_x(i,arma::span::all) - masked_center , 2) ); 
      }
    }
    //penalized distance between x_i and mu_kk
    Q_umat = dist_mat + lambda*( missing*( pow(cmat,2).t() ) );
    cluster = arma::index_min(Q_umat,1);
    
    //update cmat
    cmat = arma::zeros(k,p);
    for (int kk = 0; kk < k; ++kk){
      numerator = sum( Pomega_x.rows( arma::find(cluster == kk) ) , 0);
      denominator = sum( coefmat.rows( arma::find(cluster == kk) ) , 0);
      cmat.row(kk) = numerator/denominator;
    }
    
    //record
    ab_mat = cmat.rows(cluster);
    term1 = arma::accu( pow( (Pomega_x - ab_mat)%observed ,2) );
    term2 = lambda*( arma::accu( pow( ab_mat%missing ,2) ) );
    obj_new = term1 + term2;
    
    //stopping
    if(iter > 0){
      if( abs(obj_new-obj_old)/abs(obj_old) < 1e-5 ){
        obj_list(iter)=obj_new;
        label_list.col(iter)=cluster;
        cmat_list.slice(iter)=cmat;
        iter_end=iter;
        break;}
      }
    //update
    obj_old=obj_new;
    iter_end=iter;
    obj_list(iter)=obj_new;
    label_list.col(iter)=cluster;
    cmat_list.slice(iter)=cmat;
  }
  
  
  List result;
  result["cluster"] = cluster+1;
  result["centers"] = cmat;
  result["obj_val"] = obj_new;
  result["term1"] = term1;
  result["term2"] = term2;
  result["obj_list"] = obj_list.subvec(0,iter_end);
  result["label_list"] = label_list.cols(0,iter_end)+1;
  result["cmat_list"] = cmat_list.slices(0,iter_end);
  return(result);
}
  
  
  