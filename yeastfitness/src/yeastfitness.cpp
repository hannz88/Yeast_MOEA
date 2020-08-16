#include <algorithm>
#include <parallel/algorithm>
#include <float.h>
#include <math.h>
#include <Rcpp.h>
#include <gsl/gsl_statistics_double.h>
#include <iostream>

// [[Rcpp::plugins(cpp17)]]
using namespace Rcpp;


// [[Rcpp::export]]
double sd_gsl(const NumericMatrix & mx) {
    NumericVector vec = as<NumericVector>(mx);
    double tmp_data[vec.size()+1];
    double result;
    std::copy(vec.begin(), vec.end(), tmp_data);
    result = gsl_stats_sd(tmp_data, 1, vec.size());
    return result;
}

// [[Rcpp::export]]
double one_if_zero(const double & x) {
    if (x == 0.0) {
        return 1.0;
    } else {
        return x;
    }
}

// [[Rcpp::export]]
double sd_cpp(const NumericVector & x){
    double result;
    result = sd(x) * sqrt(((float) x.size() - 1)/ (float) x.size());
    return result;
}

// [[Rcpp::export]]
NumericMatrix standardize_matrix_cpp(const NumericMatrix & matx) {
    NumericVector gene;
    double tmp_mean;
    double tmp_sd;
    double tmp_val;
    int n = matx.nrow();
    int k = matx.ncol();
    NumericVector row_mean = NumericVector(n);
    NumericVector row_sd = NumericVector(n);
    NumericMatrix result = NumericMatrix(n, k);
    for (int i = 0; i < n; i++) { // Not trivially parallelizable
        gene = matx(i, _);
        tmp_mean = mean(gene);
        tmp_sd = sd_cpp(gene);
        row_mean[i] = tmp_mean;
        row_sd[i] = tmp_sd;
    }
    std::for_each(row_sd.begin(), row_sd.end(), one_if_zero);
    for (int i = 0; i < n; i++) { // Not trivially parallelizable
        for (int j = 0; j < k ; j++) {
            tmp_val = matx(i,j) - row_mean[i];
            tmp_val = tmp_val/row_sd[i];
            result(i,j) = tmp_val;
        }
    }
    return result;
}


// [[Rcpp::export]]
double spearman_rho_cpp(const NumericMatrix & mx, int i, int j) {
    NumericMatrix::ConstRow mx_i = mx(i, _);
    NumericMatrix::ConstRow mx_j = mx(j, _);
    double tmp_row1[mx_i.size()+1];
    double tmp_row2[mx_j.size()+1];
    double workspace[2*(mx_i.size()+1)];
    double result;
    std::copy(mx_i.begin(), mx_i.end(), tmp_row1);
    std::copy(mx_j.begin(), mx_j.end(), tmp_row2);
    result = gsl_stats_spearman(tmp_row1, 1, tmp_row2, 1, mx_i.size(), workspace);
    return result;
}

// [[Rcpp::export]]
double pearson_correlation_cpp(const NumericMatrix & mx, int i, int j) {
    NumericMatrix::ConstRow mx_i = mx(i, _);
    NumericMatrix::ConstRow mx_j = mx(j, _);
    double tmp_row1[mx_i.size()+1];
    double tmp_row2[mx_j.size()+1];
    double result;
    std::copy(mx_i.begin(), mx_i.end(), tmp_row1);
    std::copy(mx_j.begin(), mx_j.end(), tmp_row2);
    result = gsl_stats_correlation(tmp_row1, 1, tmp_row2, 1, mx_i.size());
    return result;
}

// [[Rcpp::export]]
double get_rho_cpp(const NumericMatrix & mx) {
    NumericMatrix _mx = NumericMatrix(mx.nrow(), mx.ncol());
    std::copy(mx.begin(), mx.end(), _mx.begin());
    double rho_row_sum = 0;
    int end = _mx.nrow();
    double rho = 0;
    double denominator = 0;
    for (int i = 0; i < end-1; i++){ // Cannot parallelize
        for (int j = i+1; j < end; j++){
            rho = spearman_rho_cpp(_mx, i, j);
            if (isnan(rho))
            {
                Rcout << "get_rho_cpp - Found a NaN " << std::endl;
                Rcout << "i: " << i  << " j:" << j << std::endl;
                _mx(i, 0) += 0.000000000000001; // Unbiased correction for either row being a flat line
                _mx(j, 0) += 0.000000000000001; // Because (Xi-X_mean) == 0 and this causes a Zero Division Error
                rho = spearman_rho_cpp(_mx, i, j);
                if (isnan(rho))
                {
                    Rcout << "get_rho_cpp - Could not correct a NaN " << std::endl;
                }
            }
            rho_row_sum += rho;
        }
    }
    //Rcout << end << std::endl;
    //Rcout << rho_row_sum << std::endl;
    denominator = _mx.nrow() * (_mx.nrow() - 1);
    //Rcout << denominator << std::endl;
    return rho_row_sum/denominator;
}

// [[Rcpp::export]]
double get_sub_SCS_cpp(const NumericMatrix& mx) {
    double min_scs = DBL_MAX;
    int end = mx.nrow();
    double tmp_sum, tmp_scs, row_pearson = 0;
    double abs_sum;
    for (int i = 0;  i < end; i++) { // Cannot parallelize
        tmp_sum = 0;
        for (int j = 0; j < end; j++) {
            if (i==j){
                continue;
            }
            row_pearson = pearson_correlation_cpp(mx, i, j);
            abs_sum = std::abs(row_pearson);
            //Rcout << "abs_sum: " << abs_sum << std::endl;
            tmp_sum += abs_sum;
            //Rcout << "tmp_sum: " << tmp_sum << std::endl;
        }
        tmp_scs = 1 - (tmp_sum/((float)end-1));
        //Rcout << "tmp_scs: " << tmp_scs << std::endl;
        if (tmp_scs < min_scs){
            min_scs = tmp_scs;
        }
        //Rcout << "min_scs: " << min_scs << std::endl;
    }
    return min_scs;
}



// [[Rcpp::export]]
double get_SCS_cpp(const NumericMatrix& mx){
    int n = mx.nrow();
    int k = mx.ncol();
    if ((n <= 1) || (k <= 1) ){
        return 10000;
    }
    double row_scs = get_sub_SCS_cpp(mx);
    NumericMatrix transposed_mx = transpose(mx);
    double col_scs = get_sub_SCS_cpp(transposed_mx);
    double result = std::min(row_scs, col_scs);
    if (isnan(result))
    {
        Rcout << "get_SCS_cpp - Found a NaN " << std::endl;
    }
    return result;
}

// [[Rcpp::export]]
double get_ASR_cpp(const NumericMatrix & mx) {
    int n = mx.nrow();
    int k = mx.ncol();
    if ((n <= 1) || (k <= 1) ){
        return 0;
    }
    //Rcout << "computing row_rho" << std::endl;
    double row_rho = get_rho_cpp(mx);
    //Rcout << "transposing mx" << std::endl;
    NumericMatrix transposed_mx = transpose(mx);
    //Rcout << "computing col_rho" << std::endl;
    double col_rho = get_rho_cpp(transposed_mx);
    //Rcout << "computing max_rho" << std::endl;
    double max_rho = std::max(row_rho, col_rho);
    //Rcout << "computing results" << std::endl;
    double results = std::abs(2*max_rho);
    if (isnan(results))
    {
        Rcout << "get_ASR_cpp - Found a NaN " << std::endl;
    }
    return results;
}

struct minus_mean_divide_by_sd {
    double _mean = 0;
    double _sdev = 0;
    minus_mean_divide_by_sd(double & mean, double & sdev){
        _mean = mean;
        _sdev = sdev;
    }
    void operator()(double & x) const{
        x = (x-_mean)/_sdev;
    }
};

// [[Rcpp::export]]
double get_VET_cpp(const NumericMatrix & matx) {
    NumericMatrix mx = transpose(matx);
    int n = mx.nrow();
    int k = mx.ncol();
    if ((n <= 1) || (k <= 1) ){
        return 10000;
    }
    NumericVector rho = NumericVector(k);
    NumericVector gene;
    NumericMatrix bic_hat = NumericMatrix(n,k);
    double ve, tmp, denominator = 0;
    for (int j = 0; j < k ; j++) { // Not trivially parallelizable
        gene = mx(_, j);
        rho[j] = mean(gene);
    }
    double rho_sd = sd_cpp(rho);
    double rho_mean = mean(rho);
//    Rcout << "rho: " << rho << std::endl;
    rho_sd = one_if_zero(rho_sd);
//    Rcout << "rho after one_if_zero: " << rho << std::endl;
    std::for_each(rho.begin(), rho.end(), minus_mean_divide_by_sd(rho_mean, rho_sd));
//    Rcout << "rho after (x-mean)/sd: " << rho << std::endl;
    bic_hat = standardize_matrix_cpp(mx);
    for (int i = 0; i < n ; i++) {  // Cannot parallelize
        for (int j = 0; j < k; j++) {
            tmp = abs(bic_hat(i,j) - rho[j]);
            ve += tmp;
        }
    }
    denominator = n*k;
    ve = ve/denominator;
    return ve;
}

// [[Rcpp::export]]
double get_vol_seridi(const NumericMatrix & matrix, const NumericMatrix & data_mt, double alpha, double beta){
    int n = matrix.nrow();
    int k = matrix.ncol();
    int x = data_mt.nrow();
    int y = data_mt.ncol();
    double row_coef = alpha * ((double)n/(double)x);
    double col_coef = beta * ((double)k/(double)y);
    double result = row_coef + col_coef;
    return result;
}

// [[Rcpp::export]]
double get_vol_pontes(const NumericMatrix & matrix, double genes_weight, double cond_weight){
    double row_coef;
    double col_coef;
    int n = matrix.nrow();
    int k = matrix.ncol();
    row_coef = -log(n)/(log(n)+genes_weight);
    col_coef = -log(k)/(log(k)+cond_weight);
    return row_coef + col_coef;
}

// [[Rcpp::export]]
String print_test(){
    String s("Hello, World!");
    return s;
}

// [[Rcpp::export]]
double test_mean(const NumericMatrix & mx){
    double result = mean(mx);
    return result;
}

// [[Rcpp::export]]
double test_sd(const NumericMatrix & mx) {
    double sdev = sd(mx);
    return sdev;
}
