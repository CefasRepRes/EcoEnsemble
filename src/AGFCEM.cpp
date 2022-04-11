#include <Rcpp.h>

using namespace Rcpp;

//' @title Dummy
//' @description This is a dummy function.
//' @param input a numerical input
//' @details This is a dummy function.
//' @return Nothing.
//' @examples
//' AGFCEM(3.0)
//' @export
// [[Rcpp::export]]
void AGFCEM(double input){
  if(input > 1){
    Rprintf("We are the Boro!!!");
  }else{
    Rprintf("What happened to the frog that broke down? He got TOAD away");
  }
  return;
}
