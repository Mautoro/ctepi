#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List CTEprobcpp(int n1, int n2, double a1, double b1, double a2, double b2, NumericVector p) {
  
  int num_rows = (n1 + 1) * (n2 + 1);
  NumericMatrix results(num_rows, 5 + p.size());

  // Inicializar las columnas fijas (w1, w2, u1, u2, z)
  int row_index = 0;
  for (int w1 = 0; w1 <= n1; w1++) {
    for (int w2 = 0; w2 <= n2; w2++) {
      double u1 = a1 + b1 * w1;
      double u2 = a2 + b2 * w2;
      double diff = u1 - u2;

      results(row_index, 0) = w1;
      results(row_index, 1) = w2;
      results(row_index, 2) = u1;
      results(row_index, 3) = u2;
      results(row_index, 4) = diff;

      for (int i = 0; i < p.size(); i++) {
        // Calcular dbinom para cada valor de p
        double prob_w1 = R::dbinom(w1, n1, p[i], 0);
        double prob_w2 = R::dbinom(w2, n2, p[i], 0);

        double prob = prob_w1 * prob_w2;

        results(row_index, 5 + i) = prob;
      }

      row_index++;
    }
  }

  CharacterVector namesresults(5 + p.size());
  namesresults[0] = "w1";
  namesresults[1] = "w2";
  namesresults[2] = "u1";
  namesresults[3] = "u2";
  namesresults[4] = "z";
  for (int i = 0; i < p.size(); i++) {
    namesresults[5 + i] = "p" + std::to_string(p[i]);
  }

  return List::create(Named("results") = results,
                      Named("namesresults") = namesresults);
}
