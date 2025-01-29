#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export("massCpp")]]
NumericMatrix massCpp(NumericMatrix matriz) {
  int m = matriz.nrow();
  int n = matriz.ncol();
  
  // Crear una matriz de resultados
  NumericMatrix result(m - 1, n - 1);
  
  // Convertir la matriz a long double para cálculos de mayor precisión
  for (int i = 0; i < m - 1; ++i) {
    for (int j = 0; j < n - 1; ++j) {
      long double a = static_cast<long double>(matriz(i, j));
      long double b = static_cast<long double>(matriz(i + 1, j + 1));
      long double c = static_cast<long double>(matriz(i + 1, j));
      long double d = static_cast<long double>(matriz(i, j + 1));
      
      result(i, j) = a + b - c - d;
    }
  }
  
  return result;
}
