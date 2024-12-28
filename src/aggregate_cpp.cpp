#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix aggregate_cpp(NumericVector values, NumericVector groups, NumericVector levelgroup) {
  // Map para almacenar las sumas de valores para cada nivel en levelgroup
  std::unordered_map<double, double> group_sums;

  // Iterar sobre los vectores y acumular las sumas
  for (int i = 0; i < values.size(); ++i) {
    double group = groups[i];
    double value = values[i];
    group_sums[group] += value;
  }
  
  // Crear una lista para almacenar los resultados
  std::vector<double> levelgroup_sorted(levelgroup.begin(), levelgroup.end());
  std::sort(levelgroup_sorted.begin(), levelgroup_sorted.end()); // Ordenar los niveles de levelgroup

  // Crear la matriz de resultados
  NumericMatrix result(levelgroup_sorted.size(), 2);

  // Inicializar el índice para la matriz de resultados
  int idx = 0;
  for (double level : levelgroup_sorted) {
    result(idx, 0) = level;
    result(idx, 1) = group_sums[level];
    ++idx;
  }
  
  return result;
}
