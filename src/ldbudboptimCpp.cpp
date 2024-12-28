#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ldbudboptimCpp(NumericVector z_values, 
                      Function F1, 
                      Function F2, 
                      Function F3, 
                      Function F4, 
                      NumericVector x_limit, 
                      NumericVector y_limit, 
                      double delta = 0.1,
                      Nullable<NumericVector> sopY = R_NilValue) {
  
  // Encontrar mínimo y máximo de x_limit y y_limit
  double min_x_limit = x_limit[0];
  double max_x_limit = x_limit[1];
  double min_y_limit = y_limit[0];
  double max_y_limit = y_limit[1];

  int num_z = z_values.size();
  
  // Vectores para almacenar los resultados de cotainf y cotasup para cada z
  NumericVector cotainf(num_z);
  NumericVector cotasup(num_z);

  // Iterar sobre los valores de z
  for (int j = 0; j < num_z; ++j) {
    double z = z_values[j];

    // Generar valores de u basados en delta, x_limit, y_limit y z
    double min_val = std::min(min_x_limit, min_y_limit + z) - 3 * delta;
    double max_val = std::max(max_x_limit, max_y_limit + z) + 3 * delta;
    int n_u = static_cast<int>(std::ceil((max_val - min_val) / delta)) + 1;

    // Crear un vector de valores de u
    NumericVector u(n_u);
    for (int i = 0; i < n_u; ++i) {
      u[i] = min_val + i * delta;
    }
  
    // Si sopY no es NULL, agregar sus valores a u
    if (sopY.isNotNull()) {
      NumericVector sopY_vec(sopY);
      u = union_(u, sopY_vec); // Usamos la función union_ de Rcpp para garantizar unicidad y eficiencia
    }

    // Calcular shifted_u
    NumericVector shifted_u = u - z;

    // Evaluar F1 y F3 sobre el vector u
    NumericVector F1_u = F1(u);
    NumericVector F3_u = F3(u);

    // Evaluar F2 y F4 sobre shifted_u
    NumericVector F2_shifted_u = F2(shifted_u);
    NumericVector F4_shifted_u = F4(shifted_u);

    // Calcular val12 y val34 para todos los valores de u
    NumericVector val12 = F1_u - F2_shifted_u;
    NumericVector val34 = F3_u - F4_shifted_u;

    // Añadir 0 a val12 y val34
    val12.push_back(0.0);
    val34.push_back(0.0);

    // Actualizar cotainf y cotasup para el valor actual de z
    double current_cotainf = max(val12);
    double current_cotasup = 1 + min(val34);

    // Almacenar los resultados para el valor actual de z
    cotainf[j] = current_cotainf;
    cotasup[j] = current_cotasup;
  }
  
  // Devolver resultados
  return List::create(Named("z") = z_values,
                      Named("ldb") = cotainf, 
                      Named("udb") = cotasup);
}
