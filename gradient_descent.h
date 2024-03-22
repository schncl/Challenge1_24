#ifndef GRADIENT_DESCENT_H
#define GRADIENT_DESCENT_H

#include <vector>
#include <functional>
#include <string>

//some typedef for convenience
using functions=std::function<double(const std::vector<double>&)>;
using vec_function=std::function<std::vector<double>(const std::vector<double>&)>;
using Point=std::vector<double>;

// Enum to specify step size strategy
enum class StepSizeStrategy {
    ExponentialDecay,
    InverseDecay,
    Armijo,
    Constant
};

//Enum to speify optimization strat
enum class OptimizationStrategy {
    gradient_descent,
    gradient_descent_with_momentum,
    nesterov
};


// Struct to hold optimization parameters
struct OptimizationParameters {
    std::vector<double> initial_point;
    double epsilon_r;
    double epsilon_s;
    int max_iterations;
    double initial_step;
    functions function;
    vec_function gradient;
    StepSizeStrategy step_size_strategy;
    OptimizationStrategy method;
    double mu;
    double sigma;
};
//function that calculate the l2 norm of the difference of v1 and v2
double vectorNormdiff(const std::vector<double>& v1,const std::vector<double>& v2);
//function that computes the l2 norm of x
double vectorNorm(const std::vector<double>& x);
// Various method for solving the optimization problem. It takes an imput a struct
//that holds various parameters, and where you can choose the update strategy for the learning rate
std::vector<Point> gradient_descent(const OptimizationParameters& params);
std::vector<Point> gradient_descent_w_momentum(const OptimizationParameters& params);
std::vector<Point> nesterov(const OptimizationParameters& params);
//Various strategies to update the learning rate;
double ArmijoStepSize(const OptimizationParameters& params,const double a0,const Point& xk);
double ExpDecay(const OptimizationParameters& params,const double a0,const Point& xk);
double Inverse_decay(const OptimizationParameters& params,const double a0,const Point& xk);
double update_ak(const OptimizationParameters& params,const double a0,const double ak,const Point& xk,std::size_t k);
//Function that collect al of this, you can choose method and l_r update strategy.
std::vector<Point> optimize(const OptimizationParameters& params);

#endif
