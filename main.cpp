#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include "gradient_descent.h"
#include <fstream>
#include "GetPot"
/////////////////////////DEFINING FUNCTIONS; I CAN TRY AND DO IT VIA MUPARSERX ..../////////////////////////////////////////////
auto f = [](const std::vector<double>& x) {
    return x[0] * x[1] + 4 * std::pow(x[0], 4) + x[1] * x[1] + 3 * x[0];
};

// Define the gradient of the function f(x) as a lambda function
auto grad_f = [](const std::vector<double>& x) -> std::vector<double> {
    std::vector<double> gradient(x.size());
    gradient[0] = x[1] + 16 * std::pow(x[0], 3) + 3;
    gradient[1] = x[0] + 2 * x[1];
    return gradient;
};
///////////////////////////MAIN FUNCTION, WE EXECUTE OUR TEST///////////////////////////////////////////////
int main(int argc, char **argv) {
//Read parameters from a GetPot file. We still have to define:
//f,grad_f,the optimization strategy and the update strategy on ak.




    GetPot command_line(argc, argv);
    const std::string filename = command_line.follow("data", 2, "-f", "--file");
    GetPot datafile(filename.c_str());
  
    // Initial point
    double x=datafile("x0",0.);
    double y=datafile("y0",0.);
    Point x0{x,y};
    // Tolerances
    double epsilon_r = datafile("eps_r",1e-6);
    double epsilon_s = datafile("eps_s",1e-6);
    // Maximum number of iterations
    int max_iterations = datafile("N",1000);
    // Initial step size
    double alpha0 = datafile("a0",1e-2);
    // Create a struct to hold all parameters
    double mu = datafile("mu",0.01);
    double sigma = datafile("sigma",0.125);


    //TO BE IMPLEMENTED: MUPARSER AND MUPARSERX 


    //Aggregating all parameters into a struct
    OptimizationParameters params = {
            x0,
            epsilon_r,
            epsilon_s,
            max_iterations,
            alpha0,
            f,
            grad_f,
            StepSizeStrategy::Armijo,// Use Armijo rule for step size selection
            OptimizationStrategy::gradient_descent,
            mu,
            sigma
    };


   

    // Perform optimization and printing out minimum point and minimum value at that point.
    std::vector<Point> min_point = optimize(params);
    std::cout<<min_point.back()[0]<<" "<<min_point.back()[1]<<std::endl;
    std::cout<<params.function(min_point.back())<<std::endl;

/*
    //Another Test
    OptimizationParameters params2 = {
            x0,
            epsilon_r,
            epsilon_s,
            max_iterations,
            alpha0,
            f,
            grad_f,
            StepSizeStrategy::ExponentialDecay, // Use Armijo rule for step size selection
            OptimizationStrategy::nesterov,
            mu,
            sigma
    };

    std::vector<Point> min_point2 = optimize(params2);
    std::cout<<min_point2.back()[0]<<" "<<min_point2.back()[1]<<std::endl;
    std::cout<<params.function(min_point2.back())<<std::endl;

	//Yet another test
   	 OptimizationParameters params3 = {
            x0,
            epsilon_r,
            epsilon_s,
            max_iterations,
            alpha0,
            f,
            grad_f,
            StepSizeStrategy::Armijo, // Use Armijo rule for step size selection
            OptimizationStrategy::gradient_descent,
            mu,
            sigma
    };

    std::vector<Point> min_point3 = optimize(params3);
    std::cout<<min_point3.back()[0]<<" "<<min_point3.back()[1]<<std::endl;
    std::cout<<params.function(min_point3.back())<<std::endl;

    std::vector<Point> min_point3 = gradient_descent(params);
    std::cout<<min_point3.back()[0]<<" "<<min_point3.back()[1]<<std::endl;
    std::cout<<params.function(min_point3.back())<<std::endl;

  */



    std::ofstream file("results.dat");
    if (file.is_open()) {
        // Write the header
        file << "iteration\tx_k\ty_k\tf(x_k,y_k)\n";
        // Write the results
        for (int i = 0; i < min_point.size(); ++i) {
            file <<i<< "\t" << min_point[i][0] << "\t" <<  min_point[i][1]<<"\t"<<f(min_point[i])<< "\n";
        }
        file.close();
        std::cout << "Results written to file.\n";
    }
    else {
        std::cerr << "Unable to open file.\n";
    }


    return 0;
}
