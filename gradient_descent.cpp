#include "gradient_descent.h"
#include <iostream>
#include <cmath>

using function=std::function<double(const std::vector<double>&)>;
using vec_function =std::function<std::vector<double>(const std::vector<double>&)>;
using Point=std::vector<double>;

double vectorNormdiff(const std::vector<double>& v1,const std::vector<double>& v2) {
    double norm = 0.0;
    for (auto i=0;i<v1.size();++i){
        norm += (v1[i]-v2[i])*(v1[i]-v2[i]);
    }
    return std::sqrt(norm);
}
double vectorNorm(const std::vector<double>& x){
    double norm=0.;
    for(auto i=0;i<x.size();++i){
        norm+=x[i]*x[i];
    }
    return std::sqrt(norm);
}
bool condition(function f,const vec_function grad_f, std::vector<double> xk,const double a0,double sigma){
    std::vector<double> new_point;
    new_point.resize(xk.size());
    auto grad=grad_f(xk);
    for (auto i = 0; i < xk.size(); ++i)
        new_point[i] = xk[i] - a0 * grad[i];

    double lhs=-f(new_point)+f(xk);
    double rhs=sigma*a0*vectorNorm(grad_f(xk))*vectorNorm(grad_f(xk));
    return lhs>=rhs;
}
double ArmijoStepSize(const OptimizationParameters&params,const double a0,const Point& xk) {
    double ak=a0;
    if(params.sigma<0 || params.sigma>0.5){
        std::cerr<<"Sigma must be between 0 and 0.5"<<std::endl;
    }
    bool flag=true;
    while(flag){
        if(condition(params.function,params.gradient,xk,ak,params.sigma)){
            flag=false;
        }
        ak=ak/2;
    }
    return ak;
}
double Inverse_decay(const OptimizationParameters& params, const double a0, const Point& xk,std::size_t k) {
    double alpha = a0;
    alpha /= (1 + params.mu * k);  
    return alpha;
}
double ExpDecay(const OptimizationParameters& params, const double a0, const Point& xk,std::size_t k) {
    double alpha = a0;
    alpha *= std::exp(-params.mu * k); 
    return alpha;
}
double update_ak(const OptimizationParameters& params,const double a0,const double ak,const Point& xk,std::size_t k){
    switch(params.step_size_strategy) {
        case StepSizeStrategy::Armijo:
            return ArmijoStepSize(params, ak, xk);
        case StepSizeStrategy::ExponentialDecay:
            return ExpDecay(params,a0,xk,k);
        case StepSizeStrategy::InverseDecay:
            return Inverse_decay(params,a0,xk,k);
        case StepSizeStrategy::Constant:
            return ak;    
        default:
            return ak;

    }
}
/////////////////////GRADIENT DESCENT//////////////////////////
std::vector<Point> gradient_descent(const OptimizationParameters& params){
    std::vector<Point> ret;
    Point x_new{params.initial_point};
    Point x_old{params.initial_point};
    ret.emplace_back(x_old);
    bool flag=true;
    auto k=0;
    auto a0=params.initial_step;
    double ak=a0;
    for(k=0;k<params.max_iterations && flag;++k){
        auto grad=params.gradient(x_old);
        for (auto i=0;i<params.initial_point.size();++i){
            x_new[i]=x_old[i]-ak*grad[i];
        }
        ret.emplace_back(x_new);
        auto f_old=params.function(x_old);
        auto f_new=params.function(x_new);
        if(vectorNormdiff(x_new,x_old)<params.epsilon_s || std::abs(f_new-f_old)<params.epsilon_r){
            flag=false;
        }
        ak= update_ak(params,ak,a0,x_new,k);
        x_old=x_new;
    }
    if(k==params.max_iterations){
        std::cout<<"The algorithm did not converge..."<<std::endl;
    }
    std::cout<<"Reached in:"<<k<<" iterations"<<std::endl;
    std::cout<<"The minimum is. "<<ret.back()[0]<<" "<<ret.back()[1]<<"."<<std::endl;
    return ret;
}
/////////////////////////////////////////////////////////////////////////////////
///////////////////////////GRADIENT DESCENT WITH MOMENTUM//////////////////////////////////////////////////////
std::vector<Point> gradient_descent_w_momentum(const OptimizationParameters& params){
    std::vector<Point> ret;

    Point x_new{params.initial_point};
    Point x_old{params.initial_point};

    ret.emplace_back(x_old);

    bool flag=true;
    auto k=0;
    auto a0=params.initial_step;
    double ak=a0;
    auto grad=params.gradient(x_old);
    std::vector<double> d_old;
    d_old.resize(x_old.size());
    for(auto i=0;i<d_old.size();++i){
        d_old[i]=-1*ak*grad[i];
    }
    std::vector<double> d_new{d_old};

    for(k=0;k<params.max_iterations && flag;++k) {

        for (auto i=0;i<params.initial_point.size();++i){
            x_new[i]=x_old[i]+d_old[i];
        }
        grad=params.gradient(x_new);
        for(auto i=0;i<d_old.size();++i){
            d_new[i]=params.mu*d_old[i]-ak*grad[i];
        }
        ret.emplace_back(x_new);
        auto f_old=params.function(x_old);
        auto f_new=params.function(x_new);
        if(vectorNormdiff(x_new,x_old)<params.epsilon_s || std::abs(f_new-f_old)<params.epsilon_r){
            flag=false;
        }
        x_old=x_new;
        d_old=d_new;
        ak=update_ak(params,ak,a0,x_new,k);
    }
    if(k==params.max_iterations){
        std::cout<<"The algorithm did not converge..."<<std::endl;
    }
    else {
        std::cout << "Reached in:" << k << " iterations" << std::endl;
        std::cout<<"The minimum is. "<<ret.back()[0]<<" "<<ret.back()[1]<<"."<<std::endl;
    }


    return ret;
}
/////////////////////////////////////////////////////////////////////////////////
///////////////////////////NESTEROV//////////////////////////////////////////////////////
std::vector<Point> nesterov(const OptimizationParameters& params){
    std::vector<Point> ret;

    Point x_old{params.initial_point};
    Point x_new;
    x_new.resize(x_old.size());

    ret.emplace_back(x_old);

    bool flag=true;
    auto k=0;
    double a0=params.initial_step;
    auto ak=a0;
    auto grad=params.gradient(x_old);
    for(auto i=0;i<x_new.size();++i){
        x_new[i]=x_old[i]-ak*grad[i];
    }
    ret.emplace_back(x_new);
    std::vector<double> y;
    y.resize(x_old.size());
    for(k=1;k<params.max_iterations && flag;++k) {
        for(auto i=0;i<y.size();++i){
            y[i]=x_new[i]+params.mu*(x_new[i]-ret[k-1][i]);
        }
        grad=params.gradient(y);
        for (auto i=0;i<params.initial_point.size();++i){
            x_new[i]=y[i]-ak*grad[i];
        }
        ret.emplace_back(x_new);
        auto f_old=params.function(x_old);
        auto f_new=params.function(x_new);
        if(vectorNormdiff(x_new,x_old)<params.epsilon_s || std::abs(f_new-f_old)<params.epsilon_r){
            flag=false;
        }
        ak=update_ak(params,ak,a0,x_new,k);
        x_old=x_new;
    }
    if(k==params.max_iterations){
        std::cout<<"The algorithm did not converge..."<<std::endl;
    }
    else {
        std::cout << "Reached in:" << k << " iterations" << std::endl;
        std::cout<<"The minimum is. "<<ret.back()[0]<<" "<<ret.back()[1]<<"."<<std::endl;
    }

    return ret;
}
///////////////////////////////////////////////////////////////////////////////
std::vector<Point> optimize(const OptimizationParameters& params){
    switch(params.method) {
        case OptimizationStrategy::gradient_descent:
            return gradient_descent(params);
        case OptimizationStrategy::gradient_descent_with_momentum:
            return gradient_descent_w_momentum(params);
        case OptimizationStrategy::nesterov:
            return nesterov(params);
        default:
            std::cerr<<"Use a valid method"<<std::endl;
            break;

    }
}

