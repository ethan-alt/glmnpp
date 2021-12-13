// Generated by rstantools.  Do not edit by hand.

/*
    glmnpp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glmnpp is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glmnpp.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_glm_pp_fixed_poisson_post_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_glm_pp_fixed_poisson_post");
    reader.add_event(2, 2, "include", "/functions/link.stan");
    reader.add_event(2, 0, "start", "/functions/link.stan");
    reader.add_event(25, 23, "end", "/functions/link.stan");
    reader.add_event(25, 3, "restart", "model_glm_pp_fixed_poisson_post");
    reader.add_event(88, 64, "end", "model_glm_pp_fixed_poisson_post");
    return reader;
}
template <typename T0__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
linkinv(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& eta,
            const int& mu_link, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        current_statement_begin__ = 11;
        if (as_bool(logical_eq(mu_link, 1))) {
            current_statement_begin__ = 11;
            return stan::math::promote_scalar<fun_return_scalar_t__>(eta);
        } else if (as_bool(logical_eq(mu_link, 2))) {
            current_statement_begin__ = 12;
            return stan::math::promote_scalar<fun_return_scalar_t__>(stan::math::exp(eta));
        } else if (as_bool(logical_eq(mu_link, 3))) {
            current_statement_begin__ = 13;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv_logit(eta));
        } else if (as_bool(logical_eq(mu_link, 4))) {
            current_statement_begin__ = 14;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv(eta));
        } else if (as_bool(logical_eq(mu_link, 5))) {
            current_statement_begin__ = 15;
            return stan::math::promote_scalar<fun_return_scalar_t__>(Phi(eta));
        } else if (as_bool(logical_eq(mu_link, 6))) {
            current_statement_begin__ = 16;
            return stan::math::promote_scalar<fun_return_scalar_t__>(add(divide(stan::math::atan(eta), stan::math::pi()), 0.5));
        } else if (as_bool(logical_eq(mu_link, 7))) {
            current_statement_begin__ = 17;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv_cloglog(eta));
        } else if (as_bool(logical_eq(mu_link, 8))) {
            current_statement_begin__ = 18;
            return stan::math::promote_scalar<fun_return_scalar_t__>(square(eta));
        } else if (as_bool(logical_eq(mu_link, 9))) {
            current_statement_begin__ = 19;
            return stan::math::promote_scalar<fun_return_scalar_t__>(inv(stan::math::sqrt(eta)));
        }
        current_statement_begin__ = 20;
        return stan::math::promote_scalar<fun_return_scalar_t__>(eta);
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct linkinv_functor__ {
    template <typename T0__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& eta,
            const int& mu_link, std::ostream* pstream__) const {
        return linkinv(eta, mu_link, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_glm_pp_fixed_poisson_post
  : public stan::model::model_base_crtp<model_glm_pp_fixed_poisson_post> {
private:
        int N;
        int N0;
        int p;
        std::vector<int> y;
        std::vector<int> y0;
        matrix_d X;
        matrix_d X0;
        vector_d beta0;
        matrix_d Sigma0;
        int link;
        int incl_offset;
        vector_d offset;
        vector_d offset0;
        double a0;
public:
    model_glm_pp_fixed_poisson_post(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_glm_pp_fixed_poisson_post(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_glm_pp_fixed_poisson_post_namespace::model_glm_pp_fixed_poisson_post";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 30;
            context__.validate_dims("data initialization", "N", "int", context__.to_vec());
            N = int(0);
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            N = vals_i__[pos__++];
            check_greater_or_equal(function__, "N", N, 0);
            current_statement_begin__ = 31;
            context__.validate_dims("data initialization", "N0", "int", context__.to_vec());
            N0 = int(0);
            vals_i__ = context__.vals_i("N0");
            pos__ = 0;
            N0 = vals_i__[pos__++];
            check_greater_or_equal(function__, "N0", N0, 0);
            current_statement_begin__ = 32;
            context__.validate_dims("data initialization", "p", "int", context__.to_vec());
            p = int(0);
            vals_i__ = context__.vals_i("p");
            pos__ = 0;
            p = vals_i__[pos__++];
            check_greater_or_equal(function__, "p", p, 0);
            current_statement_begin__ = 33;
            validate_non_negative_index("y", "N", N);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(N));
            y = std::vector<int>(N, int(0));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_k_0_max__ = N;
            for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                y[k_0__] = vals_i__[pos__++];
            }
            size_t y_i_0_max__ = N;
            for (size_t i_0__ = 0; i_0__ < y_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "y[i_0__]", y[i_0__], 0);
            }
            current_statement_begin__ = 34;
            validate_non_negative_index("y0", "N0", N0);
            context__.validate_dims("data initialization", "y0", "int", context__.to_vec(N0));
            y0 = std::vector<int>(N0, int(0));
            vals_i__ = context__.vals_i("y0");
            pos__ = 0;
            size_t y0_k_0_max__ = N0;
            for (size_t k_0__ = 0; k_0__ < y0_k_0_max__; ++k_0__) {
                y0[k_0__] = vals_i__[pos__++];
            }
            size_t y0_i_0_max__ = N0;
            for (size_t i_0__ = 0; i_0__ < y0_i_0_max__; ++i_0__) {
                check_greater_or_equal(function__, "y0[i_0__]", y0[i_0__], 0);
            }
            current_statement_begin__ = 35;
            validate_non_negative_index("X", "N", N);
            validate_non_negative_index("X", "p", p);
            context__.validate_dims("data initialization", "X", "matrix_d", context__.to_vec(N,p));
            X = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N, p);
            vals_r__ = context__.vals_r("X");
            pos__ = 0;
            size_t X_j_2_max__ = p;
            size_t X_j_1_max__ = N;
            for (size_t j_2__ = 0; j_2__ < X_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X_j_1_max__; ++j_1__) {
                    X(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 36;
            validate_non_negative_index("X0", "N0", N0);
            validate_non_negative_index("X0", "p", p);
            context__.validate_dims("data initialization", "X0", "matrix_d", context__.to_vec(N0,p));
            X0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(N0, p);
            vals_r__ = context__.vals_r("X0");
            pos__ = 0;
            size_t X0_j_2_max__ = p;
            size_t X0_j_1_max__ = N0;
            for (size_t j_2__ = 0; j_2__ < X0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < X0_j_1_max__; ++j_1__) {
                    X0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 37;
            validate_non_negative_index("beta0", "p", p);
            context__.validate_dims("data initialization", "beta0", "vector_d", context__.to_vec(p));
            beta0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(p);
            vals_r__ = context__.vals_r("beta0");
            pos__ = 0;
            size_t beta0_j_1_max__ = p;
            for (size_t j_1__ = 0; j_1__ < beta0_j_1_max__; ++j_1__) {
                beta0(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 38;
            validate_non_negative_index("Sigma0", "p", p);
            validate_non_negative_index("Sigma0", "p", p);
            context__.validate_dims("data initialization", "Sigma0", "matrix_d", context__.to_vec(p,p));
            Sigma0 = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(p, p);
            vals_r__ = context__.vals_r("Sigma0");
            pos__ = 0;
            size_t Sigma0_j_2_max__ = p;
            size_t Sigma0_j_1_max__ = p;
            for (size_t j_2__ = 0; j_2__ < Sigma0_j_2_max__; ++j_2__) {
                for (size_t j_1__ = 0; j_1__ < Sigma0_j_1_max__; ++j_1__) {
                    Sigma0(j_1__, j_2__) = vals_r__[pos__++];
                }
            }
            current_statement_begin__ = 39;
            context__.validate_dims("data initialization", "link", "int", context__.to_vec());
            link = int(0);
            vals_i__ = context__.vals_i("link");
            pos__ = 0;
            link = vals_i__[pos__++];
            check_greater_or_equal(function__, "link", link, 1);
            check_less_or_equal(function__, "link", link, 9);
            current_statement_begin__ = 40;
            context__.validate_dims("data initialization", "incl_offset", "int", context__.to_vec());
            incl_offset = int(0);
            vals_i__ = context__.vals_i("incl_offset");
            pos__ = 0;
            incl_offset = vals_i__[pos__++];
            check_greater_or_equal(function__, "incl_offset", incl_offset, 0);
            check_less_or_equal(function__, "incl_offset", incl_offset, 1);
            current_statement_begin__ = 41;
            validate_non_negative_index("offset", "N", N);
            context__.validate_dims("data initialization", "offset", "vector_d", context__.to_vec(N));
            offset = Eigen::Matrix<double, Eigen::Dynamic, 1>(N);
            vals_r__ = context__.vals_r("offset");
            pos__ = 0;
            size_t offset_j_1_max__ = N;
            for (size_t j_1__ = 0; j_1__ < offset_j_1_max__; ++j_1__) {
                offset(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 42;
            validate_non_negative_index("offset0", "N0", N0);
            context__.validate_dims("data initialization", "offset0", "vector_d", context__.to_vec(N0));
            offset0 = Eigen::Matrix<double, Eigen::Dynamic, 1>(N0);
            vals_r__ = context__.vals_r("offset0");
            pos__ = 0;
            size_t offset0_j_1_max__ = N0;
            for (size_t j_1__ = 0; j_1__ < offset0_j_1_max__; ++j_1__) {
                offset0(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 43;
            context__.validate_dims("data initialization", "a0", "double", context__.to_vec());
            a0 = double(0);
            vals_r__ = context__.vals_r("a0");
            pos__ = 0;
            a0 = vals_r__[pos__++];
            check_greater_or_equal(function__, "a0", a0, 0);
            check_less_or_equal(function__, "a0", a0, 1);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 48;
            validate_non_negative_index("beta", "p", p);
            num_params_r__ += p;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_glm_pp_fixed_poisson_post() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 48;
        if (!(context__.contains_r("beta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable beta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("beta");
        pos__ = 0U;
        validate_non_negative_index("beta", "p", p);
        context__.validate_dims("parameter initialization", "beta", "vector_d", context__.to_vec(p));
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            beta(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(beta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable beta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 48;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> beta;
            (void) beta;  // dummy to suppress unused var warning
            if (jacobian__)
                beta = in__.vector_constrain(p, lp__);
            else
                beta = in__.vector_constrain(p);
            // model body
            {
            current_statement_begin__ = 53;
            validate_non_negative_index("eta", "((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N : 0 )", ((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N : 0 ));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta(((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N : 0 ));
            stan::math::initialize(eta, DUMMY_VAR__);
            stan::math::fill(eta, DUMMY_VAR__);
            current_statement_begin__ = 54;
            validate_non_negative_index("eta0", "((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N0 : 0 )", ((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N0 : 0 ));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> eta0(((primitive_value(logical_neq(link, 2)) || primitive_value(logical_eq(incl_offset, 1))) ? N0 : 0 ));
            stan::math::initialize(eta0, DUMMY_VAR__);
            stan::math::fill(eta0, DUMMY_VAR__);
            current_statement_begin__ = 55;
            validate_non_negative_index("mu", "(logical_neq(link, 2) ? N : 0 )", (logical_neq(link, 2) ? N : 0 ));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu((logical_neq(link, 2) ? N : 0 ));
            stan::math::initialize(mu, DUMMY_VAR__);
            stan::math::fill(mu, DUMMY_VAR__);
            current_statement_begin__ = 56;
            validate_non_negative_index("mu0", "(logical_neq(link, 2) ? N0 : 0 )", (logical_neq(link, 2) ? N0 : 0 ));
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> mu0((logical_neq(link, 2) ? N0 : 0 ));
            stan::math::initialize(mu0, DUMMY_VAR__);
            stan::math::fill(mu0, DUMMY_VAR__);
            current_statement_begin__ = 59;
            lp_accum__.add(multi_normal_log<propto__>(beta, beta0, Sigma0));
            current_statement_begin__ = 64;
            if (as_bool((primitive_value(logical_eq(link, 2)) && primitive_value(logical_eq(incl_offset, 0))))) {
                current_statement_begin__ = 65;
                lp_accum__.add((a0 * poisson_log_glm_log(y0, X0, 0, beta)));
                current_statement_begin__ = 66;
                lp_accum__.add(poisson_log_glm_log(y, X, 0, beta));
            } else {
                current_statement_begin__ = 69;
                stan::math::assign(eta, multiply(X, beta));
                current_statement_begin__ = 70;
                stan::math::assign(eta0, multiply(X0, beta0));
                current_statement_begin__ = 71;
                if (as_bool(logical_eq(incl_offset, 1))) {
                    current_statement_begin__ = 72;
                    stan::math::assign(eta, add(eta, offset));
                    current_statement_begin__ = 73;
                    stan::math::assign(eta0, add(eta0, offset0));
                }
                current_statement_begin__ = 75;
                if (as_bool(logical_eq(link, 2))) {
                    current_statement_begin__ = 76;
                    lp_accum__.add((a0 * poisson_log_log(y0, eta0)));
                    current_statement_begin__ = 77;
                    lp_accum__.add(poisson_log_log(y, eta));
                } else {
                    current_statement_begin__ = 80;
                    stan::math::assign(mu, linkinv(eta, link, pstream__));
                    current_statement_begin__ = 81;
                    stan::math::assign(mu0, linkinv(eta0, link, pstream__));
                    current_statement_begin__ = 82;
                    lp_accum__.add((a0 * poisson_log(y0, mu0)));
                    current_statement_begin__ = 83;
                    lp_accum__.add(poisson_log(y, mu));
                }
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("beta");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dims__.push_back(p);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_glm_pp_fixed_poisson_post_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        Eigen::Matrix<double, Eigen::Dynamic, 1> beta = in__.vector_constrain(p);
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            vars__.push_back(beta(j_1__));
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_glm_pp_fixed_poisson_post";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        size_t beta_j_1_max__ = p;
        for (size_t j_1__ = 0; j_1__ < beta_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "beta" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_glm_pp_fixed_poisson_post_namespace::model_glm_pp_fixed_poisson_post stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
