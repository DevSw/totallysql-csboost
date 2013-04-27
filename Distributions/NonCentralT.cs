using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class non_central_T_distribution : distribution
    {
        double m_df, m_lambda;

        public non_central_T_distribution(double i, double lambda)
        {
            m_df = i;
            m_lambda = lambda;
            check_parameters();
        }

        public override void check_parameters()
        {
            if (m_df <= 0 || double.IsInfinity(m_df)) throw new ArgumentException(string.Format("Degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df));
            if (m_lambda < 0 || double.IsInfinity(m_lambda)) throw new ArgumentException(string.Format("Non-centrality argument must be a finite number >= 0 (got {0:G}).", m_lambda));
        }

        public override bool discrete() { return false; }

        public override bool symmetric() { return m_lambda == 0; }

        public double degrees_of_freedom()
        {
            return m_df;
        }

        public double non_centrality()
        {
            return m_lambda;
        }

        public override XMath.pair<double> range()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override XMath.pair<double> support()
        { // Range of permissible values for random variable x.
            return new XMath.pair<double>(-double.MaxValue, double.MaxValue);
        }

        public override double pdf(double t)
        {
            double v = degrees_of_freedom();
            double l = non_centrality();

            return non_central_t_pdf(v, l, t);
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double pdf_inv(double p, bool RHS)
        {
            base.pdf_inv(p, RHS);
            if (p == 0) return RHS ? double.MaxValue : -double.MaxValue;
            double lbound, ubound, result;
            if (RHS)
            {
                lbound = mode();
                ubound = support().v2;
                result = find_pdf_inv(p, lbound, ubound, false);
            }
            else
            {
                lbound = support().v1;
                ubound = mode();
                result = find_pdf_inv(p, lbound, ubound, true);
            }
            return result;
        }

        public override double cdf(double t)
        {
            base.cdf(t);
            double v = degrees_of_freedom();
            double l = non_centrality();

            return non_central_t_cdf(v, l, t, false);
        }

        public override double cdfc(double t)
        {
            base.cdfc(t);
            double v = degrees_of_freedom();
            double l = non_centrality();

            return non_central_t_cdf(v, l, t, true);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            double v = degrees_of_freedom();
            double l = non_centrality();
            return non_central_t_quantile(v, l, p, 1-p);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            double v = degrees_of_freedom();
            double l = non_centrality();
            return non_central_t_quantile(v, l, 1-q, q);
        }

        public override double mean()
        {
            return mean_imp(m_df, m_lambda);
        }

        public override double variance()
        {
            return variance_imp(m_df, m_lambda);
        }

        public override double mode()
        {
            double v = degrees_of_freedom();
            double l = non_centrality();
            double m = v < 3 ? 0 : mean_imp(v, l);
            double var = v < 4 ? 1 : variance_imp(v, l);

            return generic_find_mode(m, Math.Sqrt(var));
        }

        public override double skewness()
        {
            return skewness_imp(m_df, m_lambda);
        }

        public override double kurtosis_excess()
        {
            return kurtosis_excess_imp(m_df, m_lambda);
        }

        #region implementation
        
         double non_central_t2_p(double n, double delta, double x, double y, double init_val)
         {
            //
            // Variables come first:
            //
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)d2;
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative(k+1, d2)  * XMath.gamma_delta_ratio(k + 1, 0.5) * delta / XMath.root_two;
            if(pois == 0) return init_val;
            // Recurance term:
            double xterm = 0;
            // Starting beta term:
            double beta = x < y ? XMath.ibeta_imp((k + 1), (n / 2), x, false, true, ref xterm) : XMath.ibeta_imp((n / 2), (k + 1), y, true, true, ref xterm);

            xterm *= y / (n / 2 + k);
            double poisf = pois, betaf = beta, xtermf = xterm;
            double sum = init_val;
            if((xterm == 0) && (beta == 0)) return init_val;

            //
            // Backwards recursion first, this is the stable
            // direction for recursion:
            //
            int count = 0;
            for(int i = k; i >= 0; --i)
            {
               double term = beta * pois;
               sum += term;
               if(Math.Abs(term/sum) < errtol) break;
               pois *= (i + 0.5) / d2;
               beta += xterm;
               xterm *= (i) / (x * (n / 2 + i - 1));
               ++count;
            }
            for(int i = k + 1; ; ++i)
            {
               poisf *= d2 / (i + 0.5);
               xtermf *= (x * (n / 2 + i - 1)) / (i);
               betaf -= xtermf;
               double term = poisf * betaf;
               sum += term;
               if(Math.Abs(term/sum) < errtol) break;
               ++count;
               if(count > max_iter) throw new Exception("Non-central T distribution: CDF: series did not converge");
            }
            return sum;
         }

         double non_central_t2_q(double n, double delta, double x, double y, double init_val)
         {
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)(d2);
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative(k+1, d2) * XMath.gamma_delta_ratio(k + 1, 0.5) * delta / XMath.root_two;
            if(pois == 0) return init_val;
            // Recurance term:
            double xterm = 0;
            // Starting beta term:
            double beta = x < y ? XMath.ibeta_imp((k + 1), (n / 2), x, true, true, ref xterm) : XMath.ibeta_imp((n / 2), (k + 1), y, false, true, ref xterm);

            xterm *= y / (n / 2 + k);
            double poisf = pois, betaf = beta, xtermf = xterm;
            double sum = init_val;
            if((xterm == 0) && (beta == 0)) return init_val;

            //
            // Forward recursion first, this is the stable direction:
            //
            int count = 0;
            for(int i = k + 1; ; ++i)
            {
               poisf *= d2 / (i + 0.5);
               xtermf *= (x * (n / 2 + i - 1)) / (i);
               betaf += xtermf;

               double term = poisf * betaf;
               sum += term;
               if(Math.Abs(term/sum) < errtol) break;
               if(count > max_iter) throw new Exception("Non-central T distribution: CDF: series did not converge");
               ++count;
            }
            //
            // Backwards recursion:
            //
            for(int i = k; i >= 0; --i)
            {
               double term = beta * pois;
               sum += term;
               if(Math.Abs(term/sum) < errtol) break;
               pois *= (i + 0.5) / d2;
               beta -= xterm;
               xterm *= (i) / (x * (n / 2 + i - 1));
               ++count;
               if(count > max_iter)  throw new Exception("Non-central T distribution: CDF: series did not converge");
            }
            return sum;
         }

         double non_central_t_cdf(double n, double delta, double t, bool invert)
         {
            //
            // For t < 0 we have to use reflect:
            //
            if(t < 0)
            {
               t = -t;
               delta = -delta;
               invert = !invert;
            }
            //
            // x and y are the corresponding random
            // variables for the noncentral beta distribution,
            // with y = 1 - x:
            //
            double x = t * t / (n + t * t);
            double y = n / (n + t * t);
            double d2 = delta * delta;
            double a = 0.5;
            double b = n / 2;
            double c = a + b + d2 / 2;
            //
            // Crossover point for calculating p or q is the same
            // as for the noncentral beta:
            //
            double cross = 1 - (b / c) * (1 + d2 / (2 * c * c));
            double result;
            if(x < cross)
            {
               //
               // Calculate p:
               //
               if(x != 0)
               {
                  result = non_central_beta_distribution.non_central_beta_p(a, b, d2, x, y, 0);
                  result = non_central_t2_p(n, delta, x, y, result);
                  result /= 2;
               }
               else
                  result = 0;
               result +=  new normal_distribution(0, 1).cdf(-delta);
            }
            else
            {
               //
               // Calculate q:
               //
               invert = !invert;
               if(x != 0)
               {
                  result = non_central_beta_distribution.non_central_beta_q(a, b, d2, x, y, 0);
                  result = non_central_t2_q(n, delta, x, y, result);
                  result /= 2;
               }
               else
                  result = new normal_distribution(0, 1).cdfc(-delta);
            }
            if(invert)
               result = 1 - result;
            return result;
         }

         double non_central_t_quantile(double v, double delta, double p, double q)
         {
            double guess = 0;
            if(v > 3)
            {
               double mean = delta * Math.Sqrt(v / 2) * XMath.gamma_delta_ratio((v - 1) * 0.5, (0.5));
               double var = ((delta * delta + 1) * v) / (v - 2) - mean * mean;
               if(p < q)
                  guess = new normal_distribution(mean, var).quantile(p);
               else
                  guess = new normal_distribution(mean, var).quantilec(q);
            }
            //
            // We *must* get the sign of the initial guess correct, 
            // or our root-finder will fail, so double check it now:
            //
            double pzero = non_central_t_cdf(v, delta, 0, !(p<q));
            int s;
            if(p < q)
               s = Math.Sign(p - pzero);
            else
               s = Math.Sign(pzero - q);
            if(s != Math.Sign(guess))
            {
               guess = s;
            }

            double result = generic_quantile(p < q ? p : q, guess, (p >= q));
            return result;
         }

         double non_central_t2_pdf(double n, double delta, double x, double y, double init_val)
         {
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double d2 = delta * delta / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)(d2);
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative((k+1), d2) * XMath.gamma_delta_ratio((k + 1), (0.5)) * delta / XMath.root_two;
            // Starting beta term:
            double xterm = x < y ? XMath.ibeta_derivative((k + 1), n / 2, x) : XMath.ibeta_derivative(n / 2, (k + 1), y);
            double poisf = pois, xtermf = xterm;
            double sum = init_val;
            if((pois == 0) || (xterm == 0)) return init_val;

            //
            // Backwards recursion first, this is the stable
            // direction for recursion:
            //
            int count = 0;
            for(int i = k; i >= 0; --i)
            {
               double term = xterm * pois;
               sum += term;
               if((Math.Abs(term/sum) < errtol) || (term == 0))
                  break;
               pois *= (i + 0.5) / d2;
               xterm *= (i) / (x * (n / 2 + i));
               ++count;
               if(count > max_iter)   throw new Exception("Non-central T distribution: PDF: series did not converge");
            }
            for(int i = k + 1; ; ++i)
            {
               poisf *= d2 / (i + 0.5);
               xtermf *= (x * (n / 2 + i)) / (i);
               double term = poisf * xtermf;
               sum += term;
               if((Math.Abs(term/sum) < errtol) || (term == 0))
                  break;
               ++count;
               if(count > max_iter)   throw new Exception("Non-central T distribution: PDF: series did not converge");
            }
            return sum;
         }

         double non_central_t_pdf(double n, double delta, double t)
         {
            //
            // For t < 0 we have to use the reflection formula:
            //
            if(t < 0)
            {
               t = -t;
               delta = -delta;
            }
            if(t == 0)
            {
               return XMath.gamma_delta_ratio(n / 2 + 0.5, (0.5)) * Math.Sqrt(n / Math.PI) * Math.Exp(-delta * delta / 2) / 2;
            }
            //
            // x and y are the corresponding random
            // variables for the noncentral beta distribution,
            // with y = 1 - x:
            //
            double x = t * t / (n + t * t);
            double y = n / (n + t * t);
            double a = 0.5;
            double b = n / 2;
            double d2 = delta * delta;
            //
            // Calculate pdf:
            //
            double dt = n * t / (n * n + 2 * n * t * t + t * t * t * t);
            double result = non_central_beta_distribution.non_central_beta_pdf(a, b, d2, x, y);
            double tol = XMath.epsilon * result * 500;
            result = non_central_t2_pdf(n, delta, x, y, result);
            if(result <= tol) result = 0;
            result *= dt;
            return result;
         }

         double mean_imp(double v, double delta)
         {
            return delta * Math.Sqrt(v / 2) * XMath.gamma_delta_ratio((v - 1) * 0.5, (0.5));
         }

         double variance_imp(double v, double delta)
         {
            double result = ((delta * delta + 1) * v) / (v - 2);
            double m = mean_imp(v, delta);
            result -= m * m;
            return result;
         }

         double skewness_imp(double v, double delta)
         {
            double mean = mean_imp(v, delta);
            double l2 = delta * delta;
            double var = ((l2 + 1) * v) / (v - 2) - mean * mean;
            double result = -2 * var;
            result += v * (l2 + 2 * v - 3) / ((v - 3) * (v - 2));
            result *= mean;
            result /= Math.Pow(var, (1.5));
            return result;
         }

         double kurtosis_excess_imp(double v, double delta)
         {
            double mean = mean_imp(v, delta);
            double l2 = delta * delta;
            double var = ((l2 + 1) * v) / (v - 2) - mean * mean;
            double result = -3 * var;
            result += v * (l2 * (v + 1) + 3 * (3 * v - 5)) / ((v - 3) * (v - 2));
            result *= -mean * mean;
            result += v * v * (l2 * l2 + 6 * l2 + 3) / ((v - 4) * (v - 2));
            result /= var * var;
            return result;
         }

        #endregion

    }
}
