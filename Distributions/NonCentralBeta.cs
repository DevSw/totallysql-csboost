using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class non_central_beta_distribution : distribution
    {
        double m_alpha, m_beta, m_lambda, m_mode, m_antimode;
        beta_distribution betadist;


        public non_central_beta_distribution(double alpha, double beta, double lambda)
        {
            m_alpha = alpha;
            m_beta = beta;
            m_lambda = lambda;
            m_antimode = m_mode = double.NaN;
            check_parameters();
            if (lambda == 0) betadist = new beta_distribution(alpha, beta);
        }

        public override void check_parameters()
        {
            if (m_alpha <= 0 || double.IsInfinity(m_alpha)) throw new ArgumentException(string.Format("Alpha argument must be a finite number > 0 (got {0:G}).", m_alpha));
            if (m_beta <= 0 || double.IsInfinity(m_beta)) throw new ArgumentException(string.Format("Beta argument must be a finite number > 0 (got {0:G}).", m_beta));
            if (m_lambda < 0 || double.IsInfinity(m_lambda)) throw new ArgumentException(string.Format("Lambda argument must be a finite number >= 0 (got {0:G}).", m_lambda));
        }

        public override bool discrete() { return false; }

        public override bool inverted()
        {
            return m_alpha < 1 && m_beta < 1;
        }

        public override bool unimodal()
        {
            return (m_alpha >= 1 && m_beta >= 1);
        }

        public double alpha()
        {
            return m_alpha;
        }

        public double beta()
        {
            return m_beta;
        }

        public double non_centrality()
        {
            return m_lambda;
        }

        public override bool LHS() { return !(m_alpha < 1 && m_beta >= 1 || m_alpha == 1 && m_beta > 1); }

        public override bool RHS() { return !(m_beta < 1 && m_alpha >= 1 || m_beta == 1 && m_alpha > 1); }

        public override bool DCL() { return m_alpha < 1; }

        public override bool DCR() { return m_beta < 1; }

        public override bool tail_right() { return false; }

        public override bool tail_left() { return false; }

        public override bool symmetric() { return m_alpha == m_beta && m_lambda == 0; }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, 1);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, 1);
        }

        public override double pdf(double x)
        {
            if (x == 0 && m_alpha > 1) return 0;
            if (x == 1 && m_beta > 1) return 0;
            if (x == 0 && m_alpha < 1) return double.MaxValue;
            if (x == 1 && m_beta < 1) return double.MaxValue;
            if (m_lambda == 0) return betadist.pdf(x);
            base.pdf(x);
            return nc_beta_pdf(x);
        }

        public override double max_pdf()
        {
            if (m_lambda == 0) return betadist.max_pdf();
            if (m_alpha < 1 || m_beta < 1) return double.MaxValue;
            return pdf(mode());
        }

        public override double min_pdf()
        {
            if (m_alpha > 1 || m_beta > 1) return 0;
            if (m_alpha == 1 && m_beta == 1 && m_lambda == 0) return 1;
            if (m_alpha == 1 && m_beta == 1 && m_lambda > 0) return pdf(0);
            return pdf(antimode());
        }

        public override double pdf_inv(double p, bool rhs)
        {
            if (m_lambda == 0) return betadist.pdf_inv(p, rhs);
            if (rhs && !RHS()) throw new Exception(string.Format("The non-central beta distribution has no RHS when beta {0:G} 1 and alpha {1:G} 1.", m_beta < 1 ? "<" : "==", m_beta < 1 ? ">=" : ">"));
            if (!rhs && !LHS()) throw new Exception(string.Format("The non-central beta distribution has no LHS when alpha {0:G} 1 and beta {1:G} 1.", m_alpha < 1 ? "<" : "==", m_alpha < 1 ? ">=" : ">"));
            base.pdf_inv(p, rhs);
            double result = 0, midpoint;
            if (m_alpha == 1 && m_beta == 1 && m_lambda == 0) throw new Exception("Inverse pdf is not defined for this distribution when alpha = beta = 1 and lambda = 0");
            if (m_alpha < 1 && m_beta < 1)
            {
                midpoint = antimode();
                double pmin = pdf(midpoint);
                if (p < pmin) throw new ArgumentException(string.Format("Minimum pdf value for this distribution is {0:G} (got {1:G}).", pmin, p));
                if (rhs) result = find_pdf_inv(p, midpoint, support().v2, true);
                else result = find_pdf_inv(p, support().v1, midpoint, false);
            }
            else if (m_alpha > 1 && m_beta > 1)
            {
                midpoint = mode();
                if (rhs) result = find_pdf_inv(p, midpoint, support().v2, false);
                else result = find_pdf_inv(p, support().v1, midpoint, true);
            }
            else if (m_alpha < 1 || (m_alpha == 1 && m_beta >= 1)) result = find_pdf_inv(p, support().v1, support().v2, false);
            else result = find_pdf_inv(p, support().v1, support().v2, true);
            return result;
        }          

        public override double cdf(double x)
        {
            if (m_lambda == 0) return betadist.cdf(x);
            base.cdf(x);
            if (x == 0) return 0;
            else if (x == 1) return 1;
            return non_central_beta_cdf(x, 1-x, m_alpha, m_beta, m_lambda, false);
        }

        public override double cdfc(double x)
        {
            if (m_lambda == 0) return betadist.cdfc(x);
            base.cdfc(x);
            if (x == 0) return 1;
            else if (x == 1) return 0;
            return non_central_beta_cdf(x, 1 - x, m_alpha, m_beta, m_lambda, true);
        }

        public override double quantile(double p)
        {
            if (m_lambda == 0) return betadist.quantile(p);
            base.quantile(p);
            if (p == 0) return 0;
            if (p == 1) return 1;
            return nc_beta_quantile(p, false);
        }

        public override double quantilec(double p)
        {
            if (m_lambda == 0) return betadist.quantilec(p);
            base.quantilec(p);
            if (p == 0) return 1;
            if (p == 1) return 0;
            return nc_beta_quantile(p, true);
        }

        public override double mean()
        {
            if (m_lambda == 0) return betadist.mean();
            else throw new NotImplementedException();
        }

        public override double variance()
        {
            if (m_lambda == 0) return betadist.variance();
            else throw new NotImplementedException();
        }

        public override double antimode()
        {
            if (m_lambda == 0) return betadist.antimode();
            if (m_alpha > 1) throw new Exception(string.Format("Non-central beta: antimode undefined for alpha = {0:G}, must be <= 1!", m_alpha));
            if (m_beta > 1) throw new Exception(string.Format("Non-central beta: antimode undefined for beta = {0:G}, must be <= 1!", m_beta));

            if (double.IsNaN(m_antimode))
            {
                double c = m_alpha + m_beta + m_lambda / 2;
                double mean = 1 - (m_beta / c) * (1 + m_lambda / (2 * c * c));
                m_antimode = generic_find_antimode_01(mean);
            }
            return m_antimode;  //Cached for better performance on multiple calls
        }

        public override double mode()
        {
            if (m_lambda == 0) return betadist.mode();

            if (m_alpha < 1) throw new Exception(string.Format("Non-central beta: mode undefined for alpha = {0:G}, must be >= 1!", m_alpha));
            if (m_beta < 1) throw new Exception(string.Format("Non-central beta: mode undefined for beta = {0:G}, must be >= 1!", m_beta));

            if (double.IsNaN(m_mode))
            {
                double c = m_alpha + m_beta + m_lambda / 2;
                double mean = 1 - (m_beta / c) * (1 + m_lambda / (2 * c * c));
                m_mode = generic_find_mode_01(mean);
            }
            return m_mode;  //Cached for better performance on multiple calls
        }

        //median supplied by base class

        public override double skewness()
        {
            if (m_lambda == 0) return betadist.skewness();
            else throw new NotImplementedException();
        }

        public override double kurtosis_excess()
        {
            if (m_lambda == 0) return betadist.kurtosis();
            else throw new NotImplementedException();
        }

        //kurtosis supplied by base class

        public override void setup_bars()
        {
            if (m_lambda == 0) betadist.setup_bars();
            else base.setup_bars();
        }

        public override double random()
        {
            if (m_lambda == 0) return betadist.random();
            return(base.random());
        }

        #region implementation

        internal static double non_central_beta_p(double a, double b, double lam, double x, double y, double init_val)
        {
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double l2 = lam / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)(l2);
            if(k == 0)
               k = 1;
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative(k+1.0, l2);
            if(pois == 0) return init_val;
            // recurance term:
            double xterm = 0;
            // Starting beta term:
            double beta = x < y ? XMath.ibeta_imp(a + k, b, x, false, true, ref xterm) : XMath.ibeta_imp(b, a + k, y, true, true, ref xterm);

            xterm *= y / (a + b + k - 1);
            double poisf = pois, betaf = beta, xtermf = xterm;
            double sum = init_val;

            if((beta == 0) && (xterm == 0)) return init_val;

            //
            // Backwards recursion first, this is the stable
            // direction for recursion:
            //
            double last_term = 0;
            int count = k;
            for(int i = k; i >= 0; --i)
            {
               double term = beta * pois;
               sum += term;
               if(((Math.Abs(term/sum) < errtol) && (last_term >= term)) || (term == 0))
               {
                  count = k - i;
                  break;
               }
               pois *= i / l2;
               beta += xterm;
               xterm *= (a + i - 1) / (x * (a + b + i - 2));
               last_term = term;
            }
            for(int i = k + 1; ; ++i)
            {
               poisf *= l2 / i;
               xtermf *= (x * (a + b + i - 2)) / (a + i - 1);
               betaf -= xtermf;

               double term = poisf * betaf;
               sum += term;
               if((Math.Abs(term/sum) < errtol) || (term == 0))
               {
                  break;
               }
               if(count + i - k > max_iter) throw new Exception("Non-Central Beta Distribution: CDF: Series did not converge");
            }
            return sum;
        }

        internal static double non_central_beta_q(double a, double b, double lam, double x, double y, double init_val)
        {
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double l2 = lam / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)(l2);
            if(k == 0)
               k = 1;
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative(k+1, l2);
            if(pois == 0) return init_val;
            // recurance term:
            double xterm = 0;
            // Starting beta term:
            double beta = x < y ? XMath.ibeta_imp(a + k, b, x, true, true, ref xterm) : XMath.ibeta_imp(b, a + k, y, false, true, ref xterm);

            xterm *= y / (a + b + k - 1);
            double poisf = pois, betaf = beta, xtermf = xterm;
            double sum = init_val;
            if((beta == 0) && (xterm == 0)) return init_val;
            //
            // Forwards recursion first, this is the stable
            // direction for recursion, and the location
            // of the bulk of the sum:
            //
            double last_term = 0;
            int count = 0;
            for(int i = k + 1; ; ++i)
            {
               poisf *= l2 / i;
               xtermf *= (x * (a + b + i - 2)) / (a + i - 1);
               betaf += xtermf;

               double term = poisf * betaf;
               sum += term;
               if((Math.Abs(term/sum) < errtol) && (last_term >= term))
               {
                  count = i - k;
                  break;
               }
               if(i - k > max_iter)  throw new Exception("Non-Central Beta Distribution: CDF: Series did not converge");
               last_term = term;
            }
            for(int i = k; i >= 0; --i)
            {
               double term = beta * pois;
               sum += term;
               if(Math.Abs(term/sum) < errtol)
               {
                  break;
               }
               if(count + k - i > max_iter)  throw new Exception("Non-Central Beta Distribution: CDF: Series did not converge");
               pois *= i / l2;
               beta -= xterm;
               xterm *= (a + i - 1) / (x * (a + b + i - 2));
            }
            return sum;
        }

        internal double non_central_beta_cdf(double x, double y, double a, double b, double l, bool invert)
        {
            if(x == 0) return invert ? 1.0 : 0.0;
            if(y == 0) return invert ? 0.0 : 1.0;
            double result;
            double c = a + b + l / 2;
            double cross = 1 - (b / c) * (1 + l / (2 * c * c));
            if(l == 0) result = cdf(x);
            else if(x > cross) 
            {
                result = non_central_beta_q(a, b, l, x, y, invert ? 0: -1);               // Complement is the smaller of the two:
                invert = !invert;
            }
            else result = non_central_beta_p(a, b, l, x, y, invert ? -1: 0);               // Complement is the smaller of the two:
            if(invert) result = -result;
            return result;
        }

        class nc_beta_quantile_functor : XMath.iterand<double, double>
        {
            distribution dist;
            double target;
            bool comp;

            public nc_beta_quantile_functor(distribution d, double t, bool c)
            {
                dist = d;
                target = t;
                comp = c;
            }

            public override double  next(double x)
            {
                return comp ? target - dist.cdfc(x) : dist.cdf(x) - target;
            }
        }

        XMath.pair<double> bracket_and_solve_root_01(XMath.iterand<double, double> f, double guess, double factor, bool rising, XMath.eps_tolerance tol, ref int max_iter)
        {
            double a = guess;
            double b = a;
            double fa = f.next(a);
            double fb = fa;
            //
            // Set up invocation count:
            //
            int count = max_iter - 1;

            if((fa < 0) == (guess < 0 ? !rising : rising))
            {
               //
               // Zero is to the right of b, so walk upwards
               // until we find it:
               //
               while(Math.Sign(fb) == Math.Sign(fa))
               {
                  if(count == 0) throw new Exception("Non-Central Beta Distribution: Unable to bracket root");
                  //
                  // Heuristic: every 20 iterations we double the growth factor in case the
                  // initial guess was *really* bad !
                  //
                  if((max_iter - count) % 20 == 0)
                     factor *= 2;
                  //
                  // Now go ahead and move are guess by "factor",
                  // we do this by reducing 1-guess by factor:
                  //
                  a = b;
                  fa = fb;
                  b = 1 - ((1 - b) / factor);
                  fb = f.next(b);
                  --count;
               }
            }
            else
            {
               //
               // Zero is to the left of a, so walk downwards
               // until we find it:
               //
               while(Math.Sign(fb) == Math.Sign(fa))
               {
                  if(Math.Abs(a) < double.Epsilon)
                  {
                     // Escape route just in case the answer is zero!
                     max_iter -= count;
                     max_iter += 1;
                     return a > 0 ? new XMath.pair<double>(0, a) : new XMath.pair<double>(a, 0);
                  }
                  if(count == 0) throw new Exception("Non-Central Beta Distribution: Unable to bracket root");
                  //
                  // Heuristic: every 20 iterations we double the growth factor in case the
                  // initial guess was *really* bad !
                  //
                  if((max_iter - count) % 20 == 0)
                     factor *= 2;
                  //
                  // Now go ahead and move are guess by "factor":
                  //
                  b = a;
                  fb = fa;
                  a /= factor;
                  fa = f.next(a);
                  --count;
               }
            }
            max_iter -= count;
            max_iter += 1;
            XMath.pair<double> r = XMath.toms748_solve(f, (a < 0 ? b : a), (a < 0 ? a : b), (a < 0 ? fb : fa), (a < 0 ? fa : fb), tol, ref count);
            max_iter += count;
            return r;
        }

        double nc_beta_quantile(double p, bool comp)
        {
            double a = alpha();
            double b = beta();
            double l = non_centrality();
            //
            // Special cases first:
            //
            if(p == 0) return comp ? 1.0 : 0.0;
            if(p == 1) return !comp ? 1.0 : 0.0;

            double c = a + b + l / 2;
            double mean = 1 - (b / c) * (1 + l / (2 * c * c));
            double guess = mean;
            nc_beta_quantile_functor f = new nc_beta_quantile_functor(this, p, comp);
            XMath.eps_tolerance tol = new XMath.eps_tolerance(53);
            int max_iter = XMath.max_root_iterations;

            XMath.pair<double> ir = bracket_and_solve_root_01(f, guess, 2.5, true, tol, ref max_iter);
            double result = ir.v1 + (ir.v2 - ir.v1) / 2;

            if(max_iter >= XMath.max_root_iterations) throw new Exception("Non-Central Beta Distribution: Unable to locate solution in a reasonable time");
            return result;
         }

         internal static double non_central_beta_pdf(double a, double b, double lam, double x, double y)
         {
            //
            // Variables come first:
            //
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double l2 = lam / 2;
            //
            // k is the starting point for iteration, and is the
            // maximum of the poisson weighting term:
            //
            int k = (int)(l2);
            // Starting Poisson weight:
            double pois = XMath.gamma_p_derivative(k+1, l2);
            // Starting beta term:
            double beta = x < y ? XMath.ibeta_derivative(a + k, b, x) : XMath.ibeta_derivative(b, a + k, y);
            double sum = 0;
            double poisf = pois;
            double betaf = beta;

            //
            // Stable backwards recursion first:
            //
            int count = k;
            for(int i = k; i >= 0; --i)
            {
               double term = beta * pois;
               sum += term;
               if((Math.Abs(term/sum) < errtol) || (term == 0))
               {
                  count = k - i;
                  break;
               }
               pois *= i / l2;
               beta *= (a + i - 1) / (x * (a + i + b - 1));
            }
            for(int i = k + 1; ; ++i)
            {
               poisf *= l2 / i;
               betaf *= x * (a + b + i - 1) / (a + i - 1);

               double term = poisf * betaf;
               sum += term;
               if((Math.Abs(term/sum) < errtol) || (term == 0))
               {
                  break;
               }
               if(count + i - k > max_iter)  throw new Exception("Non-Central Beta Distribution: Unable to locate solution in a reasonable time");
            }
            return sum;
         }

         double nc_beta_pdf(double x)
         {
            double a = alpha();
            double b = beta();
            double l = non_centrality();
            return non_central_beta_pdf(a, b, l, x, 1 - x);
         }


        #endregion

    }
}
