using System;
using System.Collections.Generic;

using System.Text;

namespace CSBoost.Distributions
{
    public class non_central_chisquared_distribution : distribution
    {
        double m_df, m_lambda, m_mode;

        public non_central_chisquared_distribution(double df, double lambda)
        {
            m_df = df;
            m_lambda = lambda;
            check_parameters();
            m_mode = double.NaN;
        }

        public override void check_parameters()
        {
            if (m_df <= 0 || double.IsInfinity(m_df)) throw new ArgumentException(string.Format("Non-Central Chi-Squared Distribution: degrees of freedom argument must be a finite number > 0 (got {0:G}).", m_df));
            if (m_lambda < 0 || double.IsInfinity(m_lambda)) throw new ArgumentException(string.Format("Non-Central Chi-Squared Distribution: non-centrality argument must be a finite number >= 0 (got {0:G}).", m_lambda));
        }

        public override bool discrete() { return false; }

        public double degrees_of_freedom() { return m_df; }

        public double non_centrality() { return m_lambda; }

        public override bool LHS() { return m_df > 2; }

        public override bool DCL()
        {
            return false;
        }  

        public override bool tail_left() { return false; }

        public override bool unimodal()
        {
            return m_df > 2;
        }

        public override XMath.pair<double> range()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override XMath.pair<double> support()
        {
            return new XMath.pair<double>(0, double.MaxValue);
        }

        public override double pdf(double x)
        {
            if (m_lambda == 0) return new chisquared_distribution(m_df).pdf(x);
            base.pdf(x);
            double r;

            if (x == 0) return 0;
            if (m_lambda > 50) r = pdf_imp(x, m_df, m_lambda);
            else
            {
                r = Math.Log(x / m_lambda) * (m_df / 4 - 0.5) - (x + m_lambda) / 2;
                if (Math.Abs(r) >= XMath.log_max_value / 4) r = pdf_imp(x, m_df, m_lambda);
                else
                {
                    r = Math.Exp(r);
                    r = 0.5 * r * XMath.cyl_bessel_i(m_df / 2 - 1, Math.Sqrt(m_lambda * x));
                }
            }
            return r;
        }

        public override double min_pdf()
        {
            return 0;
        }

        public override double max_pdf()
        {
            if (m_lambda == 0) return new chisquared_distribution(m_df).max_pdf();
            return pdf(mode());
        }

        public override double pdf_inv(double p, bool RHS)
        {
            if (mode() == 0 && !RHS) throw new Exception("Non-Central Chi-Squared Distribution: with the parameters supplied this distribution has no LHS - so no LHS inverse pdf is available.");
            base.pdf_inv(p, RHS);
            if (p == 0)
            {
                if (RHS) return double.MaxValue;
                if (pdf(0) == 0) return 0;
                throw new Exception(string.Format("Non-Central Chi-Squared Distribution: pdf_inv: there is no solution for pdf(x) = {0:G} on the LHS.", p));
            }
            double ubound;
            double lbound;
            if (RHS)
            {
                lbound = mode();
                ubound = support().v2;
            }
            else
            {
                lbound = 0;
                ubound = mode();
            }
            return find_pdf_inv(p, lbound, ubound, true);
        }       

        public override double cdf(double x)
        {
            base.cdf(x);
            return cdf_imp(x, m_df, m_lambda, false);
        }

        public override double cdfc(double x)
        {
            base.cdfc(x);
            return cdf_imp(x, m_df, m_lambda, true);
        }

        public override double quantile(double p)
        {
            base.quantile(p);
            return nccs_quantile(p, false);
        }

        public override double quantilec(double q)
        {
            base.quantilec(q);
            return nccs_quantile(q, true);
        }

        public override double mean()
        {
            return m_df + m_lambda;
        }

        public override double variance()
        {
            return 2 * (2 * m_lambda + m_df);
        }

        public override double mode()        
        {
            if (m_lambda == 0) return new chisquared_distribution(m_df).mode();
            if (double.IsNaN(m_mode)) m_mode = generic_find_mode(m_df + m_lambda, 0);
            return m_mode;
        }

        //Median supplied by base class

        public override double skewness()
        {
            return Math.Pow(2 / (m_df + 2 * m_lambda), 3.0 / 2) * (m_df + 3 * m_lambda);
        }

        //Kurtosis supplied by base class

        public override double kurtosis_excess()
        {
            return 12 * (m_df + 4 * m_lambda) / ((m_df + 2 * m_lambda) * (m_df + 2 * m_lambda));
        }

        #region implementation functions

        private static double pdf_imp(double x, double n, double lambda)
        {
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double x2 = x / 2;
            double n2 = n / 2;
            double l2 = lambda / 2;
            double sum = 0;
            int k = itrunc(l2);
            double pois = XMath.gamma_p_derivative(k + 1, l2) * XMath.gamma_p_derivative(n2 + k, x2);
            if (pois == 0) return 0;
            double poisb = pois;
            for (int i = k; ; ++i)
            {
                sum += pois;
                if (pois / sum < errtol)
                    break;
                if ((i - k) >= max_iter)
                    throw new Exception(string.Format("Non-Central Chi-Squared Distribution: pdf: series did not converge, closest value was {0:G}", sum));
                pois *= l2 * x2 / ((i + 1) * (n2 + i));
            }
            for (int i = k - 1; i >= 0; --i)
            {
                poisb *= (i + 1) * (n2 + i) / (l2 * x2);
                sum += poisb;
                if (poisb / sum < errtol)
                    break;
            }
            return sum / 2;
        }

        private static int itrunc(double x)
        {
            if (x > int.MaxValue || x < int.MinValue) throw new OverflowException();
            return (int)x;
        }

        private static long ltrunc(double x)
        {
            if (x > long.MaxValue || x < long.MinValue) throw new OverflowException();
            return (long)x;
        }

        double cdf_imp(double x, double k, double l, bool invert)
        {
            double result;
            if (l == 0) result = new chisquared_distribution(k).cdf(x);
            else if (x > k + l)
            {
                // Complement is the smaller of the two:
                result = non_central_chi_square_q(x, k, l, invert ? 0 : -1);
                invert = !invert;
            }
            else if (l < 200) result = non_central_chi_square_p_ding(x, k, l, invert ? -1 : 0);
            else result = non_central_chi_square_p(x, k, l, invert ? -1 : 0);
            if (invert) result = -result;
            return result;
        }

        double nccs_quantile(double p, bool comp)
        {

            double k = degrees_of_freedom();
            double l = non_centrality();

            double b = (l * l) / (k + 3 * l);
            double c = (k + 3 * l) / (k + 2 * l);
            double ff = (k + 2 * l) / (c * c);
            double guess;
            if (comp)
                guess = b + c * new chisquared_distribution(ff).cdfc(p);
            else
                guess = b + c * new chisquared_distribution(ff).cdf(p);

            if (guess < 0) guess = double.Epsilon;

            double result = generic_quantile(p, guess, comp);

            return result;
        }

        double non_central_chi_square_q(double x, double f, double theta, double init_sum)
        {
            if(x == 0) return 1;

            double lambda = theta / 2;
            double del = f / 2;
            double y = x / 2;
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double sum = init_sum;
            int k = (int)Math.Round(lambda);
            double poisf = XMath.gamma_p_derivative(1 + k, lambda);
            double poisb = poisf * k / lambda;
            double gamf = XMath.gamma_q(del + k, y);
            double xtermf = XMath.gamma_p_derivative(del + 1 + k, y);
            double xtermb = xtermf * (del + k) / y;
            double gamb = gamf - xtermb;

            int i;
            for(i = k; (i-k) < max_iter; ++i)
            {
               double term = poisf * gamf;
               sum += term;
               poisf *= lambda / (i + 1);
               gamf += xtermf;
               xtermf *= y / (del + i + 1);
               if(((sum == 0) || (Math.Abs(term / sum) < errtol)) && (term >= poisf * gamf)) break;
            }
            if((i-k) >= max_iter) throw new Exception(string.Format("Non-central chi-squared distribution: cdf: series did not converge, closest value was {0:G}", sum));
            for(i = k - 1; i >= 0; --i)
            {
               double term = poisb * gamb;
               sum += term;
               poisb *= i / lambda;
               xtermb *= (del + i) / y;
               gamb -= xtermb;
               if((sum == 0) || (Math.Abs(term / sum) < errtol)) break;
            }

            return sum;
        }

        double non_central_chi_square_p_ding(double x, double f, double theta, double init_sum)
        {
            if(x == 0) return 0;
            double tk = XMath.gamma_p_derivative(f/2 + 1, x/2);
            double lambda = theta / 2;
            double vk = Math.Exp(-lambda);
            double uk = vk;
            double sum = init_sum + tk * vk;
            if(sum == 0) return sum;

            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;

            int i;
            double lterm = 0, term = 0;
            for(i = 1; i < max_iter; ++i)
            {
               tk = tk * x / (f + 2 * i);
               uk = uk * lambda / i;
               vk = vk + uk;
               lterm = term;
               term = vk * tk;
               sum += term;
               if((Math.Abs(term / sum) < errtol) && (term <= lterm)) break;
            }
            //Error check:
            if(i >= max_iter) throw new Exception(string.Format("Non-central chi-squared distribution: cdf: series did not converge, closest value was {0:G}", sum));
            return sum;
        }

        double non_central_chi_square_p(double y, double n, double lambda, double init_sum)
        {
            if(y == 0) return 0;
            int max_iter = XMath.max_series_iterations;
            double errtol = XMath.epsilon;
            double errorf = 0, errorb = 0;

            double x = y / 2;
            double del = lambda / 2;
            int k = (int)Math.Round(del);
            double a = n / 2 + k;
            double gamkf = XMath.gamma_p(a, x);

            if(lambda == 0) return gamkf;
            double gamkb = gamkf;
            double poiskf = XMath.gamma_p_derivative(k+1, del);
            double poiskb = poiskf;
            double xtermf = XMath.gamma_p_derivative(a, x);
            double xtermb = xtermf * x / a;
            double sum = init_sum + poiskf * gamkf;
            if(sum == 0) return sum;

            int i = 1;
            while(i <= k)
            {
               xtermb *= (a - i + 1) / x;
               gamkb += xtermb;
               poiskb = poiskb * (k - i + 1) / del;
               errorf = errorb;
               errorb = gamkb * poiskb;
               sum += errorb;
               if((Math.Abs(errorb / sum) < errtol) && (errorb <= errorf)) break;
               ++i;
            }

            i = 1;
            do
            {
               xtermf = xtermf * x / (a + i - 1);
               gamkf = gamkf - xtermf;
               poiskf = poiskf * del / (k + i);
               errorf = poiskf * gamkf;
               sum += errorf;
               ++i;
            }
            while(Math.Abs(errorf / sum) > errtol && i < max_iter);

            if(i >= max_iter) throw new Exception(string.Format("Non-central chi-squared distribution: cdf: series did not converge, closest value was {0:G}", sum));

            return sum;
        }

        #endregion
    }
}
